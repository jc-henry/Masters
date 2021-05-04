library('phangorn')
library('strex')
library('tools')
library('tidyverse')

wd <- getwd()

# set file path of input directory
files <- '../mltree_support/'

dir.create('reject_bs_trees')
dir.create('good_trees')
dir.create('good_nolen_trees')
rej_bs <- paste(wd,'reject_bs_trees/', sep='/')
gd_tre <- paste(wd, 'good_trees/', sep='/')
gd_nolen <- paste(wd,'good_nolen_trees/', sep='/')

# get the list of supported ML trees from RAxML-NG 
nwk_files <- list.files(path=files, pattern="*.raxml.support", full.names=TRUE, recursive=FALSE)

# initialise a list for adding good trees
all_good_trees_list <- c()

for(i in 1:length(nwk_files)) {
  input <-nwk_files[[i]]
  tre <- read.tree(input)
  tre_nam <- str_before_first(basename(input),'\\.') # get the cluster name
  bt <- as.integer(tre$node.label[2]) # get the central split bootstrap value
  if(bt < 70) {
    write.tree(tre, file = paste(rej_bs, tre_nam, '.fa.nwk', sep = '')) # write the tree to the "reject_bs_trees" dir
   }else{
    write.tree(tre, file = paste(gd_tre, tre_nam, '.fa.nwk', sep = '')) # write the tree to the "good_trees" dir
    el <- length(tre$edge.length) # get the number of edges
    tre$edge.length <- rep(1.0,times = el) # force the edge lengths to one for the midpointing
    mp_tre <- midpoint(tre, node.labels = "support") # need to midpoint root to interrogate the newick file
    mp_tre$edge.length <- NULL # remove edge lengths
    mp_tre$node.label <- NULL # remove bootstrap values
    write.tree(mp_tre, file = paste(gd_nolen, tre_nam, '.fa.nwk', sep = ''))
    all_good_trees_list <- append(all_good_trees_list, paste(tre_nam, '.nwk', sep='')) 
  }}

write.table(all_good_trees_list, 'all_accepted_trees.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')


pats = c('_1|_2|_3|\\(|\\)') # patterns to remove with str_remove_all() later
inbrac <- ')' # set the close bracket for comparison later
outbrac <- '(' # set the open bracket for comparison later

# opbrac <- '(' # set the open bracket for extraction later
# clbrac <- ')' # set the close bracket for extraction later

dir.create('cluster_reject_topology_trees')
dir.create('cluster_reject_taxonomy_trees')
dir.create('cluster_porters_crithmifolius_trees')
dir.create('cluster_porters_hutt_trees')
dir.create('cluster_porters_lobulatus_trees')
clus_rejtop <- paste(wd,'cluster_reject_topology_trees/', sep='/')
clus_rejtax <- paste(wd,'cluster_reject_taxonomy_trees/', sep='/')
clus_por_cri <- paste(wd, 'cluster_porters_crithmifolius_trees/', sep='/')
clus_por_hut <- paste(wd,'cluster_porters_hutt_trees/', sep='/')
clus_por_lob <- paste(wd,'cluster_porters_lobulatus_trees/', sep='/')

porters_monroi_crithmifolius_list <- c() # initialise empty lists to append to in the loop
porters_monroi_hutt_monroi_list <- c()
porters_monroi_lobulatus_list <- c()
reject_tree_list <- c()

nk_files <- list.files(path=gd_nolen, pattern="*.nwk", full.names=TRUE, recursive=FALSE) # get the list of newick files

for(i in 1:length(nk_files)) {
  input = nk_files[[i]]
  in_tre <- read.tree(input)
  in_file <- read_file(input) # read in as a character string
  no_bs <- gsub('[[:digit:]]+', '', in_file) # need to remove the bootstrap values
  filpath <- file_path_sans_ext(input) # get filepath (without .nwk extension)
  filnam <- basename(filpath) #extract file/cluster name
  first_half <- str_before_nth(in_file, ',', 6) # extract string before the 6th comma
  second_half <- str_after_nth(in_file, ',', 6)
  last_char <- str_sub(first_half, -1) # extract last character
  first_char <- str_sub(second_half, 0, 1)
  left_brac_cnt <- str_count(first_half, '\\)')
  right_brac_cnt <- str_count(second_half, '\\(')
  if(last_char != inbrac | first_char != outbrac) {
    reject_tree_list <- append(reject_tree_list, filnam) # if the 6th comma isn't placed "),(" append to reject object
    write.tree(in_tre, file = paste(clus_rejtop, filnam, '.newick', sep = ''))
  }else {
    if(left_brac_cnt != 5 | right_brac_cnt != 5) {
      reject_tree_list <- append(reject_tree_list, filnam) # another catch for polytomies
      write.tree(in_tre, file = paste(clus_rejtop, filnam, '.newick', sep = ''))
   }else {
    clean_half <- str_remove_all(first_half, pats) # remove all the unneeded characters
    porters_monroi_cnt <- str_count(clean_half, 'porters') # count the occurrences in the string
    crithmifolius_cnt <- str_count(clean_half, 'crithmifolius')
    hutt_monroi_cnt <- str_count(clean_half, 'hutt')
    lobulatus_cnt <- str_count(clean_half, 'lobulatus')
    if(porters_monroi_cnt+crithmifolius_cnt+hutt_monroi_cnt+lobulatus_cnt !=6) {
      reject_tree_list <- append(reject_tree_list, filnam)
      write.tree(in_tre, file = paste(clus_rejtax, filnam, '.newick', sep = ''))
    }else {
      if(porters_monroi_cnt+crithmifolius_cnt == 6 | hutt_monroi_cnt + lobulatus_cnt == 6) {
        porters_monroi_crithmifolius_list <- append(porters_monroi_crithmifolius_list, filnam)
        write.tree(in_tre, file = paste(clus_por_cri, filnam, '.newick', sep = ''))
      }else {
        if(porters_monroi_cnt+hutt_monroi_cnt == 6 | crithmifolius_cnt + lobulatus_cnt == 6) {
          porters_monroi_hutt_monroi_list <- append(porters_monroi_hutt_monroi_list, filnam)
          write.tree(in_tre, file = paste(clus_por_hut, filnam, '.newick', sep = ''))
        }else {
          if(porters_monroi_cnt+lobulatus_cnt == 6 | crithmifolius_cnt+hutt_monroi_cnt == 6) {
            porters_monroi_lobulatus_list <- append(porters_monroi_lobulatus_list, filnam)
            write.tree(in_tre, file = paste(clus_por_lob, filnam, '.newick', sep = ''))
          }else {
            reject_tree_list <- append(reject_tree_list, filnam)
            write.tree(in_tre, file = paste(clus_rejtax, filnam, '.newick', sep = ''))
          }}}}}}}


write.table(porters_monroi_crithmifolius_list, 'cluster_porters_monroi_crithmifolius.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
write.table(porters_monroi_hutt_monroi_list, 'cluster_porters_monroi_hutt_monroi.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
write.table(porters_monroi_lobulatus_list, 'cluster_porters_monroi_lobulatus.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
write.table(reject_tree_list, 'cluster_reject_trees.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')





