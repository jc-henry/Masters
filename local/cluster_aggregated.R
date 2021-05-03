library('phangorn')
library('strex')
library('tools')
library('tidyverse')

wd <- getwd()

# set file path of input directory
files<-'../mltree_support/'

dir.create('reject_bs_trees')
dir.create('good_trees')
dir.create('good_nolen_trees')
rej_bs <- paste(wd,'reject_bs_trees/', sep='/')
gd_tre <- paste(wd, 'good_trees/', sep='/')        # edge lengths retained for manual examination 
gd_nolen <- paste(wd,'good_nolen_trees/', sep='/') # edge lengths removed for clustering

# get the list of supported ML trees from RAxML-NG 
nwk_files <- list.files(path=files, pattern="*.raxml.support", full.names=TRUE, recursive=FALSE)

# initialise a list for adding good trees
all_good_trees_list <- c()


for(i in 1:length(nwk_files)) {
  input <-nwk_files[[i]]
  tre <- read.tree(input)
  tre_nam <- str_before_first(basename(input),'\\.') # get the cluster name
  bt <- as.integer(tre$node.label[2])
  if(bt < 70) {
    write.tree(tre, file = paste(rej_bs, tre_nam, '.fa.nwk', sep = ''))
   }else{
    write.tree(tre, file = paste(gd_tre, tre_nam, '.fa.nwk', sep = ''))
    tre$edge.length <- NULL # remove edge lengths
    tre$node.label <- NULL # remove bootstrap values
    write.tree(tre, file = paste(gd_nolen, tre_nam, '.fa.nwk', sep = ''))
    all_good_trees_list <- append(all_good_trees_list, paste(tre_nam, '.nwk', sep='')) 
  }}

  

opbrac <- '(' # set the open bracket for extraction later
clbrac <- ')' # set the close bracket for extraction later

porters_monroi_crithmifolius_list <- c() # initialise empty lists to append to in the loop
porters_monroi_hutt_monroi_list <- c()
porters_monroi_lobulatus_list <- c()
reject_tree_list <- c()

nwk_files <- list.files(path=gd_nolen, pattern="*.nwk", full.names=TRUE, recursive=FALSE) # get the list of newick files

for(i in 1:length(nwk_files)) {
  input = nwk_files[[i]]
  in_file <- read_file(input) # read newick file as a character string
  filpath <- file_path_sans_ext(input) # get filepath (without .nwk extension)
  filnam <- basename(filpath) # extract file/cluster name
  first <- str_before_nth(in_file, '\\)', 1) # extract string before the first close bracket (need to escape)
  cluster <- str_after_nth(first, '\\(', 2) # now have two clustered taxa
    porters_monroi_cnt <- str_count(cluster, 'porters') # count the occurrences in the string
    crithmifolius_cnt <- str_count(cluster, 'crithmifolius')
    hutt_monroi_cnt <- str_count(cluster, 'hutt')
    lobulatus_cnt <- str_count(cluster, 'lobulatus')
      if(porters_monroi_cnt+crithmifolius_cnt == 2 | hutt_monroi_cnt + lobulatus_cnt == 2) {
        porters_monroi_crithmifolius_list <- append(porters_monroi_crithmifolius_list, filnam)
      }else {
        if(porters_monroi_cnt+hutt_monroi_cnt == 2 | crithmifolius_cnt + lobulatus_cnt == 2) {
          porters_monroi_hutt_monroi_list <- append(porters_monroi_hutt_monroi_list, filnam)
        }else {
          if(porters_monroi_cnt+lobulatus_cnt == 2 | crithmifolius_cnt+hutt_monroi_cnt == 2) {
            porters_monroi_lobulatus_list <- append(porters_monroi_lobulatus_list, filnam)
          }else {
           reject_tree_list <- append(reject_tree_list, filnam)
          }}}}

write.table(porters_monroi_crithmifolius_list, 'cluster_porters_monroi_crithmifolius.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
write.table(porters_monroi_hutt_monroi_list, 'cluster_porters_monroi_hutt_monroi.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
write.table(porters_monroi_lobulatus_list, 'cluster_porters_monroi_lobulatus.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
write.table(reject_tree_list, 'cluster_reject_trees.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
write.table(all_good_trees_list, 'all_accepted_trees.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')




