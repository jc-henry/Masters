library('tidyverse')
library('strex')
library('tools')

wd <- getwd()

# read in the SCO TPM values
df_tpm <- read.table('sc_orthologue_TPM.tsv', header = TRUE, sep = '\t')

pats = c('_1|_2|_3|\\(|\\)') # patterns to remove with str_remove_all() later
inbrac <- ')' # set the close bracket for comparison later
outbrac <- '(' # set the open bracket for comparison later

porters_monroi_crithmifolius_list <- c() # intialise empty lists to append to in the loop
porters_monroi_hutt_monroi_list <- c()
porters_monroi_lobulatus_list <- c()
reject_list <- c()

for(i in 1:nrow(df_tpm)) {
  row <- df_tpm[i,]
  row_no <- remove_rownames(row)
  c2r <- column_to_rownames(row_no, var = "Gene") # seems necessary for transposing
  tran <- as.data.frame(t(c2r)) # transpose the df
  col_nm <- colnames(tran) # extract the colname/trinity ID
  colnames(tran) <- "x" # set the colname to "x"
  rw_nm <- rownames_to_column(tran, var='Gene') # make the rownames a column
  ord <- rw_nm[order(rw_nm$x),] # order the df (decreasing TPM)
  top <- ord[1:6,] # extract the top (lowest 6)
  bot <- ord[7:12,]
  last_top <- top$x[6] # testing that the 6th and 7th values aren't the same
  first_bot <- bot$x[1]
  if(last_top == first_bot) {
    reject_list <- append(reject_list, col_nm)
  }else{
    vec <- top$Gene # into a vector
    clean <- str_remove_all(vec, pats) # remove all the unneeded characters
    porters_monroi_cnt <- sum(str_count(clean, 'porters_monroi')) # count the occurrences in the string
    crithmifolius_cnt <- sum(str_count(clean, 'crithmifolius'))
    hutt_monroi_cnt <- sum(str_count(clean, 'hutt_monroi'))
    lobulatus_cnt <- sum(str_count(clean, 'lobulatus'))
    if(porters_monroi_cnt+crithmifolius_cnt == 6 | hutt_monroi_cnt + lobulatus_cnt == 6) {
      porters_monroi_crithmifolius_list <- append(porters_monroi_crithmifolius_list, col_nm)
    }else {
      if(porters_monroi_cnt+hutt_monroi_cnt == 6 | crithmifolius_cnt + lobulatus_cnt == 6) {
        porters_monroi_hutt_monroi_list <- append(porters_monroi_hutt_monroi_list, col_nm)
      }else {
        if(porters_monroi_cnt+lobulatus_cnt == 6 | crithmifolius_cnt+hutt_monroi_cnt == 6) {
          porters_monroi_lobulatus_list <- append(porters_monroi_lobulatus_list, col_nm)
        }else {
          reject_list <- append(reject_list, col_nm)
        }}}}}

write.table(porters_monroi_crithmifolius_list, 'cluster_porters_monroi_crithmifolius.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
write.table(porters_monroi_hutt_monroi_list, 'cluster_porters_monroi_hutt_monroi.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
write.table(porters_monroi_lobulatus_list, 'cluster_porters_monroi_lobulatus.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')
write.table(reject_list, 'cluster_reject.txt', quote = FALSE, row.names = FALSE, col.names = FALSE, sep = '\t')

































