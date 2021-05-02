### blast search with max_target_seqs set to 5 then filter by lowest E-value

library(readxl)
library(tidyverse)
library(strex)

bla <- read.csv('blastp_output.outfmt6', header = FALSE, sep = '\t')

bla_fin <- bla %>%
  group_by(V1) %>%                  # group by gene ID
  slice_min(V11) %>%                # keep only the ID lowest e-value
  distinct(V1, .keep_all = TRUE)    # retain only unique genes 

write.table(bla_fin, file = 'filtered.outfmt6', col.names = FALSE, row.names = FALSE, quote = FALSE,  sep = '\t')



### generate list of Trinity gene IDs ("_i" and ".p1" removed) with Arabidopsis annotations

bla_fin_fil <- bla_fin %>%
  select(V1,V2)                     # select the contig ID and Arabidopsis ID columns

trin <- bla_fin_fil$V1              # extract the TransDecoder ID
trin_clean <- trin %>%              
  str_before_first('\\.') %>%       # remove the ".p..."
  str_before_last('\\_') %>%        # remove the "_i..."
  as.data.frame()

ara <- as.data.frame(bla_fin_fil$V2)# extract the Arabidopsis IDs

fin_clean <- cbind(trin_clean, ara) # bind back together

write.table(fin_clean, file = 'taxon_arabidopsis.txt', col.names = FALSE, row.names = FALSE, quote = FALSE,  sep = '\t')
