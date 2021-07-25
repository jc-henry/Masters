library('tidyverse')

# get the Trinity gene ID - Araport mapping
map <- read.table('path/to/lobulatus_blast.outfmt6', header = FALSE, sep = '\t')

# get the list of background Trinity genes
bk <- read.table('path/to/background_genes.txt', header = F, sep = '\t')

# make the background with Arabidopsis IDs
bk_df <- left_join(bk, map, by='V1') %>%
  drop_na('V2')                             # drop any gene with no blast match originally

write.table(bk_df, file = 'background_genes_araport.txt', col.names = FALSE, row.names = FALSE, quote = FALSE,  sep = '\t')

# get the Trinity gene ID for the hutt-crith cluster 
hutcri_id <- read.table('../geneIDonly_cluster_hutt_crithmifolius_vs_porters_lobulatus.txt', header = FALSE, sep = '\t')

# Join together and drop rows with NA in V2 (no blast match originally)
hutcri_df <- left_join(hutcri_id, map, by='V1') %>%
  drop_na('V2')

write.table(hutcri_df, file = 'cluster_hutt_crithmifolius_arabidopsis_annotated_geneID.txt', col.names = FALSE, row.names = FALSE, quote = FALSE,  sep = '\t')

# get the Trinity gene ID for the hutt-lob cluster 
hutlob_id <- read.table('../geneIDonly_cluster_hutt_lobulatus_vs_porters_crithmifolius', header = FALSE, sep = '\t')

# Join together and drop rows with NA in V2 (no blast match originally)
hutlob_df <- left_join(hutlob_id, map, by='V1') %>%
  drop_na('V2')

write.table(hutlob_df, file = 'cluster_hutt_lobulatus_arabidopsis_annotated_geneID.txt', col.names = FALSE, row.names = FALSE, quote = FALSE,  sep = '\t')

# get the Trinity gene ID for the hutt-porters cluster 
hutpor_id <- read.table('../geneIDonly_cluster_porters_hutt_vs_lobulatus_crithmifolius.txt', header = FALSE, sep = '\t')

# Join together and drop rows with NA in V2 (no blast match originally)
hutpor_df <- left_join(hutpor_id, map, by='V1') %>%
  drop_na('V2')

write.table(hutpor_df, file = 'cluster_porters_hutt_arabidopsis_annotated_geneID.txt', col.names = FALSE, row.names = FALSE, quote = FALSE,  sep = '\t')
