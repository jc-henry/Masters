library('tidyverse')
library('strex')

# get the Trinity gene ID - Araport mapping
# use lobulatus_blast.outfmt6 for individual assemblies
map <- read.table('path/to/lobulatus_blast.outfmt6', header = FALSE, sep = '\t')

# get OrthoFinder cluster to lobulatus TransDecoder mapping
# use lobulatus_1 ortholist for individual assemblies
df <- read.table('path/to/ortholist_lobulatus.txt', header = F, sep = '\t')

# make a dataframe of OrthoFinder to Trinity mapping
clus <- str_before_first(df$V1, '_')      # extract the OrthoFinder cluster/gene IDs
gene <- str_before_last(df$V1, "_") %>%   # extract the Trinity gene IDs
  str_after_nth('_', 2)
lob_map <- data.frame(clus, gene)         # make into a dataframe with OrthoFinder and Trinity IDs

# annotate clusters with Araport gene ID
clus_df <- left_join(lob_map, map, by=c('gene' = 'V1')) %>%
  drop_na('V2')
head(clus_df)
dim(clus_df)

write.table(clus_df, file = 'background_genes_araport_expression.txt', col.names = FALSE, row.names = FALSE, quote = FALSE,  sep = '\t')

# get the cluster ID for the hutt-crith cluster 
hutcri_id <- read.table('../cluster_porters_monroi_lobulatus.txt', header = FALSE, sep = '\t')
head(hutcri_id)
dim(hutcri_id)

# Join together and drop rows with NA in V2 (no blast match originally)
hutcri_df <- left_join(hutcri_id, clus_df, by=c('V1' = 'clus')) %>%
  drop_na('V2')
head(hutcri_df)
dim(hutcri_df)

write.table(hutcri_df, file = 'cluster_hutt_crithmifolius_arabidopsis_annotated_geneID.txt', col.names = FALSE, row.names = FALSE, quote = FALSE,  sep = '\t')

# get the cluster ID for the hutt-lob cluster 
hutlob_id <- read.table('../cluster_porters_monroi_crithmifolius.txt', header = FALSE, sep = '\t')
head(hutlob_id)
dim(hutlob_id)

# Join together and drop rows with NA in V2 (no blast match originally)
hutlob_df <- left_join(hutlob_id, clus_df, by=c('V1' = 'clus')) %>%
  drop_na('V2')
head(hutlob_df)
dim(hutlob_df)

write.table(hutlob_df, file = 'cluster_hutt_lobulatus_arabidopsis_annotated_geneID.txt', col.names = FALSE, row.names = FALSE, quote = FALSE,  sep = '\t')

# get the cluster ID for the hutt-porters cluster 
hutpor_id <- read.table('../cluster_porters_monroi_hutt_monroi.txt', header = FALSE, sep = '\t')
head(hutpor_id)
dim(hutpor_id)

# Join together and drop rows with NA in V2 (no blast match originally)
hutpor_df <- left_join(hutpor_id, clus_df, by=c('V1' = 'clus')) %>%
  drop_na('V2')
head(hutpor_df)
dim(hutpor_df)

write.table(hutpor_df, file = 'cluster_porters_hutt_arabidopsis_annotated_geneID.txt', col.names = FALSE, row.names = FALSE, quote = FALSE,  sep = '\t')
