library('tidyverse')
library('readxl')

# reading in GO terms from uniprot
tair <- read.csv('path/to/uniprot_GO_annotations_for_arabidopsis_20210430.txt', header = TRUE, sep = '\t')

# keep only the Arabidopsis IDs with the most GO terms (in case of duplicates)
tair_len <- tair %>%
  subset(select = c(yourlist, Gene.ontology.IDs)) %>% 
  rename(go_id = Gene.ontology.IDs) %>%
  mutate(golen=str_length(go_id)) %>% # make a new column of string (go term) lengths
  group_by(yourlist) %>%              # group by gene ID
  slice_max(golen) %>%                # keep only the ID with the most GO terms
  filter(golen != 0) %>%              # drop any row with no GO terms
  distinct() %>%                      # retain only unique rows
  select(-golen)                      # drop the length column

# write this out and separate the go_id in excel!
write.table(tair_len, file = 'GO_linear.txt', col.names = F, row.names = F, quote = F,  sep = '\t')

# bring back in as tab separated
in_tair <- read.csv('GO_linear.txt', header = F, sep = '\t')

library("reshape")

# use the Melt function to move all the column values to rows.
# retains the information as to which column the GO terms came from.
go_df_fil <- melt(in_tair, id.vars = c('V1')) %>%
  select(-variable) %>%
  mutate(golen=str_length(value)) %>% # make a new column of string (go term) lengths (zero when no mapping)
  filter(golen != 0) %>%             # drop any row with no GO terms
  select(-golen)

# read in the background (with Araport annotations)
bk <- read.table('background_genes_araport.txt', header = FALSE, sep = '\t')

# remove the isoform identifier from the background output
bk_spl <- separate(bk, V2, into = c('Locus', 'iso'), sep = '\\.') %>%
  select(-iso)
head(bk_spl)

# left join the background and go mapping
bk_fin <- left_join(bk_spl, go_df_fil, by = c('Locus' = 'V1')) %>%
  select(-Locus) %>%
  drop_na()

# write out the background with GO annotations
write.table(bk_fin, 'background_genes_go.txt', col.names = F, row.names = F, quote = F,  sep = '\t')

# read in the hutt_monroi - crithmifolius cluster
hut_cri <- read.table('cluster_hutt_crithmifolius_arabidopsis_annotated_geneID.txt', header = FALSE, sep = '\t')

# remove the Arabidopsis IDs
hut_cri_cl <- hut_cri %>%
  select(V1)

# join the cluster to the background but keep both cols
merg <- left_join(hut_cri_cl, bk_fin, by = 'V1') %>%
  drop_na()
head(merg)

write.table(merg,'GO_annotations_for_AgriGO_hutt_crith.txt', col.names = F, row.names = F, quote = F, sep = '\t')

# read in the hutt_monroi - porters_monroi cluster
hut_por <- read.table('cluster_porters_hutt_arabidopsis_annotated_geneID.txt', header = FALSE, sep = '\t')
head(hut_por)

# remove the Arabidopsis IDs
hut_por_cl <- hut_por %>%
  select(V1)

# join the cluster to the background but keep both cols
mer <- left_join(hut_por_cl, bk_fin, by = 'V1') %>%
  drop_na()
head(mer)

write.table(mer,'GO_annotations_for_AgriGO_hutt_porters.txt', col.names = F, row.names = F, quote = F, sep = '\t')

# read in the hutt_monroi - lobulatus cluster
hut_lob <- read.table('cluster_hutt_lobulatus_arabidopsis_annotated_geneID.txt', header = FALSE, sep = '\t')

# remove the Arabidopsis IDs
hut_lob_cl <- hut_lob %>%
  select(V1)

# join the cluster to the background but keep both cols
me <- left_join(hut_lob_cl, bk_fin, by = 'V1') %>%
  drop_na()

write.table(me,'GO_annotations_for_AgriGO_hutt_lobulatus.txt', col.names = F, row.names = F, quote = F, sep = '\t')




