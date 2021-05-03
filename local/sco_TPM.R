library('tidyverse')
library('gplots')

# set directory paths
dir_tpm <- "path/to/allranunculus_rsem_mapping_taxon"           # the sample gene-level TPM files
dir_list <- "path/to/ortholists"                                # the lists of SCOs per taxon 

# make objects of the TPM files
porters_monroi_1 <- paste(dir_tpm, "porters_monroi_1_geneTPM.tsv", sep="/")
porters_monroi_2 <- paste(dir_tpm, "porters_monroi_2_geneTPM.tsv", sep="/")
porters_monroi_3 <- paste(dir_tpm, "porters_monroi_3_geneTPM.tsv", sep="/")
crithmifolius_1 <- paste(dir_tpm, "crithmifolius_1_geneTPM.tsv", sep="/")
crithmifolius_2 <- paste(dir_tpm, "crithmifolius_2_geneTPM.tsv", sep="/")
crithmifolius_3 <- paste(dir_tpm, "crithmifolius_3_geneTPM.tsv", sep="/")
hutt_monroi_1 <- paste(dir_tpm, "hutt_monroi_1_geneTPM.tsv", sep="/")
hutt_monroi_2 <- paste(dir_tpm, "hutt_monroi_2_geneTPM.tsv", sep="/")
hutt_monroi_3 <- paste(dir_tpm, "hutt_monroi_3_geneTPM.tsv", sep="/")
lobulatus_1 <- paste(dir_tpm, "lobulatus_1_geneTPM.tsv", sep="/")
lobulatus_2 <- paste(dir_tpm, "lobulatus_2_geneTPM.tsv", sep="/")
lobulatus_3 <- paste(dir_tpm, "lobulatus_3_geneTPM.tsv", sep="/")

# make objects of the ortholist files
porters_monroi_list <- paste(dir_list, "ortholist_porters_monroi.txt", sep="/")
crithmifolius_list <- paste(dir_list, "ortholist_crithmifolius.txt", sep="/")
hutt_monroi_list <- paste(dir_list, "ortholist_hutt_monroi.txt", sep="/")
lobulatus_list <- paste(dir_list, "ortholist_lobulatus.txt", sep="/")

# read in the TPM files
porters_monroi_1_tpm <- read.table(file=porters_monroi_1, sep = "\t", header = TRUE)
porters_monroi_2_tpm <- read.table(file=porters_monroi_2, sep = "\t", header = TRUE)
porters_monroi_3_tpm <- read.table(file=porters_monroi_3, sep = "\t", header = TRUE)
crithmifolius_1_tpm <- read.table(file=crithmifolius_1, sep = "\t", header = TRUE)
crithmifolius_2_tpm <- read.table(file=crithmifolius_2, sep = "\t", header = TRUE)
crithmifolius_3_tpm <- read.table(file=crithmifolius_3, sep = "\t", header = TRUE)
hutt_monroi_1_tpm <- read.table(file=hutt_monroi_1, sep = "\t", header = TRUE)
hutt_monroi_2_tpm <- read.table(file=hutt_monroi_2, sep = "\t", header = TRUE)
hutt_monroi_3_tpm <- read.table(file=hutt_monroi_3, sep = "\t", header = TRUE)
lobulatus_1_tpm <- read.table(file=lobulatus_1, sep = "\t", header = TRUE)
lobulatus_2_tpm <- read.table(file=lobulatus_2, sep = "\t", header = TRUE)
lobulatus_3_tpm <- read.table(file=lobulatus_3, sep = "\t", header = TRUE)

# read in the ortholist files
porters_monroi_orth <- read.table(file=porters_monroi_list, sep = '\t', header = FALSE)
crithmifolius_orth <- read.table(file=crithmifolius_list, sep = '\t', header = FALSE)
hutt_monroi_orth <- read.table(file=hutt_monroi_list, sep = '\t', header = FALSE)
lobulatus_orth <- read.table(file=lobulatus_list, sep = '\t', header = FALSE)

# make a TPM object with the "_i" attached (needed for grep to be specific since no -x flag in R) - just in rare case more than 10 "genes"
porters_monroi_1_tpm_p <- porters_monroi_1_tpm %>%
  mutate(Hold = paste(Gene, '_i', sep = ''))
porters_monroi_2_tpm_p <- porters_monroi_2_tpm %>%
  mutate(Hold = paste(Gene, '_i', sep = ''))
porters_monroi_3_tpm_p <- porters_monroi_3_tpm %>%
  mutate(Hold = paste(Gene, '_i', sep = ''))
crithmifolius_1_tpm_p <- crithmifolius_1_tpm %>%
  mutate(Hold = paste(Gene, '_i', sep = ''))
crithmifolius_2_tpm_p <- crithmifolius_2_tpm %>%
  mutate(Hold = paste(Gene, '_i', sep = ''))
crithmifolius_3_tpm_p <- crithmifolius_3_tpm %>%
  mutate(Hold = paste(Gene, '_i', sep = ''))
hutt_monroi_1_tpm_p <- hutt_monroi_1_tpm %>%
  mutate(Hold = paste(Gene, '_i', sep = ''))
hutt_monroi_2_tpm_p <- hutt_monroi_2_tpm %>%
  mutate(Hold = paste(Gene, '_i', sep = ''))
hutt_monroi_3_tpm_p <- hutt_monroi_3_tpm %>%
  mutate(Hold = paste(Gene, '_i', sep = ''))
lobulatus_1_tpm_p <- lobulatus_1_tpm %>%
  mutate(Hold = paste(Gene, '_i', sep = ''))
lobulatus_2_tpm_p <- lobulatus_2_tpm %>%
  mutate(Hold = paste(Gene, '_i', sep = ''))
lobulatus_3_tpm_p <- lobulatus_3_tpm %>%
  mutate(Hold = paste(Gene, '_i', sep = ''))

# search for the Trinity gene ID in the ortholist then bind matches together
porters_monroi_1_idx2 <- sapply(porters_monroi_1_tpm_p$Hold, grep, porters_monroi_orth$V1, fixed=TRUE) # create an index 
porters_monroi_1_idx1 <- sapply(seq_along(porters_monroi_1_idx2), function(i) rep(i, length(porters_monroi_1_idx2[[i]]))) #duplicate the original indices so the results align
porters_monroi_1_merge <- cbind(porters_monroi_1_tpm_p[unlist(porters_monroi_1_idx1),,drop=F], porters_monroi_orth[unlist(porters_monroi_1_idx2),,drop=F]) # merge the datasets with cbind aligned on the new indices

porters_monroi_2_idx2 <- sapply(porters_monroi_2_tpm_p$Hold, grep, porters_monroi_orth$V1, fixed=TRUE) 
porters_monroi_2_idx1 <- sapply(seq_along(porters_monroi_2_idx2), function(i) rep(i, length(porters_monroi_2_idx2[[i]]))) 
porters_monroi_2_merge <- cbind(porters_monroi_2_tpm_p[unlist(porters_monroi_2_idx1),,drop=F], porters_monroi_orth[unlist(porters_monroi_2_idx2),,drop=F]) 

porters_monroi_3_idx2 <- sapply(porters_monroi_3_tpm_p$Hold, grep, porters_monroi_orth$V1, fixed=TRUE) 
porters_monroi_3_idx1 <- sapply(seq_along(porters_monroi_3_idx2), function(i) rep(i, length(porters_monroi_3_idx2[[i]]))) 
porters_monroi_3_merge <- cbind(porters_monroi_3_tpm_p[unlist(porters_monroi_3_idx1),,drop=F], porters_monroi_orth[unlist(porters_monroi_3_idx2),,drop=F]) 

crithmifolius_1_idx2 <- sapply(crithmifolius_1_tpm_p$Hold, grep, crithmifolius_orth$V1, fixed=TRUE) 
crithmifolius_1_idx1 <- sapply(seq_along(crithmifolius_1_idx2), function(i) rep(i, length(crithmifolius_1_idx2[[i]]))) 
crithmifolius_1_merge <- cbind(crithmifolius_1_tpm_p[unlist(crithmifolius_1_idx1),,drop=F], crithmifolius_orth[unlist(crithmifolius_1_idx2),,drop=F]) 

crithmifolius_2_idx2 <- sapply(crithmifolius_2_tpm_p$Hold, grep, crithmifolius_orth$V1, fixed=TRUE) 
crithmifolius_2_idx1 <- sapply(seq_along(crithmifolius_2_idx2), function(i) rep(i, length(crithmifolius_2_idx2[[i]]))) 
crithmifolius_2_merge <- cbind(crithmifolius_2_tpm_p[unlist(crithmifolius_2_idx1),,drop=F], crithmifolius_orth[unlist(crithmifolius_2_idx2),,drop=F]) 

crithmifolius_3_idx2 <- sapply(crithmifolius_3_tpm_p$Hold, grep, crithmifolius_orth$V1, fixed=TRUE) 
crithmifolius_3_idx1 <- sapply(seq_along(crithmifolius_3_idx2), function(i) rep(i, length(crithmifolius_3_idx2[[i]]))) 
crithmifolius_3_merge <- cbind(crithmifolius_3_tpm_p[unlist(crithmifolius_3_idx1),,drop=F], crithmifolius_orth[unlist(crithmifolius_3_idx2),,drop=F]) 

hutt_monroi_1_idx2 <- sapply(hutt_monroi_1_tpm_p$Hold, grep, hutt_monroi_orth$V1, fixed=TRUE) 
hutt_monroi_1_idx1 <- sapply(seq_along(hutt_monroi_1_idx2), function(i) rep(i, length(hutt_monroi_1_idx2[[i]]))) 
hutt_monroi_1_merge <- cbind(hutt_monroi_1_tpm_p[unlist(hutt_monroi_1_idx1),,drop=F], hutt_monroi_orth[unlist(hutt_monroi_1_idx2),,drop=F]) 

hutt_monroi_2_idx2 <- sapply(hutt_monroi_2_tpm_p$Hold, grep, hutt_monroi_orth$V1, fixed=TRUE) 
hutt_monroi_2_idx1 <- sapply(seq_along(hutt_monroi_2_idx2), function(i) rep(i, length(hutt_monroi_2_idx2[[i]]))) 
hutt_monroi_2_merge <- cbind(hutt_monroi_2_tpm_p[unlist(hutt_monroi_2_idx1),,drop=F], hutt_monroi_orth[unlist(hutt_monroi_2_idx2),,drop=F]) 

hutt_monroi_3_idx2 <- sapply(hutt_monroi_3_tpm_p$Hold, grep, hutt_monroi_orth$V1, fixed=TRUE) 
hutt_monroi_3_idx1 <- sapply(seq_along(hutt_monroi_3_idx2), function(i) rep(i, length(hutt_monroi_3_idx2[[i]]))) 
hutt_monroi_3_merge <- cbind(hutt_monroi_3_tpm_p[unlist(hutt_monroi_3_idx1),,drop=F], hutt_monroi_orth[unlist(hutt_monroi_3_idx2),,drop=F]) 

lobulatus_1_idx2 <- sapply(lobulatus_1_tpm_p$Hold, grep, lobulatus_orth$V1, fixed=TRUE)
lobulatus_1_idx1 <- sapply(seq_along(lobulatus_1_idx2), function(i) rep(i, length(lobulatus_1_idx2[[i]]))) 
lobulatus_1_merge <- cbind(lobulatus_1_tpm_p[unlist(lobulatus_1_idx1),,drop=F], lobulatus_orth[unlist(lobulatus_1_idx2),,drop=F]) 

lobulatus_2_idx2 <- sapply(lobulatus_2_tpm_p$Hold, grep, lobulatus_orth$V1, fixed=TRUE)
lobulatus_2_idx1 <- sapply(seq_along(lobulatus_2_idx2), function(i) rep(i, length(lobulatus_2_idx2[[i]]))) 
lobulatus_2_merge <- cbind(lobulatus_2_tpm_p[unlist(lobulatus_2_idx1),,drop=F], lobulatus_orth[unlist(lobulatus_2_idx2),,drop=F]) 

lobulatus_3_idx2 <- sapply(lobulatus_3_tpm_p$Hold, grep, lobulatus_orth$V1, fixed=TRUE)
lobulatus_3_idx1 <- sapply(seq_along(lobulatus_3_idx2), function(i) rep(i, length(lobulatus_3_idx2[[i]]))) 
lobulatus_3_merge <- cbind(lobulatus_3_tpm_p[unlist(lobulatus_3_idx1),,drop=F], lobulatus_orth[unlist(lobulatus_3_idx2),,drop=F]) 

# drop the "Gene" and "Hold columns to leave only TPM and ortholist ID
porters_monroi_1_nogene <- subset(porters_monroi_1_merge, select = -c(Gene, Hold)) 
porters_monroi_2_nogene <- subset(porters_monroi_2_merge, select = -c(Gene, Hold)) 
porters_monroi_3_nogene <- subset(porters_monroi_3_merge, select = -c(Gene, Hold)) 
crithmifolius_1_nogene <- subset(crithmifolius_1_merge, select = -c(Gene, Hold)) 
crithmifolius_2_nogene <- subset(crithmifolius_2_merge, select = -c(Gene, Hold)) 
crithmifolius_3_nogene <- subset(crithmifolius_3_merge, select = -c(Gene, Hold)) 
hutt_monroi_1_nogene <- subset(hutt_monroi_1_merge, select = -c(Gene, Hold)) 
hutt_monroi_2_nogene <- subset(hutt_monroi_2_merge, select = -c(Gene, Hold)) 
hutt_monroi_3_nogene <- subset(hutt_monroi_3_merge, select = -c(Gene, Hold)) 
lobulatus_1_nogene <- subset(lobulatus_1_merge, select = -c(Gene, Hold)) 
lobulatus_2_nogene <- subset(lobulatus_2_merge, select = -c(Gene, Hold)) 
lobulatus_3_nogene <- subset(lobulatus_3_merge, select = -c(Gene, Hold)) 

# leave only TPM and OrthoFinder cluster ID
# splits at all "_" but discards info after the first
porters_monroi_1_split <- porters_monroi_1_nogene %>%
  separate(V1, into = c('Gene'), sep="_") 
porters_monroi_2_split <- porters_monroi_2_nogene %>%
  separate(V1, into = c('Gene'), sep="_") 
porters_monroi_3_split <- porters_monroi_3_nogene %>%
  separate(V1, into = c('Gene'), sep="_") 
crithmifolius_1_split <- crithmifolius_1_nogene %>%
  separate(V1, into = c('Gene'), sep="_") 
crithmifolius_2_split <- crithmifolius_2_nogene %>%
  separate(V1, into = c('Gene'), sep="_") 
crithmifolius_3_split <- crithmifolius_3_nogene %>%
  separate(V1, into = c('Gene'), sep="_") 
hutt_monroi_1_split <- hutt_monroi_1_nogene %>%
  separate(V1, into = c('Gene'), sep="_") 
hutt_monroi_2_split <- hutt_monroi_2_nogene %>%
  separate(V1, into = c('Gene'), sep="_") 
hutt_monroi_3_split <- hutt_monroi_3_nogene %>%
  separate(V1, into = c('Gene'), sep="_") 
lobulatus_1_split <- lobulatus_1_nogene %>%
  separate(V1, into = c('Gene'), sep="_") 
lobulatus_2_split <- lobulatus_2_nogene %>%
  separate(V1, into = c('Gene'), sep="_") 
lobulatus_3_split <- lobulatus_3_nogene %>%
  separate(V1, into = c('Gene'), sep="_") 

# reorder so "Gene" (cluster ID) is the first column
porters_monroi_1_reord <- porters_monroi_1_split[c('Gene', 'TPM')]
porters_monroi_2_reord <- porters_monroi_2_split[c('Gene', 'TPM')]
porters_monroi_3_reord <- porters_monroi_3_split[c('Gene', 'TPM')]
crithmifolius_1_reord <- crithmifolius_1_split[c('Gene', 'TPM')]
crithmifolius_2_reord <- crithmifolius_2_split[c('Gene', 'TPM')]
crithmifolius_3_reord <- crithmifolius_3_split[c('Gene', 'TPM')]
hutt_monroi_1_reord <- hutt_monroi_1_split[c('Gene', 'TPM')]
hutt_monroi_2_reord <- hutt_monroi_2_split[c('Gene', 'TPM')]
hutt_monroi_3_reord <- hutt_monroi_3_split[c('Gene', 'TPM')]
lobulatus_1_reord <- lobulatus_1_split[c('Gene', 'TPM')]
lobulatus_2_reord <- lobulatus_2_split[c('Gene', 'TPM')]
lobulatus_3_reord <- lobulatus_3_split[c('Gene', 'TPM')]

# rename the TPM columns with sample ID
porters_monroi_1_name <- porters_monroi_1_reord %>%
  rename(porters_monroi_1 = TPM) 
porters_monroi_2_name <- porters_monroi_2_reord %>%
  rename(porters_monroi_2 = TPM) 
porters_monroi_3_name <- porters_monroi_3_reord %>%
  rename(porters_monroi_3 = TPM) 
crithmifolius_1_name <- crithmifolius_1_reord %>%
  rename(crithmifolius_1 = TPM) 
crithmifolius_2_name <- crithmifolius_2_reord %>%
  rename(crithmifolius_2 = TPM) 
crithmifolius_3_name <- crithmifolius_3_reord %>%
  rename(crithmifolius_3 = TPM) 
hutt_monroi_1_name <- hutt_monroi_1_reord %>%
  rename(hutt_monroi_1 = TPM) 
hutt_monroi_2_name <- hutt_monroi_2_reord %>%
  rename(hutt_monroi_2 = TPM) 
hutt_monroi_3_name <- hutt_monroi_3_reord %>%
  rename(hutt_monroi_3 = TPM) 
lobulatus_1_name <- lobulatus_1_reord %>%
  rename(lobulatus_1 = TPM) 
lobulatus_2_name <- lobulatus_2_reord %>%
  rename(lobulatus_2 = TPM) 
lobulatus_3_name <- lobulatus_3_reord %>%
  rename(lobulatus_3 = TPM) 

# make the dataframe with all TPM values against the renamed genes
tpm_df <- Reduce(merge,list(porters_monroi_1_name,porters_monroi_2_name,porters_monroi_3_name,
                            crithmifolius_1_name,crithmifolius_2_name,crithmifolius_3_name,
                            hutt_monroi_1_name,hutt_monroi_2_name,hutt_monroi_3_name,
                            lobulatus_1_name,lobulatus_2_name,lobulatus_3_name))


# export the file
write.table(tpm_df, 'sc_orthologue_TPM.tsv', quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

### make a heatmap of all ortho gene TPM
# setup with log scaling (necessary for making heatmap).
copy_df <- tpm_df
rownames(copy_df) <- NULL
name_df <- column_to_rownames(copy_df, var = "Gene")
fin_mat <- data.matrix(name_df)
Log_LC <- log10(fin_mat+1) # add 1 to values to avoid decimals which will have large negative log values
SC_LC <- scale(Log_LC,scale = TRUE)
y <- data.matrix(SC_LC)

## Row- and column-wise clustering 
hr <- hclust(as.dist(1-cor(t(y), method="pearson")), method="complete") # row clustering (t transposes the data since cor correlates on columns)
hc <- hclust(as.dist(1-cor(y, method="spearman")), method="complete") # column clustering

mycol <- colorpanel(75, "#56B4E9", "black", "#F0E442") 
heatmap.2(
  y,
  Rowv = as.dendrogram(hr),
  Colv = as.dendrogram(hc),
  col = mycol,
  density.info = "none",
  trace = "none",
  dendrogram = "both",
  scale = "row",
  labRow = FALSE,
  labCol = NULL,
  srtCol=45,  adjCol = c(1,1),
  cexCol = 1,
  margins = c(10,7)
)

