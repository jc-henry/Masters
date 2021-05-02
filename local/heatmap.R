library(tidyverse)
library(gplots)
library(matrixTests)

# set directory paths
dir_tpm <- "path/to/allranunculus_rsem_mapping_taxon"                    # the sample gene-level TPM files

# make objects of the TPM files
porters_monroi_1 <- paste(dir_tpm, "porters_monroi_1_rsem.txt", sep="/")
porters_monroi_2 <- paste(dir_tpm, "porters_monroi_2_rsem.txt", sep="/")
porters_monroi_3 <- paste(dir_tpm, "porters_monroi_3_rsem.txt", sep="/")
crithmifolius_1 <- paste(dir_tpm, "crithmifolius_1_rsem.txt", sep="/")
crithmifolius_2 <- paste(dir_tpm, "crithmifolius_2_rsem.txt", sep="/")
crithmifolius_3 <- paste(dir_tpm, "crithmifolius_3_rsem.txt", sep="/")
hutt_monroi_1 <- paste(dir_tpm, "hutt_monroi_1_rsem.txt", sep="/")
hutt_monroi_2 <- paste(dir_tpm, "hutt_monroi_2_rsem.txt", sep="/")
hutt_monroi_3 <- paste(dir_tpm, "hutt_monroi_3_rsem.txt", sep="/")
lobulatus_1 <- paste(dir_tpm, "lobulatus_1_rsem.txt", sep="/")
lobulatus_2 <- paste(dir_tpm, "lobulatus_2_rsem.txt", sep="/")
lobulatus_3 <- paste(dir_tpm, "lobulatus_3_rsem.txt", sep="/")

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

# retain the "gene_id" and "TPM" columns
porters_monroi_1_trim <- subset(porters_monroi_1_tpm, select = c(gene_id, TPM)) 
porters_monroi_2_trim <- subset(porters_monroi_2_tpm, select = c(gene_id, TPM)) 
porters_monroi_3_trim <- subset(porters_monroi_3_tpm, select = c(gene_id, TPM))
hutt_monroi_1_trim <- subset(hutt_monroi_1_tpm, select = c(gene_id, TPM)) 
hutt_monroi_2_trim <- subset(hutt_monroi_2_tpm, select = c(gene_id, TPM)) 
hutt_monroi_3_trim <- subset(hutt_monroi_3_tpm, select = c(gene_id, TPM)) 
crithmifolius_1_trim <- subset(crithmifolius_1_tpm, select = c(gene_id, TPM)) 
crithmifolius_2_trim <- subset(crithmifolius_2_tpm, select = c(gene_id, TPM)) 
crithmifolius_3_trim <- subset(crithmifolius_3_tpm, select = c(gene_id, TPM))
lobulatus_1_trim <- subset(lobulatus_1_tpm, select = c(gene_id, TPM)) 
lobulatus_2_trim <- subset(lobulatus_2_tpm, select = c(gene_id, TPM)) 
lobulatus_3_trim <- subset(lobulatus_3_tpm, select = c(gene_id, TPM))

# rename the columns
porters_monroi_1_name <- porters_monroi_1_trim %>%
  rename(porters_monroi_1 = TPM) 
porters_monroi_2_name <- porters_monroi_2_trim %>%
  rename(porters_monroi_2 = TPM) 
porters_monroi_3_name <- porters_monroi_3_trim %>%
  rename(porters_monroi_3 = TPM) 
crithmifolius_1_name <- crithmifolius_1_trim %>%
  rename(crithmifolius_1 = TPM) 
crithmifolius_2_name <- crithmifolius_2_trim %>%
  rename(crithmifolius_2 = TPM) 
crithmifolius_3_name <- crithmifolius_3_trim %>%
  rename(crithmifolius_3 = TPM) 
hutt_monroi_1_name <- hutt_monroi_1_trim %>%
  rename(hutt_monroi_1 = TPM) 
hutt_monroi_2_name <- hutt_monroi_2_trim %>%
  rename(hutt_monroi_2 = TPM) 
hutt_monroi_3_name <- hutt_monroi_3_trim %>%
  rename(hutt_monroi_3 = TPM) 
lobulatus_1_name <- lobulatus_1_trim %>%
  rename(lobulatus_1 = TPM) 
lobulatus_2_name <- lobulatus_2_trim %>%
  rename(lobulatus_2 = TPM) 
lobulatus_3_name <- lobulatus_3_trim %>%
  rename(lobulatus_3 = TPM) 

# make the dataframe with all TPM values against the renamed genes
df_tpm <- Reduce(merge,list(porters_monroi_1_name,porters_monroi_2_name,porters_monroi_3_name,
                            crithmifolius_1_name,crithmifolius_2_name,crithmifolius_3_name,
                            hutt_monroi_1_name,hutt_monroi_2_name,hutt_monroi_3_name,
                            lobulatus_1_name,lobulatus_2_name,lobulatus_3_name))


# export the file
write.table(df_tpm, 'all_gene_TPM.tsv', quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

# create a vector of gene(cluster) names
gene_id <- as.data.frame(df_tpm$gene_id)

# convert dataframe to a numeric matrix
mat_no_id <- data.matrix(df_tpm[,2:13])

# perform a  Kruskal-Wallis rank sum test on each row of the matrix (group=taxa)
kw_group <- row_kruskalwallis(mat_no_id[,1:12], c("porters_monroi","porters_monroi","porters_monroi","crithmifolius","crithmifolius","crithmifolius","hutt_monroi","hutt_monroi","hutt_monroi","lobulatus","lobulatus","lobulatus"))

# merge the two dataframes
mer_df <- cbind(gene_id,  kw_group)

# rename the gene id column
mer_df_rename <- mer_df %>%
  rename(gene = 'df_tpm$gene_id')

# extract genes with pvalue < 0.05
less_05 <- subset(mer_df_rename, mer_df_rename[,6] < 0.05)

# remove all cols except Gene id
sig <- subset(less_05, select = -c(obs.tot, obs.groups, df, statistic, pvalue)) 

# export the file
write.table(sig, 'sig_diff_genes.tsv', quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')




### necessary to make heat map of only genes found to be significantly differently expressed due to amount of data


# make a TPM object with the "_i" attached (needed for grep to be specific since no -x flag in R) - just in rare case more than 10 "genes" of same cluster
sig_p <- sig %>%
  mutate(Hold = paste(gene, '_i', sep = ''))

# need to do the same as above for the input df_tpm
df_tpm_p <- df_tpm %>%
  mutate(hold = paste(gene_id, '_i', sep = ''))

# grep and merge the files
clus_idx2 <- sapply(sig_p$Hold, grep, df_tpm_p$hold) # create an index
clus_idx1 <- sapply(seq_along(clus_idx2), function(i) rep(i, length(clus_idx2[[i]]))) #duplicate the original indices so the results align
clus_merge <- cbind(sig[unlist(clus_idx1),,drop=F], df_tpm[unlist(clus_idx2),,drop=F]) # merge the datasets with cbind aligned on the new indices


# remove unnecessary "gene" column
# remove all cols except Gene id
fin_df <- subset(clus_merge, select = -c(gene))
head(fin_df)

# write out the final df if case needed later
write.table(fin_df, 'cleaned_final_sigdiff.txt', quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

### heatmap of all sigdiff genes
# setup with log scaling (necessary for making heatmap).
copy_df <- fin_df
rownames(copy_df) <- NULL
name_df <- column_to_rownames(copy_df, var = "gene_id")
fin_mat <- data.matrix(name_df)
Log_LC <- log10(fin_mat+1)
SC_LC <- scale(Log_LC,scale = TRUE)
y <- data.matrix(SC_LC)

## Row- and column-wise clustering 
hr <- hclust(as.dist(1-cor(t(y), method="pearson")), method="complete") # row clustering (t transposes the data since cor correlates on columns)
hc <- hclust(as.dist(1-cor(y, method="spearman")), method="complete") # column clustering


pdf("sigdiff_gene_heatmap.pdf")
mycol <- colorpanel(100, "#83A552", "black", "red")

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
dev.off()



