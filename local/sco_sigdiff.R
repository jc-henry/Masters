library('matrixTests')
library('tidyverse')
library('gplots')

# read in the single-copy orthogroup TPM dataframe
df_tpm <- read.table('sc_orthologue_TPM.tsv', sep = '\t', header = TRUE)

# create a vector of gene(cluster) names
gene_id <- as.data.frame(df_tpm$Gene)

# convert dataframe to a numeric matrix
mat_no_id <- data.matrix(df_tpm[,2:13])

# perform a Kruskal-Wallis rank sum test on each row of the matrix (group=taxa)
kw_group <- row_kruskalwallis(mat_no_id[,1:12], c("porters_monroi","porters_monroi","porters_monroi","crithmifolius","crithmifolius","crithmifolius","hutt_monroi","hutt_monroi","hutt_monroi","lobulatus","lobulatus","lobulatus"))

# merge the two dataframes
mer_df <- cbind(gene_id,  kw_group)

# rename the gene id column
mer_df_rename <- mer_df %>%
  rename(gene = 'df_tpm$Gene')

# extract genes with pvalue < 0.05
less_05 <- subset(mer_df_rename, mer_df_rename[,6] < 0.05)

# remove all cols except Gene id
clean_df <- subset(less_05, select = -c(obs.tot, obs.groups, df, statistic, pvalue)) 

# export the file
write.table(clean_df, 'sig_diff_genes.tsv', quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

# read in the gene TPM values
sc_tpm <- read.table('sc_orthologue_TPM.tsv', header = TRUE, sep = '\t')

# grep and merge the files
clus_idx2 <- sapply(clean_df$gene, grep, sc_tpm$Gene) # create an index 
clus_idx1 <- sapply(seq_along(clus_idx2), function(i) rep(i, length(clus_idx2[[i]]))) #duplicate the original indices so the results align
clus_merge <- cbind(clean_df[unlist(clus_idx1),,drop=F], sc_tpm[unlist(clus_idx2),,drop=F]) # merge the datasets with cbind aligned on the new indices

#remove unnecessary "gene" column
# remove all cols except Gene id
fin_df <- subset(clus_merge, select = -c(gene)) 

# in case want it later
write.table(fin_df, 'sigdiff_orthologue_TPM.txt', col.names = T, sep = '\t', quote = F)

### make a heatmap of sigdiff genes
# setup with log scaling (necessary for making heatmap).
copy_df <- fin_df
rownames(copy_df) <- NULL
name_df <- column_to_rownames(copy_df, var = "Gene")
fin_mat <- data.matrix(name_df)
Log_LC <- log10(fin_mat+1)
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

