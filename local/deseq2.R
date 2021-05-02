# describes differential gene expression analysis using DESeq2
# shows one comparison (i.e. taxons A & B compared with taxons C & D)
# need to run twice more with changed groupings (e.g. taxons A & C vs taxons B & D)

library(tximport)
library(tibble)
library(ggplot2)
library(DESeq2)

# set the path to the dir holding RSEM gene count directories for each sample
dir  <- "path/to/allranunculus_rsem_mapping_lobulatus"  

# read in the text file describing the sample relationships and comparisons (e.g. both R. monroi grouped together)
# this file is adjusted for each grouping comparison
samples <- read.table('taxonA_taxonB_cluster_rsem_samples.txt', header = TRUE)

# set path to the RSEM quantification files
files <- file.path(dir, samples$directory, "RSEM.genes.results")

# path to files based on the correct samples.txt column
names(files) <- samples$sample
all(file.exists(files))

# import the quantification data for DESeq2 using tximport
txi <- tximport(files, type = "rsem", txIn = TRUE, txOut = TRUE)

# construct the DESeqDataSet from the txi object and sample info in samples
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ cluster)

############################################
### this step only carried out for the initial analyses

# pre-filter low count genes to increase speed and reduce memory. Keep only rows with at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

###########################################
### this step only carried out for second analyses with removal of genes with any read count less than 1 

# correct for library sizes
dds <- estimateSizeFactors(dds)

# remove any rows with any zero read counts
idx <- rowSums( counts(dds, normalized=TRUE) >= 1 ) >= 12
dds <- dds[idx,]

# get the list of genes in which all samples map at least one read
plus_one <- as.data.frame(rownames(dds))
write.table(plus_one, 'genes_with_all_mapping.txt', col.names = F, row.names = F, quote = F)

###########################################
### the remainder is common to all analyses

# set the control condition
dds$cluster <- relevel(dds$cluster, ref = "taxonC_taxonD")

# run the differential expression analysis
dds <- DESeq(dds)

# get only genes with significance below FDR 0.05
res05 <- results(dds, alpha = 0.05)

# show a summary of results such as number of up and downregulated genes
summary(res05)

# subset the results to filter by padj. This filters out all the transcripts with padj value of NA
res05Sig <- subset(res05, padj < 0.05)

# export the results as tsv. And export just the gene IDs
res05Sig_df <- as.data.frame(res05Sig)
res05Sig_ids <- rownames_to_column(res05Sig_df, var = 'gene_id')
write.table(res05Sig_ids, file="cluster_cluster_taxonA_taxonB_vs_taxonC_taxonD_padj0.5_nofc.txt", sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)
write.table(res05Sig_ids$gene_id, file="geneIDonly_cluster_taxonA_taxonB_vs_taxonC_taxonD_padj0.5_nofc.txt", sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)



### generate a volcano plot of expression
res05df <- as.data.frame(res05)
# add a column for regulation
res05df$diffexpressed <- "No change"
# if log2Foldchange < 0 and FDR < 0.05, set as "DOWN"
res05df$diffexpressed[res05df$log2FoldChange < 0 & res05df$padj < 0.05] <- "Downregulated"
# if log2Foldchange > 0 and FDR < 0.05, set as "UP"
res05df$diffexpressed[res05df$log2FoldChange > 0 & res05df$padj < 0.05] <- "Upregulated"

ggplot(data=res05df, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed)) +
  geom_point() +
  theme_minimal() +
  guides(fill = guide_legend(reverse = TRUE)) +
  scale_color_manual(values=c("#5FA1F7", "black", "#9B1F1A")) +
  theme(legend.title = element_blank()) +
  guides(color = guide_legend(reverse = TRUE)) + 
  theme(text = element_text(size=15)) +
  labs(y='Significance (-log10 adjusted p-value)', x='Fold change (log2)')



