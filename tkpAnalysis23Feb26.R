##############################
############ RNA SEQ #########
##############################

# Install DESeq2 if not already installed
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("DESeq2")
#install.packages("DBI")
# Load the necessary libraries
library(DESeq2)
library(ggplot2)
library("data.table")

setwd("/home/aa/pigmentation/tkpAnalysis/petalRNASeq/")

# Define function to read STAR count files (reverse strand data)
read_star_counts <- function(filepath) {
  data <- read.table(filepath, header = FALSE, skip = 4)
  counts <- data[, c(1, 4)]  # Extracting gene ID and reverse strand counts (column 4)
  colnames(counts) <- c("GeneID", sub("_ReadsPerGene.out.tab", "", basename(filepath)))
  return(counts)
}
# Define paths to all STAR output files
file_paths <- list.files("/home/aa/pigmentation/tkpAnalysis/petalRNASeq/data/starOutput/", pattern = "_ReadsPerGene.out.tab$", full.names = TRUE)
# Read all count files
count_list <- lapply(file_paths, read_star_counts)
# Merge all count data into a single data frame
counts_data <- Reduce(function(x, y) merge(x, y, by = "GeneID"), count_list)
# Set rownames as Gene IDs and remove GeneID column
rownames(counts_data) <- counts_data$GeneID
counts_data <- counts_data[, -1]

# Load the metadata
metadata <- read.table("data/SampleMetadata.txt", header = TRUE, sep = "\t")
# Ensure that metadata sample IDs match count matrix column names
metadata <- metadata[match(colnames(counts_data), metadata$sampleID), ]
# Check if all samples in the metadata have corresponding counts
if (any(is.na(metadata$sampleID))) {
  stop("Some sample IDs in metadata do not match the read counts.")}

# Prepare for DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = metadata, design = ~ color) #  design = ~ ploidy + ecotype + color

# Keep genes with at least some counts across all samples
hist(rowSums(counts(dds) == 0),breaks = 30)
hist(log10(rowSums(counts(dds))),breaks = 100)
dds <- dds[rowSums(counts(dds)) >= 1200, ]
hist(log10(rowSums(counts(dds))),breaks = 100)
dds <- dds[rowSums(counts(dds) <= 5) <= 6, ]
hist(log10(rowSums(counts(dds))),breaks = 100)


# Normalize counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
hist(log10(rowSums(normalized_counts)),breaks = 100)

pink_samples <- which(metadata$color == "pink")
white_samples <- which(metadata$color == "white")
# Calculate variance for pink and white groups separately
pink_variances <- apply(normalized_counts[, pink_samples], 1, var)
white_variances <- apply(normalized_counts[, white_samples], 1, var)
hist(log10(pink_variances),breaks = 100)
hist(log10(white_variances),breaks = 100)
# Define threshold for filtering high-variance genes
variance_threshold <- quantile(c(pink_variances, white_variances), 0.95)
# Filter out genes with variance higher than the threshold in either group
high_var_genes <- (pink_variances > variance_threshold) | (white_variances > variance_threshold)
dds <- dds[!high_var_genes, ]
hist(log10(rowSums(counts(dds))),breaks = 100)

# Run the DESeq2 pipeline
dds <- DESeq(dds)
resultsNames(dds)

#resLFC <- lfcShrink(dds, coef="color_white_vs_pink", type="apeglm")
explan<-as.data.frame(mcols(mcols(dds), use.names=TRUE))
###### The end of the main DeSeq2 analysis 

# Plot PCA to check overall sample clustering
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "color")
# Nicer PCA
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "color", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
# Add sample labels from the metadata
pcaData$sampleID <- metadata$sampleID  # Assuming metadata has a sampleID column
# Create a PCA plot with ggplot2, including sample labels
pdf("results/PCAbyColor4OctFilt.pdf", width = 8.5, height = 5)
ggplot(pcaData, aes(PC1, PC2, color = color, label = sampleID)) +
  geom_point(size = 5) +
  geom_text(vjust = 1.0, hjust = -0.3) +  # Add text labels
  scale_color_manual(values=c("pink"="magenta", "white"="grey70")) +  
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("Gene expression in petals")
dev.off()

# Results with shrinkage
res_color <- lfcShrink(dds, coef="color_white_vs_pink", type="apeglm")

# Alternative: extract the results for your variable of interest (e.g., color)
#res_color <- results(dds, contrast = c("color", "pink", "white"))
# Order results by adjusted p-value
res_ordered <- res_color[order(res_color$padj), ]
# Show summary of significant results
summary(res_color)
sum(res_color$padj < 0.001, na.rm=TRUE)




###### Extract significant
# Extract significant DEGs (adjusted p-value < 0.05)
sig_DEGs <- subset(res_color, padj < 0.001)
sig_DEGs_df <- as.data.frame(sig_DEGs)
sig_DEGs_df$absLog2FC <- abs(sig_DEGs_df$log2FoldChange)
sig_DEGs_df<-sig_DEGs_df[order(sig_DEGs_df$absLog2FC,decreasing = T), ] ### REORDER HERE for plotting 
sig_DEGs_df$ID <- substr(row.names(sig_DEGs_df),1,9)
write.csv(sig_DEGs_df, file = "results/significantColor4OctFilt_0.001.csv", row.names = F)

# ### Or:
# # Order genes by adjusted p-value
# top_genes <- res_color[order(abs(res_color$log2FoldChange), decreasing = T), ][1:500, ]
# # Convert to a data frame for easier export
# top_genes_df <- as.data.frame(top_genes)
# top_genes_df$absLog2FC <- abs(top_genes_df$log2FoldChange)
# top_genes_df$ID <- substr(row.names(top_genes_df),1,9)
# write.csv(top_genes_df, file = "results/topLogFCColor_500.csv", row.names = F)

#### Annotation
# Load data
dict <- fread("/home/aa/Desktop/references/lyrataV2/functions/ALATdict.txt")
dict$AT <- substr(dict$AT, 1, 9)
ann <- fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR10orth_des_20150927.txt", h = T, quote = "")
ann2018 <- fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR11orth_des_20181231.txt", h = T, quote = "")

# Initialize an empty list to collect the rows
result_list <- list()
# Loop through sig_DEGs_df IDs
for (id in sig_DEGs_df$ID) { # id = "AL2G24170"
  # Find AT code from dict
  if (nrow(subset(dict, dict$AL == id)) > 0) {
    atcode <- subset(dict, dict$AL == id)$AT
    # Extract info from ann and ann2018
    an15 <- subset(ann, ann$`Version-2 A.t. Ortholog` == atcode)[, c(2,4,6,7)]
    an18 <- subset(ann2018, substr(ann2018$AT, 1, 9) == atcode)[1, 4:6]
    # If an15 or an18 are empty, replace with NA
    if (nrow(an15) == 0) {
      an15 <- data.table(`Version-2 A.t. Ortholog` = NA, `Gene Model Description` = NA, `Primary Gene Symbol` = NA, `All Gene Symbols` = NA)  # Adjust column names as per your ann file
    }
    if (nrow(an18) == 0) {
      an18 <- data.table(annShort = NA, annCurator = NA, annComputer = NA)  # Adjust column names as per your ann2018 file
    }
    # Merge info with sig_DEGs_df row
    result_list[[id]] <- cbind(sig_DEGs_df[sig_DEGs_df$ID == id, ], an15, an18)
  }
}
# Convert the list to a data frame
sig_DEGs_df_ann <- rbindlist(result_list, use.names = TRUE, fill = TRUE)

write.csv(sig_DEGs_df_ann, file = "results/significantColor4OctFilt_0.001_annotated.csv", row.names = F)


###### Visualize
# interesting interactive thing: not always needed
plotMA(res_color, ylim=c(-5.5,5.5))
#idx <- identify(res_color$baseMean, res_color$log2FoldChange)
#rownames(res_color)[idx]
#plotMA(resLFC, ylim=c(-8,15))




# candidate genes
plotCounts(dds, gene=which.min(res_color$padj), intgroup="color")
plotCounts(dds, gene="AL2G24200.v2.1", intgroup="color")


sig_DEGs_df_ann<-fread("results/significantColor4OctFilt_0.001_annotated.csv")

pdf("results/dotplot_significantColor4OctFilt_0.001_annotated.pdf", width = 4, height = 8)
for (id in sig_DEGs_df_ann$ID) {
  # Extract normalized counts for the specific gene
  gene_counts <- plotCounts(dds, gene=paste0(id,".v2.1"), intgroup="color", returnData=TRUE)
  # Check if gene_counts is valid
  if (!is.null(gene_counts) && nrow(gene_counts) > 0) {
    # Add the sample IDs from your metadata
    gene_counts$sampleID <- rownames(gene_counts)
    # Create a ggplot with sample IDs as labels
    p <- ggplot(gene_counts, aes(x=color, y=count, label=sampleID)) +
      geom_point(size=5, aes(color=color)) +  # Plot points
      geom_text(vjust=0.5, hjust=-0.5) +  # Add sample ID labels
      scale_color_manual(values=c("pink"="magenta", "white"="grey80")) +  
      labs(title=paste0(id," ",subset(sig_DEGs_df_ann, sig_DEGs_df_ann$ID == id)[,8][[1]]), # 9 if shrinkage not used
           x=sub(";.*", "", subset(sig_DEGs_df_ann, sig_DEGs_df_ann$ID == id)[,14][[1]]), # 15 if shrinkage not used
           y="Normalized Counts") +
      theme_minimal() #+ 
    #      scale_y_log10()  # Optionally use log scale for better visualization
    # Print the plot
    print(p)
  }
}
dev.off()






##############################################
############# ALTERNATIVE WITH WILCOXON'S TEST
# Load the necessary libraries
library(DESeq2)
library(ggplot2)
library("data.table")

setwd("/home/aa/pigmentation/tkpAnalysis/petalRNASeq/")
# Define function to read STAR count files (reverse strand data)
read_star_counts <- function(filepath) {
  data <- read.table(filepath, header = FALSE, skip = 4)
  counts <- data[, c(1, 4)]  # Extracting gene ID and reverse strand counts (column 4)
  colnames(counts) <- c("GeneID", sub("_ReadsPerGene.out.tab", "", basename(filepath)))
  return(counts)
}
# Define paths to all STAR output files
file_paths <- list.files("/home/aa/pigmentation/tkpAnalysis/petalRNASeq/data/starOutput/", pattern = "_ReadsPerGene.out.tab$", full.names = TRUE)
# Read all count files
count_list <- lapply(file_paths, read_star_counts)
# Merge all count data into a single data frame
counts_data <- Reduce(function(x, y) merge(x, y, by = "GeneID"), count_list)
# Set rownames as Gene IDs and remove GeneID column
rownames(counts_data) <- counts_data$GeneID
counts_data <- counts_data[, -1]

# Load the metadata
metadata <- read.table("data/SampleMetadata.txt", header = TRUE, sep = "\t")
# Ensure that metadata sample IDs match count matrix column names
metadata <- metadata[match(colnames(counts_data), metadata$sampleID), ]
# Check if all samples in the metadata have corresponding counts
if (any(is.na(metadata$sampleID))) {
  stop("Some sample IDs in metadata do not match the read counts.")
}

# Prepare for DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = metadata, design = ~ color) #  design = ~ ploidy + ecotype + color

# Filter
# Keep genes with at least some counts across all samples
hist(rowSums(counts(dds) == 0),breaks = 30)
hist(log10(rowSums(counts(dds))),breaks = 100)
dds <- dds[rowSums(counts(dds)) >= 240, ]
hist(log10(rowSums(counts(dds))),breaks = 100)
dds <- dds[rowSums(counts(dds) <= 1) <= 15, ]
hist(log10(rowSums(counts(dds))),breaks = 100)

# Plot PCA to check overall sample clustering
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "color")
# Nicer PCA
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "color", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
# Add sample labels from the metadata
pcaData$sampleID <- metadata$sampleID  # Assuming metadata has a sampleID column
# Create a PCA plot with ggplot2, including sample labels
pdf("results/PCAbyColor4OctFilt240_15.pdf", width = 8, height = 5)
ggplot(pcaData, aes(PC1, PC2, color = color, label = sampleID)) +
  geom_point(size = 5) +
  geom_text(vjust = 1.0, hjust = -0.3) +  # Add text labels
  scale_color_manual(values=c("pink"="magenta", "white"="grey70")) +  
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("Gene expression in petals")
dev.off()

# Normalize counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
## Optional: 
write.csv(normalized_counts,"data/normalizedCounts_filtered_240RC_15nonzero_7Oct.csv")
hist(log10(rowSums(normalized_counts)),breaks = 100)
pink_samples <- which(metadata$color == "pink")
white_samples <- which(metadata$color == "white")

### Rank-Based Testing Using Wilcoxon Rank-Sum Test ###
# Perform Wilcoxon rank-sum test for each gene
wilcox_test_results <- apply(normalized_counts, 1, function(gene_counts) {
  wilcox.test(gene_counts[pink_samples], gene_counts[white_samples])$p.value
})

wilcox_test_stat <- apply(normalized_counts, 1, function(gene_counts) {
  wilcox.test(gene_counts[pink_samples], gene_counts[white_samples])$statistic
})
# Adjust p-values using Benjamini-Hochberg correction for multiple testing
adjusted_pvals <- p.adjust(wilcox_test_results, method = "BH")
# Combine results into a data frame
wilcox_res <- data.frame(GeneID = substr(rownames(normalized_counts),1,9), pvalue = wilcox_test_results, padj = adjusted_pvals, statW = wilcox_test_stat)
# Filter for significant genes (adjusted p-value < 0.001)
sig_genes <- subset(wilcox_res, padj < 0.05)
sig_genes <- sig_genes[order(sig_genes$padj), ]
# Write significant genes to a file
write.csv(wilcox_res, file = "results/rankBasedColor_all_20601_4OctFilt.csv", row.names = FALSE)
write.csv(sig_genes, file = "results/rankBasedColor_p0.05_4OctFilt.csv", row.names = FALSE)

#### Annotation
# Load data
dict <- fread("/home/aa/Desktop/references/lyrataV2/functions/ALATdict.txt")
dict$AT <- substr(dict$AT, 1, 9)
ann <- fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR10orth_des_20150927.txt", h = T, quote = "")
ann2018 <- fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR11orth_des_20181231.txt", h = T, quote = "")

# Initialize an empty list to collect the rows
result_list <- list()
# Loop through sig_DEGs_df IDs
for (id in sig_genes$GeneID) { # id = "AL2G24170"
  # Find AT code from dict
  if (nrow(subset(dict, dict$AL == id)) > 0) {
    atcode <- subset(dict, dict$AL == id)$AT
    # Extract info from ann and ann2018
    an15 <- subset(ann, ann$`Version-2 A.t. Ortholog` == atcode)[, c(2,4,6,7)]
    an18 <- subset(ann2018, substr(ann2018$AT, 1, 9) == atcode)[1, 4:6]
    # If an15 or an18 are empty, replace with NA
    if (nrow(an15) == 0) {
      an15 <- data.table(`Version-2 A.t. Ortholog` = NA, `Gene Model Description` = NA, `Primary Gene Symbol` = NA, `All Gene Symbols` = NA)  # Adjust column names as per your ann file
    }
    if (nrow(an18) == 0) {
      an18 <- data.table(annShort = NA, annCurator = NA, annComputer = NA)  # Adjust column names as per your ann2018 file
    }
    # Merge info with sig_DEGs_df row
    result_list[[id]] <- cbind(sig_genes[sig_genes$GeneID == id, ], an15, an18)
  }
}
# Convert the list to a data frame
sig_DEGs_df_ann <- rbindlist(result_list, use.names = TRUE, fill = TRUE)
rangeW<-mean(range(sig_DEGs_df_ann$statW))
sig_DEGs_df_ann$adjusted_W <- abs(sig_DEGs_df_ann$statW - rangeW)
sig_DEGs_df_ann <- sig_DEGs_df_ann[order(sig_DEGs_df_ann$adjusted_W,decreasing = T), ]

write.csv(sig_DEGs_df_ann, file = "results/rankBasedColor_p0.05_4OctFilt_annotated.csv", row.names = F)


###### Visualize
# candidate genes
sig_DEGs_df_ann<-fread("results/rankBasedColor_p0.05_4OctFilt_annotated.csv")

pdf("results/rankBasedColor_p0.05_4OctFilt_annotated.pdf", width = 4, height = 8)
for (id in sig_DEGs_df_ann$GeneID) {
  # Extract normalized counts for the specific gene
  gene_counts <- plotCounts(dds, gene=paste0(id,".v2.1"), intgroup="color", returnData=TRUE)
  # Check if gene_counts is valid
  if (!is.null(gene_counts) && nrow(gene_counts) > 0) {
    # Add the sample IDs from your metadata
    gene_counts$sampleID <- rownames(gene_counts)
    # Create a ggplot with sample IDs as labels
    p <- ggplot(gene_counts, aes(x=color, y=count, label=sampleID)) +
      geom_point(size=5, aes(color=color)) +  # Plot points
      geom_text(vjust=0.5, hjust=-0.5) +  # Add sample ID labels
      scale_color_manual(values=c("pink"="magenta", "white"="grey80")) +  
      labs(title=paste0(id," - ",subset(sig_DEGs_df_ann, sig_DEGs_df_ann$GeneID == id)[,5][[1]]), # 9 if shrinkage not used
           x=sub(";.*", "", subset(sig_DEGs_df_ann, sig_DEGs_df_ann$GeneID == id)[,11][[1]]), # 15 if shrinkage not used
           y="Normalized Counts") +
      theme_minimal() #+ 
    #      scale_y_log10()  # Optionally use log scale for better visualization
    # Print the plot
    print(p)
  }
}
dev.off()


############# ECOTYPIC DIVERGENCE IN ZEP-SUB## FOR JACHYM
##############################################
############# ALTERNATIVE WITH WILCOXON'S TEST
# Load the necessary libraries
library(DESeq2)
library(ggplot2)
library("data.table")

setwd("/home/aa/pigmentation/tkpAnalysis/petalRNASeq/")
# Define function to read STAR count files (reverse strand data)
read_star_counts <- function(filepath) {
  data <- read.table(filepath, header = FALSE, skip = 4)
  counts <- data[, c(1, 4)]  # Extracting gene ID and reverse strand counts (column 4)
  colnames(counts) <- c("GeneID", sub("_ReadsPerGene.out.tab", "", basename(filepath)))
  return(counts)
}
# Define paths to all STAR output files
file_paths <- list.files("/home/aa/pigmentation/tkpAnalysis/petalRNASeq/data/starOutput/", pattern = "_ReadsPerGene.out.tab$", full.names = TRUE)
# Read all count files
count_list <- lapply(file_paths, read_star_counts)
# Merge all count data into a single data frame
counts_data <- Reduce(function(x, y) merge(x, y, by = "GeneID"), count_list)
# Set rownames as Gene IDs and remove GeneID column
rownames(counts_data) <- counts_data$GeneID
counts_data <- counts_data[, -1]

# Load the metadata
metadata <- read.table("data/SampleMetadata.txt", header = TRUE, sep = "\t")
# Ensure that metadata sample IDs match count matrix column names
metadata <- metadata[match(colnames(counts_data), metadata$sampleID), ]
# Check if all samples in the metadata have corresponding counts
if (any(is.na(metadata$sampleID))) {
  stop("Some sample IDs in metadata do not match the read counts.")
}

### Subset only diploids
metadata <- metadata[metadata$ploidy == "2x", ]
# Subset the gene expression data to keep only columns corresponding to 2x individuals
counts_data <- counts_data[, colnames(counts_data) %in% metadata_subs$sampleID]

# Prepare for DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = metadata, design = ~ color) #  design = ~ ploidy + ecotype + color

# Filter
# Keep genes with at least some counts across all samples
hist(rowSums(counts(dds) <= 1),breaks = 30)
hist(log10(rowSums(counts(dds))),breaks = 100)
dds <- dds[rowSums(counts(dds)) >= 140, ]
hist(log10(rowSums(counts(dds))),breaks = 100)
dds <- dds[rowSums(counts(dds) <= 1) <= 8, ]
hist(log10(rowSums(counts(dds))),breaks = 100)

# Plot PCA to check overall sample clustering
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "color")
# Nicer PCA
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "color", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
# Add sample labels from the metadata
pcaData$sampleID <- metadata$sampleID  # Assuming metadata has a sampleID column
# Create a PCA plot with ggplot2, including sample labels
pdf("results/PCAbyColor_2x.pdf", width = 8, height = 5)
ggplot(pcaData, aes(PC1, PC2, color = color, label = sampleID)) +
  geom_point(size = 5) +
  geom_text(vjust = 1.0, hjust = -0.3) +  # Add text labels
  scale_color_manual(values=c("pink"="magenta", "white"="grey70")) +  
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("Gene expression in petals")
dev.off()

# Normalize counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
## Optional: 
write.csv(normalized_counts,"data/normalizedCounts_2x.csv")
hist(log10(rowSums(normalized_counts)),breaks = 100)
pink_samples <- which(metadata$color == "pink")
white_samples <- which(metadata$color == "white")

### Rank-Based Testing Using Wilcoxon Rank-Sum Test ###
# Perform Wilcoxon rank-sum test for each gene
wilcox_test_results <- apply(normalized_counts, 1, function(gene_counts) {
  wilcox.test(gene_counts[pink_samples], gene_counts[white_samples])$p.value
})

wilcox_test_stat <- apply(normalized_counts, 1, function(gene_counts) {
  wilcox.test(gene_counts[pink_samples], gene_counts[white_samples])$statistic
})
# Adjust p-values using Benjamini-Hochberg correction for multiple testing
adjusted_pvals <- p.adjust(wilcox_test_results, method = "BH")
# Combine results into a data frame
wilcox_res <- data.frame(GeneID = substr(rownames(normalized_counts),1,9), pvalue = wilcox_test_results, padj = adjusted_pvals, statW = wilcox_test_stat)
# Filter for significant genes (adjusted p-value < 0.001)
sig_genes <- subset(wilcox_res, padj < 0.05)
sig_genes <- sig_genes[order(sig_genes$padj), ]
# Write significant genes to a file
write.csv(sig_genes, file = "results/rankBasedColor_p0.01_2x.csv", row.names = FALSE)

#### Annotation
# Load data
dict <- fread("/home/aa/Desktop/references/lyrataV2/functions/ALATdict.txt")
dict$AT <- substr(dict$AT, 1, 9)
ann <- fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR10orth_des_20150927.txt", h = T, quote = "")
ann2018 <- fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR11orth_des_20181231.txt", h = T, quote = "")

# Initialize an empty list to collect the rows
result_list <- list()
# Loop through sig_DEGs_df IDs
for (id in sig_genes$GeneID) { # id = "AL2G24170"
  # Find AT code from dict
  if (nrow(subset(dict, dict$AL == id)) > 0) {
    atcode <- subset(dict, dict$AL == id)$AT
    # Extract info from ann and ann2018
    an15 <- subset(ann, ann$`Version-2 A.t. Ortholog` == atcode)[, c(2,4,6,7)]
    an18 <- subset(ann2018, substr(ann2018$AT, 1, 9) == atcode)[1, 4:6]
    # If an15 or an18 are empty, replace with NA
    if (nrow(an15) == 0) {
      an15 <- data.table(`Version-2 A.t. Ortholog` = NA, `Gene Model Description` = NA, `Primary Gene Symbol` = NA, `All Gene Symbols` = NA)  # Adjust column names as per your ann file
    }
    if (nrow(an18) == 0) {
      an18 <- data.table(annShort = NA, annCurator = NA, annComputer = NA)  # Adjust column names as per your ann2018 file
    }
    # Merge info with sig_DEGs_df row
    result_list[[id]] <- cbind(sig_genes[sig_genes$GeneID == id, ], an15, an18)
  }
}
# Convert the list to a data frame
sig_DEGs_df_ann <- rbindlist(result_list, use.names = TRUE, fill = TRUE)
rangeW<-mean(range(sig_DEGs_df_ann$statW))
sig_DEGs_df_ann$adjusted_W <- abs(sig_DEGs_df_ann$statW - rangeW)
sig_DEGs_df_ann <- sig_DEGs_df_ann[order(sig_DEGs_df_ann$adjusted_W,decreasing = T), ]

write.csv(sig_DEGs_df_ann, file = "results/rankBasedColor_p0.05_2x_annotated.csv", row.names = F)


###### Visualize
# candidate genes
sig_DEGs_df_ann<-fread("results/rankBasedColor_p0.01_2x_annotated.csv")

pdf("results/rankBasedColor_p0.01_2x_annotated.pdf", width = 4, height = 8)
for (id in sig_DEGs_df_ann$GeneID) {
  # Extract normalized counts for the specific gene
  gene_counts <- plotCounts(dds, gene=paste0(id,".v2.1"), intgroup="color", returnData=TRUE)
  # Check if gene_counts is valid
  if (!is.null(gene_counts) && nrow(gene_counts) > 0) {
    # Add the sample IDs from your metadata
    gene_counts$sampleID <- rownames(gene_counts)
    # Create a ggplot with sample IDs as labels
    p <- ggplot(gene_counts, aes(x=color, y=count, label=sampleID)) +
      geom_point(size=5, aes(color=color)) +  # Plot points
      geom_text(vjust=0.5, hjust=-0.5) +  # Add sample ID labels
      scale_color_manual(values=c("pink"="magenta", "white"="grey80")) +  
      labs(title=paste0(id," - ",subset(sig_DEGs_df_ann, sig_DEGs_df_ann$GeneID == id)[,5][[1]]), # 9 if shrinkage not used
           x=sub(";.*", "", subset(sig_DEGs_df_ann, sig_DEGs_df_ann$GeneID == id)[,11][[1]]), # 15 if shrinkage not used
           y="Normalized Counts") +
      theme_minimal() #+ 
    #      scale_y_log10()  # Optionally use log scale for better visualization
    # Print the plot
    print(p)
  }
}
dev.off()





############# COLOR DIVERGENCE IN TKO ##
##############################################
############# ALTERNATIVE WITH WILCOXON'S TEST
# Load the necessary libraries
library(DESeq2)
library(ggplot2)
library("data.table")

setwd("/home/aa/pigmentation/tkpAnalysis/petalRNASeq/")
# Define function to read STAR count files (reverse strand data)
read_star_counts <- function(filepath) {
  data <- read.table(filepath, header = FALSE, skip = 4)
  counts <- data[, c(1, 4)]  # Extracting gene ID and reverse strand counts (column 4)
  colnames(counts) <- c("GeneID", sub("_ReadsPerGene.out.tab", "", basename(filepath)))
  return(counts)
}
# Define paths to all STAR output files
file_paths <- list.files("/home/aa/pigmentation/tkpAnalysis/petalRNASeq/data/starOutput/", pattern = "_ReadsPerGene.out.tab$", full.names = TRUE)
# Read all count files
count_list <- lapply(file_paths, read_star_counts)
# Merge all count data into a single data frame
counts_data <- Reduce(function(x, y) merge(x, y, by = "GeneID"), count_list)
# Set rownames as Gene IDs and remove GeneID column
rownames(counts_data) <- counts_data$GeneID
counts_data <- counts_data[, -1]

# Load the metadata
metadata <- read.table("data/SampleMetadata.txt", header = TRUE, sep = "\t")
# Ensure that metadata sample IDs match count matrix column names
metadata <- metadata[match(colnames(counts_data), metadata$sampleID), ]
# Check if all samples in the metadata have corresponding counts
if (any(is.na(metadata$sampleID))) {
  stop("Some sample IDs in metadata do not match the read counts.")
}

### Subset only tetraploids
metadata <- metadata[metadata$ploidy == "4x", ]
# Subset the gene expression data to keep only columns corresponding to 2x individuals
counts_data <- counts_data[, colnames(counts_data) %in% metadata$sampleID]

# Prepare for DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = metadata, design = ~ color) #  design = ~ ploidy + ecotype + color

# Filter
# Keep genes with at least some counts across all samples
hist(rowSums(counts(dds) <= 1),breaks = 30)
hist(log10(rowSums(counts(dds))),breaks = 100)
dds <- dds[rowSums(counts(dds)) >= 100, ]
hist(log10(rowSums(counts(dds))),breaks = 100)
dds <- dds[rowSums(counts(dds) <= 1) <= 6, ]
hist(log10(rowSums(counts(dds))),breaks = 100)

# Plot PCA to check overall sample clustering
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "color")
# Nicer PCA
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "color", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
# Add sample labels from the metadata
pcaData$sampleID <- metadata$sampleID  # Assuming metadata has a sampleID column
# Create a PCA plot with ggplot2, including sample labels
pdf("results/PCAbyColor_4x.pdf", width = 8, height = 5)
ggplot(pcaData, aes(PC1, PC2, color = color, label = sampleID)) +
  geom_point(size = 5) +
  geom_text(vjust = 1.0, hjust = -0.3) +  # Add text labels
  scale_color_manual(values=c("pink"="magenta", "white"="grey70")) +  
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("Gene expression in petals")
dev.off()

# Normalize counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
## Optional: 
write.csv(normalized_counts,"data/normalizedCounts_4x.csv")
hist(log10(rowSums(normalized_counts)),breaks = 100)
pink_samples <- which(metadata$color == "pink")
white_samples <- which(metadata$color == "white")

### Rank-Based Testing Using Wilcoxon Rank-Sum Test ###
# Perform Wilcoxon rank-sum test for each gene
wilcox_test_results <- apply(normalized_counts, 1, function(gene_counts) {
  wilcox.test(gene_counts[pink_samples], gene_counts[white_samples])$p.value
})

wilcox_test_stat <- apply(normalized_counts, 1, function(gene_counts) {
  wilcox.test(gene_counts[pink_samples], gene_counts[white_samples])$statistic
})
# Adjust p-values using Benjamini-Hochberg correction for multiple testing
adjusted_pvals <- p.adjust(wilcox_test_results, method = "BH")
# Combine results into a data frame
wilcox_res <- data.frame(GeneID = substr(rownames(normalized_counts),1,9), pvalue = wilcox_test_results, padj = adjusted_pvals, statW = wilcox_test_stat)
# Filter for significant genes (adjusted p-value < 0.001)
sig_genes <- subset(wilcox_res, pvalue < 0.05)
sig_genes <- sig_genes[order(sig_genes$pvalue), ]
# Write significant genes to a file
write.csv(sig_genes, file = "results/rankBasedColor_pnonadj0.05_4x.csv", row.names = FALSE)

#### Annotation
# Load data
dict <- fread("/home/aa/Desktop/references/lyrataV2/functions/ALATdict.txt")
dict$AT <- substr(dict$AT, 1, 9)
ann <- fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR10orth_des_20150927.txt", h = T, quote = "")
ann2018 <- fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR11orth_des_20181231.txt", h = T, quote = "")

# Initialize an empty list to collect the rows
result_list <- list()
# Loop through sig_DEGs_df IDs
for (id in sig_genes$GeneID) { # id = "AL2G24170"
  # Find AT code from dict
  if (nrow(subset(dict, dict$AL == id)) > 0) {
    atcode <- subset(dict, dict$AL == id)$AT
    # Extract info from ann and ann2018
    an15 <- subset(ann, ann$`Version-2 A.t. Ortholog` == atcode)[, c(2,4,6,7)]
    an18 <- subset(ann2018, substr(ann2018$AT, 1, 9) == atcode)[1, 4:6]
    # If an15 or an18 are empty, replace with NA
    if (nrow(an15) == 0) {
      an15 <- data.table(`Version-2 A.t. Ortholog` = NA, `Gene Model Description` = NA, `Primary Gene Symbol` = NA, `All Gene Symbols` = NA)  # Adjust column names as per your ann file
    }
    if (nrow(an18) == 0) {
      an18 <- data.table(annShort = NA, annCurator = NA, annComputer = NA)  # Adjust column names as per your ann2018 file
    }
    # Merge info with sig_DEGs_df row
    result_list[[id]] <- cbind(sig_genes[sig_genes$GeneID == id, ], an15, an18)
  }
}
# Convert the list to a data frame
sig_DEGs_df_ann <- rbindlist(result_list, use.names = TRUE, fill = TRUE)
rangeW<-mean(range(sig_DEGs_df_ann$statW))
sig_DEGs_df_ann$adjusted_W <- abs(sig_DEGs_df_ann$statW - rangeW)
sig_DEGs_df_ann <- sig_DEGs_df_ann[order(sig_DEGs_df_ann$adjusted_W,decreasing = T), ]

write.csv(sig_DEGs_df_ann, file = "results/rankBasedColor_pnonadj0.05_4x_annotated.csv", row.names = F)


###### Visualize
# candidate genes
sig_DEGs_df_ann<-fread("results/rankBasedColor_pnonadj0.05_4x_annotated.csv")

pdf("results/rankBasedColor_pnonadj0.05_4x_annotated.pdf", width = 4, height = 8)
for (id in sig_DEGs_df_ann$GeneID) {
  # Extract normalized counts for the specific gene
  gene_counts <- plotCounts(dds, gene=paste0(id,".v2.1"), intgroup="color", returnData=TRUE)
  # Check if gene_counts is valid
  if (!is.null(gene_counts) && nrow(gene_counts) > 0) {
    # Add the sample IDs from your metadata
    gene_counts$sampleID <- rownames(gene_counts)
    # Create a ggplot with sample IDs as labels
    p <- ggplot(gene_counts, aes(x=color, y=count, label=sampleID)) +
      geom_point(size=5, aes(color=color)) +  # Plot points
      geom_text(vjust=0.5, hjust=-0.5) +  # Add sample ID labels
      scale_color_manual(values=c("pink"="magenta", "white"="grey80")) +  
      labs(title=paste0(id," - ",subset(sig_DEGs_df_ann, sig_DEGs_df_ann$GeneID == id)[,5][[1]]), # 9 if shrinkage not used
           x=sub(";.*", "", subset(sig_DEGs_df_ann, sig_DEGs_df_ann$GeneID == id)[,11][[1]]), # 15 if shrinkage not used
           y="Normalized Counts") +
      theme_minimal() #+ 
    #      scale_y_log10()  # Optionally use log scale for better visualization
    # Print the plot
    print(p)
  }
}
dev.off()


######### For publication. first, get the dds object above, then run this #############
library(data.table)
library(ggplot2)
library(DESeq2)

sig_DEGs_df_ann <- fread("/home/aa/pigmentation/tkpAnalysis/petalRNASeq/stringTOM60genes.csv")[1:23,]
pdf("results/schematic_boxplots_flavonoidPathway_FINAL_v4.pdf", width = 1.4, height = 1.6)

for (id in sig_DEGs_df_ann$GeneID) {
  gene_id <- paste0(id, ".v2.1")
  gene_counts <- plotCounts(dds, gene = gene_id, intgroup = "color", returnData = TRUE)
  if (!is.null(gene_counts) && nrow(gene_counts) > 0) {
    gene_counts$kcounts <- gene_counts$count / 1000
    ymin <- floor(min(gene_counts$kcounts))
    ymax <- ceiling(max(gene_counts$kcounts))
    
    p <- ggplot(gene_counts, aes(x = color, y = kcounts, fill = color)) +
      geom_boxplot(width = 0.8, outlier.shape = NA, 
                   color = "black", linewidth = 0.6) +
      geom_jitter(width = 0.1, size = 2, shape = 21, stroke = 0.2, fill = "black") +
      scale_fill_manual(values = c("white" = "white", "pink" = "#e3a8ff")) +
      scale_y_continuous(
        limits = c(ymin, ymax),
        breaks = c(ymin, ymax),
        labels = c(ymin, ymax),
        expand = expansion(mult = c(0.05, 0.05))
      ) +
      theme_void() +
      theme(
        axis.text.y = element_text(size = 20),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.margin = unit(rep(0.1, 4), "cm"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 12)
      ) +
      ggtitle(id)
    print(p)
  }
}
dev.off()


##############################################
############# PLOIDY DIFFERENTIATED
# Load the necessary libraries
library(DESeq2)
library(ggplot2)
library("data.table")

setwd("/home/aa/pigmentation/tkpAnalysis/petalRNASeq/")
# Define function to read STAR count files (reverse strand data)
read_star_counts <- function(filepath) {
  data <- read.table(filepath, header = FALSE, skip = 4)
  counts <- data[, c(1, 4)]  # Extracting gene ID and reverse strand counts (column 4)
  colnames(counts) <- c("GeneID", sub("_ReadsPerGene.out.tab", "", basename(filepath)))
  return(counts)
}
# Define paths to all STAR output files
file_paths <- list.files("/home/aa/pigmentation/tkpAnalysis/petalRNASeq/data/starOutput/", pattern = "_ReadsPerGene.out.tab$", full.names = TRUE)
# Read all count files
count_list <- lapply(file_paths, read_star_counts)
# Merge all count data into a single data frame
counts_data <- Reduce(function(x, y) merge(x, y, by = "GeneID"), count_list)
# Set rownames as Gene IDs and remove GeneID column
rownames(counts_data) <- counts_data$GeneID
counts_data <- counts_data[, -1]

# Load the metadata
metadata <- read.table("data/SampleMetadata.txt", header = TRUE, sep = "\t")
# Ensure that metadata sample IDs match count matrix column names
metadata <- metadata[match(colnames(counts_data), metadata$sampleID), ]
# Check if all samples in the metadata have corresponding counts
if (any(is.na(metadata$sampleID))) {
  stop("Some sample IDs in metadata do not match the read counts.")
}

# Prepare for DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = counts_data, colData = metadata, design = ~ ploidy) #  design = ~ ploidy + ecotype + color

# Filter
# Keep genes with at least some counts across all samples
hist(rowSums(counts(dds) == 0),breaks = 30)
hist(log10(rowSums(counts(dds))),breaks = 100)
dds <- dds[rowSums(counts(dds)) >= 240, ]
hist(log10(rowSums(counts(dds))),breaks = 100)
dds <- dds[rowSums(counts(dds) <= 1) <= 15, ]
hist(log10(rowSums(counts(dds))),breaks = 100)

# Plot PCA to check overall sample clustering
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "ploidy")
# Nicer PCA
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "ploidy", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
# Add sample labels from the metadata
pcaData$sampleID <- metadata$sampleID  # Assuming metadata has a sampleID column
# Create a PCA plot with ggplot2, including sample labels
pdf("results/PCAbyPloidy4OctFilt240_15.pdf", width = 8, height = 5)
ggplot(pcaData, aes(PC1, PC2, color = ploidy, label = sampleID)) +
  geom_point(size = 5) +
  geom_text(vjust = 1.0, hjust = -0.3) +  # Add text labels
  scale_color_manual(values=c("2x"="blue", "4x"="red")) +  
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("Gene expression in petals")
dev.off()

# Normalize counts
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)
hist(log10(rowSums(normalized_counts)),breaks = 100)
pink_samples <- which(metadata$ploidy == "4x")
white_samples <- which(metadata$ploidy == "2x")

### Rank-Based Testing Using Wilcoxon Rank-Sum Test ###
# Perform Wilcoxon rank-sum test for each gene
wilcox_test_results <- apply(normalized_counts, 1, function(gene_counts) {
  wilcox.test(gene_counts[pink_samples], gene_counts[white_samples])$p.value
})

wilcox_test_stat <- apply(normalized_counts, 1, function(gene_counts) {
  wilcox.test(gene_counts[pink_samples], gene_counts[white_samples])$statistic
})
# Adjust p-values using Benjamini-Hochberg correction for multiple testing
adjusted_pvals <- p.adjust(wilcox_test_results, method = "BH")
# Combine results into a data frame
wilcox_res <- data.frame(GeneID = substr(rownames(normalized_counts),1,9), pvalue = wilcox_test_results, padj = adjusted_pvals, statW = wilcox_test_stat)
# Filter for significant genes (adjusted p-value < 0.001)
sig_genes <- subset(wilcox_res, padj < 0.001)
sig_genes <- sig_genes[order(sig_genes$padj), ]
# Write significant genes to a file
write.csv(sig_genes, file = "results/rankBasedPloidy_p0.001_4OctFilt.csv", row.names = FALSE)

#### Annotation
# Load data
dict <- fread("/home/aa/Desktop/references/lyrataV2/functions/ALATdict.txt")
dict$AT <- substr(dict$AT, 1, 9)
ann <- fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR10orth_des_20150927.txt", h = T, quote = "")
ann2018 <- fread("/home/aa/Desktop/references/lyrataV2/LyV2_TAIR11orth_des_20181231.txt", h = T, quote = "")

# Initialize an empty list to collect the rows
result_list <- list()
# Loop through sig_DEGs_df IDs
for (id in sig_genes$GeneID) { # id = "AL2G24170"
  # Find AT code from dict
  if (nrow(subset(dict, dict$AL == id)) > 0) {
    atcode <- subset(dict, dict$AL == id)$AT
    # Extract info from ann and ann2018
    an15 <- subset(ann, ann$`Version-2 A.t. Ortholog` == atcode)[, c(2,4,6,7)]
    an18 <- subset(ann2018, substr(ann2018$AT, 1, 9) == atcode)[1, 4:6]
    # If an15 or an18 are empty, replace with NA
    if (nrow(an15) == 0) {
      an15 <- data.table(`Version-2 A.t. Ortholog` = NA, `Gene Model Description` = NA, `Primary Gene Symbol` = NA, `All Gene Symbols` = NA)  # Adjust column names as per your ann file
    }
    if (nrow(an18) == 0) {
      an18 <- data.table(annShort = NA, annCurator = NA, annComputer = NA)  # Adjust column names as per your ann2018 file
    }
    # Merge info with sig_DEGs_df row
    result_list[[id]] <- cbind(sig_genes[sig_genes$GeneID == id, ], an15, an18)
  }
}
# Convert the list to a data frame
sig_DEGs_df_ann <- rbindlist(result_list, use.names = TRUE, fill = TRUE)
rangeW<-mean(range(sig_DEGs_df_ann$statW))
sig_DEGs_df_ann$adjusted_W <- abs(sig_DEGs_df_ann$statW - rangeW)
sig_DEGs_df_ann <- sig_DEGs_df_ann[order(sig_DEGs_df_ann$adjusted_W,decreasing = T), ]

write.csv(sig_DEGs_df_ann, file = "results/rankBasedPloidy_p0.001_4OctFilt_annotated.csv", row.names = F)


###### Visualize
# candidate genes
sig_DEGs_df_ann<-fread("results/rankBasedPloidy_p0.001_4OctFilt_annotated.csv")

pdf("results/rankBasedPloidy_p0.001_4OctFilt_annotated.pdf", width = 4, height = 8)
for (id in sig_DEGs_df_ann$GeneID) {
  # Extract normalized counts for the specific gene
  gene_counts <- plotCounts(dds, gene=paste0(id,".v2.1"), intgroup="ploidy", returnData=TRUE)
  # Check if gene_counts is valid
  if (!is.null(gene_counts) && nrow(gene_counts) > 0) {
    # Add the sample IDs from your metadata
    gene_counts$sampleID <- rownames(gene_counts)
    # Create a ggplot with sample IDs as labels
    p <- ggplot(gene_counts, aes(x=ploidy, y=count, label=sampleID)) +
      geom_point(size=5, aes(color=ploidy)) +  # Plot points
      geom_text(vjust=0.5, hjust=-0.5) +  # Add sample ID labels
      scale_color_manual(values=c("4x"="red", "2x"="blue")) +  
      labs(title=paste0(id," - ",subset(sig_DEGs_df_ann, sig_DEGs_df_ann$GeneID == id)[,5][[1]]), # 9 if shrinkage not used
           x=sub(";.*", "", subset(sig_DEGs_df_ann, sig_DEGs_df_ann$GeneID == id)[,11][[1]]), # 15 if shrinkage not used
           y="Normalized Counts") +
      theme_minimal() #+ 
    #      scale_y_log10()  # Optionally use log scale for better visualization
    # Print the plot
    print(p)
  }
}
dev.off()

