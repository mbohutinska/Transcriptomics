##############################
############ RNA SEQ #########
##############################

#Install DESeq2 if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
install.packages("DBI")
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
counts_data <- counts_data[, colnames(counts_data) %in% metadata$sampleID]

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
pcaData <- plotPCA(vsd, intgroup = "ecotype", returnData = TRUE)
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
sig_genes <- subset(wilcox_res, padj < 0.009)
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
sig_DEGs_df_ann<-sig_DEGs_df_ann[1:100,]
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






#################################
############## WCGNA ############
#################################
#install.packages("WGCNA")
setwd("/home/aa/pigmentation/tkpAnalysis/petalRNASeq/WCGNA/")
library(WGCNA)
library(dplyr)
library(tidyr)
library(readr)
library(pheatmap)

options(stringsAsFactors = FALSE);
allowWGCNAThreads() 

metadata<-read.table("../data/SampleMetadata.txt",h=T)
dataExpr<-read.csv("../data/normalizedCounts_filtered_240RC_15nonzero_7Oct.csv")
rownames(dataExpr) <- dataExpr$X
dataExpr <- dataExpr[, -1]

#### Optional, to subsample genes based on DE results
wilcox<-read.csv("/home/aa/pigmentation/tkpAnalysis/petalRNASeq/results/rankBasedColor_all_20601_4OctFilt.csv",h=T)
wilcoxSubs<-subset(wilcox,wilcox$padj<0.1)
dataExpr<-subset(dataExpr, rownames(dataExpr) %in% paste0(wilcoxSubs$GeneID,".v2.1"))
#### The end of optional

dataExpr <- as.data.frame(t(dataExpr))
# My are already filtered, so just in case
# gsg <- goodSamplesGenes(dataExpr, verbose = 3)
# dataExpr <- dataExpr[gsg$goodSamples, gsg$goodGenes]

# Perform sample clustering to detect outliers and visualize them:
sampleTree <- hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample Clustering", xlab="", sub="")


#### Optional
# Choose a soft-thresholding power (e.g., 6) based on scale-free topology fit. 
#This runs long and is only needed once. for my full data it's 8.
powers <- c(1:30)
# alternative from maize tutorial: powers = c(c(1:20), seq(from = 22, to=30, by=2))
sft <- pickSoftThreshold(dataExpr, powerVector = powers, verbose = 5)
power <- sft$powerEstimate
# Vizualise and pick the threshold
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit, signed R^2",
     main = paste("Scale independence"))
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.90, col = "red")
plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers,
     cex = cex1, col = "red")

######### The result was 8, no need to repeat ###############

# Construct the network and detect modules:
cor <- WGCNA::cor
net <- blockwiseModules(dataExpr, power = 8, TOMType = "unsigned",
                        minModuleSize = 20, maxBlockSize = 3000,
                        mergeCutHeight = 0.25, numericLabels = TRUE, 
                        pamRespectsDendro = FALSE, saveTOMs = TRUE, 
                        verbose = 3, loadTOM = F) # Or power = power, Try with min modulesize 20, reassign = 0.05, # IF TOMs GENERATED, LOAD TOM = T
cor<- stats::cor
sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)
pdf(file = "module_tree_blockwise.pdf", width = 8, height = 6);
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# First, ensure sample IDs in metadata are formatted similarly to the row names in dataExpr
metadata$sampleID <- paste0("X", as.numeric(metadata$sampleID))
# Reorder metadata based on the row names of dataExpr
metadata_reordered <- metadata[match(rownames(dataExpr), metadata$sampleID), ]
# Check if the reordering worked correctly
all.equal(rownames(dataExpr), metadata_reordered$sampleID)
# relate to colors
rownames(dataExpr)
metadata_reordered$color
# Ensure the categorical traits are converted to binary (dummy) variables
metadata_reordered_binary <- model.matrix(~ ploidy + ecotype + color, data = metadata_reordered)[, -1]
rownames(metadata_reordered_binary)<-metadata_reordered$sampleID

# Calculate correlations and p-values
MEs <- moduleEigengenes(dataExpr, net$colors)$eigengenes
correlationMatrix <- cor(MEs, metadata_reordered_binary, use = "pairwise.complete.obs")
correlationMatrixAbs<-as.data.frame(abs(correlationMatrix))
pValues <- corPvalueStudent(correlationMatrix, nrow(MEs))

# Using pheatmap to plot the correlation matrix with significance annotations
pdf("traitModuleCorAbs.pdf",height = 12,width = 6)
pheatmap(correlationMatrixAbs, cluster_rows = TRUE, cluster_cols = TRUE,display_numbers = signif(pValues, 2), number_format = "%.2f",main = "Module-Trait Relationships")
dev.off()
# Explore the results
net$colors[net$blockGenes[[1]]]
MEs
correlationMatrixAbs
pValues
colorModuleVector<-substr(rownames(subset(correlationMatrixAbs,correlationMatrixAbs$colorwhite > 0.7)),3,4)
colorModuleGenes <- net$colors[net$colors == colorModuleVector]
# Export correlations
module_df <- data.frame(
  gene_id = names(net$colors),
  modules = net$colors,
  colors = labels2colors(net$colors)
)
write.csv(module_df,"geneModules.csv")
write.csv(MEs,"modulesIndividuals.csv")
write.csv(cbind(correlationMatrix,pValues),"moduleTraitCorr.csv")

### Identify hub genes and prepare Cytoscape input
#modules_of_interest=colorModuleVector #OR: c("0","1","2")
genes_of_interest = module_df %>%
  subset(modules %in% colorModuleVector)
expr_of_interest = dataExpr[,genes_of_interest$gene_id]
expr_of_interest[1:5,1:5]
TOM = TOMsimilarityFromExpr(expr_of_interest,power = 8)
row.names(TOM) = colnames(expr_of_interest)
colnames(TOM) = colnames(expr_of_interest)
TOM[1:15,1:15]

#zoom into a candidate gene
geneCand<-TOM["AL2G24170.v2.1",]
geneCand<-geneCand[geneCand<1]
hist(geneCand,breaks = 100)
summary(geneCand)

geneCand<-TOM["AL7G25360.v2.1",]
geneCand<-geneCand[geneCand<1]
hist(geneCand,breaks = 100)
summary(geneCand)
geneCand1<-subset(geneCand,geneCand > 0.2)




### Filter network irrespective of gene of interest and find the most connected hub genes
# Set the threshold
threshold <- 0.11
# Filter the TOM matrix based on the threshold
filtered_edges <- as.data.frame(which(TOM > threshold, arr.ind = TRUE))
filtered_edges <- filtered_edges %>%
  mutate(
    correlation = TOM[cbind(row, col)],
    gene1 = row.names(TOM)[row],
    gene2 = colnames(TOM)[col]
  ) %>%
  select(-row, -col) %>%
  filter(gene1 != gene2) %>%
  mutate(
    module1 = module_df[gene1, ]$modules,
    module2 = module_df[gene2, ]$modules
  )
# Remove duplicate edges by ensuring each pair appears only once
filtered_edges <- filtered_edges %>%
  unique()

connectedness<-as.data.frame(table(filtered_edges$gene1))
hist(connectedness$Freq,breaks = 100)
hist(filtered_edges$correlation,breaks = 100)
summary(filtered_edges$correlation)
# Write the filtered edge list to a file
write_delim(filtered_edges, file = "filtered_edgelist_R2_0.11.tsv", delim = "\t")


##############The same but log transformed ######### WORKS BETTER
#################################
########## WGCNA ###############
#################################
setwd("/home/aa/pigmentation/tkpAnalysis/petalRNASeq/WCGNA/")
library(WGCNA)
library(dplyr)
library(data.table)
library(pheatmap)
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

# Load metadata and expression
metadata <- read.table("../data/SampleMetadata.txt", h = TRUE)
dataExpr <- read.csv("../data/normalizedCounts_filtered_240RC_15nonzero_7Oct.csv")
rownames(dataExpr) <- dataExpr$X
dataExpr <- dataExpr[, -1]

# Optional: filter to DE genes (padj < 0.1)
wilcox <- read.csv("/home/aa/pigmentation/tkpAnalysis/petalRNASeq/results/rankBasedColor_all_20601_4OctFilt.csv", h = TRUE)
wilcoxSubs <- subset(wilcox, wilcox$padj < 0.05)
dataExpr <- subset(dataExpr, rownames(dataExpr) %in% paste0(c(wilcoxSubs$GeneID,"AL8G18230",'AL6G19370','AL4G25370'), ".v2.1"))

## alternative for candidate gene set:
# dataExpr <- subset(dataExpr, rownames(dataExpr) %in% paste0(c('AL2G24170','AL1G13000','AL5G19690','AL7G30640','AL1G13540','AL8G29160','AL7G42320','AL7G25360','AL6G43450','AL4G33370','AL6G28270','AL4G33640','AL1G20780','AL6G24580','AL8G29370','AL6G48240','AL7G47960','AL4G27010','AL6G24080','AL6G18330','AL2G19520','AL3G35670','AL6G41210','AL1G37490','AL8G41090','AL4G29100','AL2G19700','AL7G31440','AL6G53230','AL3G23970','AL4G47060','AL5G31640','AL8G18230'), ".v2.1"))


# Transpose and log-transform
dataExpr <- as.data.frame(t(dataExpr))
dataExpr <- log2(dataExpr + 1)

# Optional: sample clustering
sampleTree <- hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample Clustering", xlab="", sub="")

cor <- WGCNA::cor
net <- blockwiseModules(dataExpr, power = 8, TOMType = "unsigned",
                        minModuleSize = 20, maxBlockSize = 3000,
                        mergeCutHeight = 0.25, numericLabels = TRUE, 
                        pamRespectsDendro = FALSE, saveTOMs = TRUE, 
                        verbose = 3, loadTOM = F)

mergedColors <- labels2colors(net$colors)
pdf("module_tree_blockwise.pdf", width = 8, height = 6)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors", dendroLabels = FALSE, hang = 0.03)
dev.off()

# Traits and module correlations
metadata$sampleID <- paste0("X", as.numeric(metadata$sampleID))
metadata_reordered <- metadata[match(rownames(dataExpr), metadata$sampleID), ]
metadata_reordered_binary <- model.matrix(~ ploidy + ecotype + color, data = metadata_reordered)[, -1]
rownames(metadata_reordered_binary) <- metadata_reordered$sampleID

MEs <- moduleEigengenes(dataExpr, net$colors)$eigengenes
correlationMatrix <- cor(MEs, metadata_reordered_binary, use = "pairwise.complete.obs")
correlationMatrixAbs <- abs(correlationMatrix)
pValues <- corPvalueStudent(correlationMatrix, nrow(MEs))

pdf("traitModuleCorAbs.pdf", height = 12, width = 6)
pheatmap(correlationMatrixAbs, cluster_rows = TRUE, cluster_cols = TRUE,
         display_numbers = signif(pValues, 2), number_format = "%.2f",
         main = "Module-Trait Relationships")
dev.off()

# Extract module and expression info
correlationMatrixAbs<-as.data.frame(correlationMatrixAbs)
colorModuleVector <- substr(rownames(subset(correlationMatrixAbs, correlationMatrixAbs$colorwhite > 0.3)), 3, 4)
module_df <- data.frame(
  gene_id = names(net$colors),
  modules = net$colors,
  colors = labels2colors(net$colors)
)
write.csv(module_df, "geneModules.csv")
write.csv(MEs, "modulesIndividuals.csv")
write.csv(cbind(correlationMatrix, pValues), "moduleTraitCorr.csv")

# Subset TOM for genes of interest
genes_of_interest <- module_df %>% filter(modules %in% colorModuleVector)
expr_of_interest <- dataExpr[, genes_of_interest$gene_id]
TOM <- TOMsimilarityFromExpr(expr_of_interest, power = 8)
rownames(TOM) <- colnames(expr_of_interest)
colnames(TOM) <- colnames(expr_of_interest)

geneCand <- TOM["AL2G24170.v2.1",]
TOM["AL2G24170.v2.1","AL4G25370.v2.1"]
TOM["AL2G24170.v2.1","AL6G19370.v2.1"]


geneCand <- geneCand[geneCand < 1]
hist(geneCand,breaks = 100) # this value is not correlation but rather neighborhood similarity. It can be much lower..
quantile(geneCand,0.9) # a lot of random posibilities for filtration criteria
threshold <- 0.075 # same here..
coexpressed_genes <- names(geneCand[geneCand > threshold])
selected_genes <- c("AL2G24170.v2.1", coexpressed_genes)
TOM_sub <- TOM[selected_genes, selected_genes]

# Network visualization
library(igraph)
TOM_sub[TOM_sub < threshold] <- 0
graph_obj <- graph_from_adjacency_matrix(TOM_sub, mode = "undirected", weighted = TRUE, diag = FALSE)
module_colors <- module_df$colors[match(V(graph_obj)$name, module_df$gene_id)]
V(graph_obj)$color <- module_colors

pdf("coexpressedWithPAP2.pdf", width = 14, height = 14)
plot(graph_obj, vertex.label = substr(V(graph_obj)$name, 1, 9),
     vertex.size = 5, edge.width = E(graph_obj)$weight * 10,
     layout = layout_with_fr)
dev.off()

# Annotate pre-selected PAP2-co-expressed genes
# Load ortholog dictionary
dict <- fread("/home/aa/Desktop/references/lyrataV2/functions/ALATdict.txt")
dict$AT <- substr(dict$AT, 1, 9)
# Get unique genes from TOM_sub
genes_to_annotate <- substr(colnames(TOM_sub), 1, 9)
genes_to_annotate <- unique(genes_to_annotate)

annotation_list <- list()
pap2_id <- "AL2G24170.v2.1"
for (id in genes_to_annotate) {
  full_id <- colnames(TOM_sub)[substr(colnames(TOM_sub), 1, 9) == id][1]
  tom_value <- if (!is.na(full_id) && full_id != pap2_id) TOM_sub[pap2_id, full_id] else NA
  match_row <- dict[AL == id]
  if (nrow(match_row) > 0) {
    atcode <- match_row$AT[1]
    annotation_list[[id]] <- data.table(gene_id = id, TAIR_ID = atcode, TOM_with_PAP2 = tom_value)
  } else {
    annotation_list[[id]] <- data.table(gene_id = id, TAIR_ID = NA, TOM_with_PAP2 = tom_value)
  }
}

# Combine and save
TOM_annotation <- rbindlist(annotation_list, use.names = TRUE, fill = TRUE)
write.csv(TOM_annotation, "annotated_TOM_sub_geneslogtr.csv", row.names = FALSE)

#### Add PAP2-TOM to the table of DEGs
tab<-fread("/home/aa/pigmentation/tkpAnalysis/petalRNASeq/results/rankBasedColor4OctFilt_0.01_annotated.csv")
pap2_id <- "AL2G24170.v2.1"
tab$TOM_with_PAP2 <- NA

for (i in 1:nrow(tab)) { #  i=1
  gene_id <- paste0(tab$GeneID[i],".v2.1")
  if (gene_id %in% colnames(TOM) && pap2_id %in% rownames(TOM)) {
    tab$TOM_with_PAP2[i] <- TOM[pap2_id, gene_id]
  }
}
fwrite(tab, "rankBasedColor4OctFilt_0.01_annotated_PAP2TOM.csv")
tab1<-subset(tab, tab$TOM_with_PAP2!=1)
plot(tab1$adjusted_W~tab1$TOM_with_PAP2)
############### THE END

