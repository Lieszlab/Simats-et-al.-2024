---
title: "Heart bulk sequencing _ Simats et al., 2023"

---



################################################################################
# Load libraries and sources
################################################################################

library(magrittr)
library(tidyverse)
library(DESeq2)

################################################################################
# Sample filtering
################################################################################

## # Expression matrix provided is already filtered. The filtering process is
## # howewever shared below:

## # Calculate the percent of mitochondrial genes per cell
## NoERCC <- PercentageFeatureSet(NoERCC,
##                                features = grep(pattern = "^MT-", x = rownames(x = NoERCC), value = TRUE),
##                                col.name = "percent.mt")

## # Calculate the percent of ribosomal genes per cell
## NoERCC <- PercentageFeatureSet(NoERCC,
##                                features = grep(pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA", x = rownames(x = NoERCC@assays$RNA@data), value = TRUE),
##                                col.name = "percent.rb")

## # Filter the mitochondrial high samples
## NoERCC_filtered <- subset(NoERCC, subset = percent.mt < 20)

## # Filter the ribosomal high samples
## NoERCC_filtered <- subset(NoERCC_filtered, subset = percent.rb < 8)

################################################################################
# Create DESeq2 object
################################################################################

# Load the provided expression matrix
NoERCC_filtered <- read.table("expression_matrix.tsv", header = TRUE, dec = ".", sep = "\t")

# Create DESeq2 object from expression matrix
pre_dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(NoERCC_filtered),
  colData = data.frame(
    "sample" = colnames(NoERCC_filtered),
    "condition" = c("stroke", "control", "stroke", "stroke", "stroke", "control", "stroke", "stroke", "stroke", "stroke", "control", "control", "control")
  ),
  design = ~ condition)

################################################################################
# Gene filtering
################################################################################

# Explore the overall gene expression
quantile(rowSums(counts(pre_dds)), seq(0, 1, by = 0.1))

threshold <- quantile(rowSums(counts(pre_dds)), seq(0, 1, length.out = length(rowSums(counts(pre_dds)))))
threshold_df <- data.frame("quant" = rev(seq(0, 1, length.out = length(threshold))), threshold)

# quantile plot of all genes by total sum
ggplot(threshold_df, aes(quant, threshold)) +
  theme_bw() +
  geom_point(size = 1) +
  scale_x_continuous(minor_breaks = seq(0, 1, by = 0.01), breaks = seq(0, 1, by = 0.05))
# same thing as a boxplot in log10
boxplot(log10(rowSums(counts(pre_dds)+1))~as.factor(rep(1, nrow(counts(pre_dds)))))
# same with histogram in log10 scale
hist(log10(rowSums(counts(pre_dds))+1))

# zoom in the quantile plot to the head up to the first 5%
ggplot(threshold_df %>% filter(quant <= 0.05), aes(quant, threshold)) +
  theme_bw() +
  geom_point(size = 1) +
  scale_x_continuous(minor_breaks = seq(0, 1, by = 0.001), breaks = seq(0, 1, by = 0.005)) +
  geom_hline(yintercept = 20000, col = "red")

sum(rowSums(counts(pre_dds)) > 60000)
# 20 genes are expressed above 60000 counts
sum(rowSums(counts(pre_dds)) > 20000)
# 74 genes are expressed above 20000 counts

rownames(pre_dds)[rowSums(counts(pre_dds)) > 60000]
rownames(pre_dds)[rowSums(counts(pre_dds)) > 20000]

# filtering very high expressed genes perturbating the normalization
pre_dds <- pre_dds[rowSums(counts(pre_dds)) < 20000, ]

pre_dds <- estimateSizeFactors(pre_dds)
plotSparsity(pre_dds)
# still a lot of single genes expressed in only one sample

# Check the distribution of gene counts per sample
quantile(counts(pre_dds), seq(0, 1, by = 0.1))
# 70% of the data is below a count of 1

# Filter in only the genes that are expressed in 3 samples and with >1 counts
pre_dds <- pre_dds[rowSums(counts(pre_dds) > 1 ) >= 3, ]
pre_dds
# we have 8319 genes remaining

# let's check the diversity
pre_dds <- estimateSizeFactors(pre_dds)
plotSparsity(pre_dds)
# highly expressed genes are still unique to several samples. Maybe raise the
# min sample threshold.
table(colData(pre_dds)$condition)
# the smallest group is control with 5 samples. Let's set up the minimum to 5.

pre_dds <- pre_dds[rowSums(counts(pre_dds) > 1 ) >= 5, ]
pre_dds
# we have 3424 genes remaining

# let's check the diversity
pre_dds <- estimateSizeFactors(pre_dds)
plotSparsity(pre_dds)

################################################################################
# Differential gene expression analysis
################################################################################

pre_dds <- DESeq(pre_dds)

resultsNames(pre_dds)

res_condition <- results(pre_dds, contrast = c("condition", "stroke", "control"))
(res_condition_filtered <- subset(res_condition, padj <= 0.05))
dim(res_condition_filtered)

# 15 genes are differentially expressed between stroke and control
write.csv(res_condition_filtered, "Data/res_condition_strokevscontrol_filtered.csv")

################################################################################
# Volcano
################################################################################

library(EnhancedVolcano)

EnhancedVolcano(res_condition,
    lab = rownames(res_condition),
    x = 'log2FoldChange',
    y = 'padj',
    pCutoff = 0.05,
    drawConnectors = TRUE)

################################################################################
# Normalisation
################################################################################

# Normalisation
vsd <- vst(pre_dds, blind=TRUE)
rld <- rlog(pre_dds, blind=TRUE)
vst <- varianceStabilizingTransformation(pre_dds, blind=TRUE)

# Gene variation per sample before normalisation
boxplot(assay(pre_dds), las=2, main="cts")

# Check if no over inflation of low expressed (vsd)
mean_vsd <- apply(assay(vsd), 1, mean)
sd_vsd <- apply(assay(vsd), 1, sd)
plot(mean_vsd, sd_vsd, asp = 1)

# Gene variation per sample after normalisation with vst
boxplot(assay(vsd), las=2, main="vsd")

# Check if no over inflation of low expressed (rld)
mean_rld <- apply(assay(rld), 1, mean)
sd_rld <- apply(assay(rld), 1, sd)
plot(mean_rld, sd_rld, asp = 1)

# Gene variation per sample after normalisation with rld
boxplot(assay(rld), las=2, main="rld")

# Check if no over inflation of low expressed (vst)
mean_vst <- apply(assay(vst), 1, mean)
sd_vst <- apply(assay(vst), 1, sd)
plot(mean_vst, sd_vst, asp = 1)

# Gene variation per sample after normalisation with vst
boxplot(assay(vst), las=2, main="vst")

# we choose rld normalisation but still remaining highly variable genes
# (plotSparsity was pointing to those already)

################################################################################
# Heatmap of DGE
################################################################################

library(pheatmap)

my_sample_col <- data.frame(condition = colData(rld)$condition)
row.names(my_sample_col) <- colnames(rld)

my_gene_row <- data.frame(regulation = c("up", "up", "down", "up", "up", "up", "down", "up", "down", "up", "up", "down", "down", "down", "down"))
rownames(my_gene_row) <- rownames(res_condition_filtered)

my_data <- assay(rld[rownames(rld) %in% rownames(res_condition_filtered), ])
my_data %<>% as.data.frame()
my_data <- select(my_data, paste0("S", c("02", "07", "15", "19", "20", "01", "03", "04", "05", "09", "11", "12", "13")))
my_data %<>% as.matrix()

dim(my_data)
my_datasorted <- my_data[order(my_gene_row$regulation), ]

my_datasorted_zscore <- t(my_datasorted) %>% scale() %>% t()

pheatmap(my_datasorted_zscore, cluster_rows = F, cluster_cols = F,
         show_rownames = T, show_colnames = T,
         fontsize = 10,
         annotation_col = my_sample_col,
         annotation_row = my_gene_row,
         color = colorRampPalette(c(rep("midnightblue", 1), rep("blue2", 1), "grey", rep("orange", 3), rep("red", 3)))(256))

################################################################################
# Boxplot of DGE
################################################################################

my_sample_col <- data.frame(condition = colData(rld)$condition)
row.names(my_sample_col) <- colnames(rld)
my_sample_col %<>% rownames_to_column("sample")

my_data_gg <- assay(rld[rownames(rld) %in% rownames(res_condition_filtered), ]) %>%
  as.data.frame() %>%
  select(., paste0("S", c("02", "07", "15", "19", "20", "01", "03", "04", "05", "09", "11", "12", "13"))) %>%
  rownames_to_column("gene") %>%
  gather("sample", "expression", -gene) %>%
  left_join(., my_sample_col, by = "sample")

ggplot(my_data_gg, aes(gene, expression, col = condition)) +
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_color_manual(values = c("black", "red"))
