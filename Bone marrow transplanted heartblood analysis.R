

###### sc RNA seq : BMT 2 HEART-BLOOD ####

library(devtools)
library(Seurat)
library(dplyr)
library(Matrix)
library(tidyr)
library(patchwork)

getwd()
setwd("/...")

test <- Read10X_h5(filename = "/.../filtered_feature_bc_matrix.h5")
table(test$`Gene Expression`["agfp", ] ==0) 
table(test$`Gene Expression`["kgfp", ] ==0) 

getwd()
matrix_dir = "/.../filtered_feature_bc_matrix/"
list.files(matrix_dir)

barcode.path <- paste0(matrix_dir, "barcodes.tsv.gz")
features.path <- paste0(matrix_dir, "features.tsv.gz")
matrix.path <- paste0(matrix_dir, "matrix.mtx.gz")
mat <- readMM(file = matrix.path)

barcode.names = read.delim(barcode.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)
feature.names = read.delim(features.path, 
                           header = FALSE,
                           stringsAsFactors = FALSE)

GEX <- mat[0:32287,]
GEX_feature.names <- feature.names[0:32287,]

HTO <- mat[32288:32289,]
HTO_feature.names <- feature.names[32288:32289,]
HTO_feature.names
colnames(GEX) = barcode.names$V1
colnames(HTO) = barcode.names$V1
rownames(GEX) = GEX_feature.names$V2
rownames(HTO) = HTO_feature.names$V1
rownames(GEX) <- make.unique(rownames(GEX))
rownames(HTO) <- paste0(rownames(HTO), "-HTO")

# Create seurat object 
s1 <- CreateSeuratObject (counts = GEX, project = "BMT2")
s1[["HTO"]] <- CreateAssayObject(counts = HTO)
s1 <- NormalizeData(s1, assay = "HTO", normalization.method = "CLR")
s1 <- HTODemux(s1, assay = "HTO", positive.quantile = 0.99)
table(s1$HTO_classification.global)
Idents(s1) <- "HTO_maxID"
RidgePlot(s1, assay = "HTO", features = rownames(s1[["HTO"]])[1:4], ncol = 2)
Idents(s1) <- "HTO_classification.global"
VlnPlot(s1, features = "nCount_RNA", pt.size = 0.1, log = TRUE)
BMT <- subset(s1, idents = "Doublet", invert = TRUE)

sample <- sapply(colnames(GetAssayData(object = BMT, slot = "counts")),
                 FUN=function(x){substr(x,18,18)})
sample <- as.numeric(as.character(sample))
names(sample) <- colnames(GetAssayData(object = BMT, slot = "counts"))
BMT <- AddMetaData(BMT, sample, "sample")

new.grouping <- c("sample2")
BMT[[new.grouping]] <- new.grouping
BMT$sample2[BMT$HTO_classification == "B0302-HTO" ] <- "Blood"
BMT$sample2[BMT$HTO_classification == "B0303-HTO" ] <- "Heart"

new.grouping <- c("group")
BMT[[new.grouping]] <- new.grouping
BMT$group[BMT$sample == "1" ] <- "Control"
BMT$group[BMT$sample == "2" ] <- "Stroke"
BMT$group[BMT$sample == "3" ] <- "Control"
BMT$group[BMT$sample == "4" ] <- "Stroke"

BMT <- subset(BMT, idents = "Negative", invert = TRUE)
BMT[["percent.mt"]] <- PercentageFeatureSet(BMT, pattern = "^mt-")
VlnPlot(object = BMT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
BMT <- subset(BMT, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
VlnPlot(object = BMT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

BMT <- NormalizeData(BMT, normalization.method = "LogNormalize", scale.factor = 10000)
BMT <- FindVariableFeatures(BMT, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(BMT), 10)
all.genes <- rownames(BMT)
BMT <- ScaleData(BMT, features = all.genes)
BMT <- RunPCA(BMT, features = VariableFeatures(object = BMT))
DimHeatmap(BMT, dims = 5:20, cells = 500, balanced = TRUE)

BMT <- FindNeighbors(BMT, dims = 1:13) 
BMT <- FindClusters(BMT, resolution = 0.8) 
BMT <- RunUMAP(BMT, dims = 1:13)

FeaturePlot(BMT, reduction = 'umap', features = c("agfp", "kgfp"), blend = TRUE)
GFP <- list(c("agfp", "kgfp"))
BMT <- AddModuleScore(object = BMT, features = GFP, name = "GFP_score")
FeaturePlot(object = BMT, features = "GFP_score1", cols = c("lightgray", "darkgreen"))

new.grouping <- c("gfp_signal")
BMT[[new.grouping]] <- new.grouping
colnames(BMT@meta.data)
BMT$gfp_signal[BMT$GFP_score1 >0 ] <- "positive"
BMT$gfp_signal[BMT$GFP_score1 <=0] <- "negative"
table(BMT$gfp_signal) 
Idents(BMT) <- "gfp_signal"
BMT_subset_gfp <- WhichCells(object = BMT, idents = c("positive"))
BMT_subset_gfp <- subset(x = BMT, cells = BMT_subset_gfp)
DimPlot(BMT_subset_gfp, group.by = "seurat_clusters", split.by = "group", label = TRUE)

Idents(BMT_subset_gfp) <- "seurat_clusters"
BMT_subset_gfp_moma <- WhichCells(object = BMT_subset_gfp, idents = c("1", "0", "3", "13", "10", "12", "4", "8"))
BMT_subset_gfp_moma <- subset(x = BMT_subset_gfp, cells = BMT_subset_gfp_moma)

Idents(BMT_subset_gfp_moma) <- "seurat_clusters"
BMTmarkers <- FindAllMarkers(object = BMT_subset_gfp_moma, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
BMTmarkers <- BMTmarkers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_log2FC)
table <- createWorkbook()
addWorksheet(table, "markers_BMT")
writeData(table, "markers_BMT", BMTmarkers, startRow = 1, startCol = 1, rowNames = TRUE)
saveWorkbook(table, file = "markers_by_cluster_X.xlsx", overwrite = TRUE)

Idents(BMT_subset_gfp_moma) <- "sample2"
Blood <- FindMarkers(object = BMT_subset_gfp_moma, ident.1 = "Control", ident.2 = "Stroke", group.by = "group", subset.ident = "Blood", only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
Heart <- FindMarkers(object = BMT_subset_gfp_moma, ident.1 = "Control", ident.2 = "Stroke", group.by = "group", subset.ident = "Heart", only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(Blood, 'Blood.tsv', sep='\t')
write.table(Heart, 'Heart.tsv', sep='\t')
Blood$p_val_adj2 = p.adjust(Blood$p_val, method='fdr')
Heart$p_val_adj2 = p.adjust(Heart$p_val, method='fdr')
Blood.degs <- row.names(Blood[which(Blood$p_val_adj2 <0.05),])
Heart.degs <- row.names(Heart[which(Heart$p_val_adj2 <0.05),])
xlist <- list(
  A = Blood.degs, 
  B = Heart.degs )

ggvenn(
  xlist, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4, 
  auto_scale = TRUE
) 


# Pseudotime analysis (Monocle3 package): same code as in the Bone Marrow analysis

