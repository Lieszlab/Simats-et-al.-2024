
###### sc RNA seq : HEART - ccr2tdt mice ####

library(devtools)
library(Seurat)
library(dplyr)
library(Matrix)
library(tidyr)
library(patchwork)

getwd()
setwd("/.../")

test <- Read10X_h5(filename = "/.../filtered_feature_bc_matrix.h5")
table(test["tdtomato", ] ==0) 

heart <- Read10X(data.dir = "/.../filtered_feature_bc_matrix")
heart <- CreateSeuratObject(counts = heart, project = "ccr2tdt-heart", min.cells = 3, min.features = 200)
table(heart$orig.ident) 

heart[["percent.mt"]] <- PercentageFeatureSet(heart, pattern = "^mt-")
head(heart)
sample <- sapply(colnames(GetAssayData(object = heart, slot = "counts")),FUN=function(x){substr(x,18,19)})

sample <- as.numeric(as.character(sample))
names(sample) <- colnames(GetAssayData(object = heart, slot = "counts"))
heart <- AddMetaData(heart, sample, "sample")
table(heart$orig.ident) 

VlnPlot(object = heart, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(heart, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(heart, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

heart_filtered <- subset(heart, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
VlnPlot(object = heart_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

heart_filtered <- NormalizeData(heart_filtered, normalization.method = "LogNormalize", scale.factor = 10000)
heart_filtered <- FindVariableFeatures(heart_filtered, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(heart_filtered), 10)
all.genes <- rownames(heart_filtered)
heart_filtered <- ScaleData(heart_filtered, features = all.genes)
heart_filtered <- RunPCA(heart_filtered, features = VariableFeatures(object = heart_filtered))
DimHeatmap(heart_filtered, dims = 5:15, cells = 500, balanced = TRUE)

heart_filtered <- FindNeighbors(heart_filtered, dims = 1:13)
heart_filtered <- FindClusters(heart_filtered, resolution = 0.4)
heart_filtered <- RunUMAP(heart_filtered, dims = 1:13)

library(RColorBrewer)
mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Set1", n = 8), brewer.pal(name="Accent", n = 6))
DimPlot(heart_filtered, reduction = "umap", label = TRUE, pt.size = 1, label.size = 8)
DimPlot(heart_filtered, reduction = "umap", split.by =  "sample", group.by = "sample")

heartok <- heart_filtered 
new.grouping <- c("group")
heartok[[new.grouping]] <- new.grouping
colnames(heartok@meta.data)

heartok$group[heartok$sample == "2"] <- "stroke"
heartok$group[heartok$sample == "1"] <- "control"
heartok$group[heartok$sample == "3"] <- "control"
heartok$group[heartok$sample == "4"] <- "stroke"
head(heartok@meta.data, 10)

markers <- FindAllMarkers(object = heartok, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)
write.table((markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC)), 'markers_by_cluster.tsv', sep='\t')

table1 <- (markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC))
table <- createWorkbook()
addWorksheet(table, "markersMOMA")
writeData(table, "markersMOMA", table1, startRow = 1, startCol = 1, rowNames = TRUE)
saveWorkbook(table, file = "markers_by_cluster_all_16feb.xlsx", overwrite = TRUE)

for (i in c(0:19)) {
  gmarker.markers <- FindMarkers(heartok, ident.1 = "stroke", ident.2 = "control", group.by = "group",
                                 subset.ident = i, verbose = FALSE, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
  assign(paste0( "markers.", i), gmarker.markers)}

markers.0$p_val_adj2 = p.adjust(markers.0$p_val, method='fdr')# adjust cluster number
wb <- createWorkbook()
addWorksheet(wb, "c0")
writeData(wb, "c0", markers.0, startRow = 1, startCol = 1, rowNames = TRUE)


Idents(heartok) <- "seurat_clusters"
cells <- WhichCells(object = heartok, idents = c("11", "7", "12")) #mono-macrophages : cluster 7, 11, 12
heartMOMA_raw <- subset(heartok, cells = cells)
DimPlot(heartMOMA_raw, split.by = "group")
head(heartMOMA_raw)
heartMOMA <- heartMOMA_raw
heartMOMA <- FindVariableFeatures(heartMOMA, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(heartMOMA), 10)
all.genes <- rownames(heartMOMA)
heartMOMA <- ScaleData(heartMOMA, features = all.genes)
heartMOMA <- RunPCA(heartMOMA, features = VariableFeatures(object = heartMOMA))
DimHeatmap(heartMOMA, dims = 1:10, cells = 500, balanced = TRUE)
heartMOMA <- FindNeighbors(heartMOMA, dims = 1:5) 
heartMOMA <- FindClusters(heartMOMA, resolution = 0.4) 
heartMOMA <- RunUMAP(heartMOMA, dims = 1:5)

markers <- FindAllMarkers(object = heartMOMA, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
table1 <- (markers %>% group_by(cluster) %>% top_n(n = 200, wt = avg_log2FC))
table <- createWorkbook()
addWorksheet(table, "markersMOMA")
writeData(table, "markersMOMA", table1, startRow = 1, startCol = 1, rowNames = TRUE)
saveWorkbook(table, file = "markers_by_cluster_MoMa_16feb.xlsx", overwrite = TRUE)

Idents(heartMOMA) <- "seurat_clusters"
cells <- WhichCells(object = heartMOMA, idents = c("0", "1", "2", "3", "4", "6", "7")) #remove cluster 5
heartMOMA_red <- subset(heartMOMA, cells = cells)
DimPlot(heartMOMA_red, split.by = "group")
head(heartMOMA_red)
heartMOMA_red <- RunUMAP(heartMOMA_red, dims = 1:5)

Idents(object = heartMOMA_red) <- "RNA_snn_res.0.4"
for (i in c(0, 1, 2, 3, 4, 6, 7)) {
  gmarker.markers <- FindMarkers(heartMOMA_red, ident.1 = "stroke", ident.2 = "control", group.by = "group",
                                 subset.ident = i, verbose = FALSE, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1 )
  assign(paste0( "markers.", i), gmarker.markers)}

markers.0$p_val_adj2 = p.adjust(markers.0$p_val, method='fdr') #ajdust for each cluster number
wb <- createWorkbook()
addWorksheet(wb, "c0")
writeData(wb, "c0", markers.0, startRow = 1, startCol = 1, rowNames = TRUE)
saveWorkbook(wb, file = "X.xlsx", overwrite = TRUE)

heartMOMA <- heartMOMA_red
FeaturePlot(heartMOMA, features = c("tdtomato"), pt.size = 2, cols = c("lightgray", "red", "darkred"))
DimPlot(heartMOMA, pt.size = 2, cols = mycolors)
tdt_expression <- GetAssayData(object = heartMOMA, assay = "RNA")["tdtomato",]
head(tdt_expression)
positive_tdt_MOMA = names(which(tdt_expression>0))
negative_tdt_MOMA = names(which(tdt_expression<=0))
pos_tdt_MOMA = subset(heartMOMA, cells=positive_tdt_MOMA)
neg_tdt_MOMA = subset(heartMOMA, cells=negative_tdt_MOMA)

Idents(object = pos_tdt_MOMA) <- "seurat_clusters"
for (i in c(0, 1, 2, 3, 4, 6, 7)) {
  gmarker.markers <- FindMarkers(pos_tdt_MOMA, ident.1 = "stroke", ident.2 = "control", group.by = "group", subset.ident = i, verbose = FALSE,
                                 only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
  assign(paste0( "markers.", i), gmarker.markers)}
markers.0$p_val_adj2 = p.adjust(markers.0$p_val, method='fdr') # ajdust for each cluster number
wb <- createWorkbook()
addWorksheet(wb, "c0")
writeData(wb, "c0", markers.0, startRow = 1, startCol = 1, rowNames = TRUE)
saveWorkbook(wb, file = "X.xlsx", overwrite = TRUE)

table(pos_tdt_MOMA$group)
Idents(pos_tdt_MOMA) <- "group"
DEG_alltdt <- FindMarkers(object = pos_tdt_MOMA, ident.1 = "control", ident.2 = "stroke", group.by = "group",  only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
DEG_alltdt$p_val_adj2 = p.adjust(DEG_alltdt$p_val, method='fdr')
write.table(DEG_alltdt, 'DEG_alltdt.tsv', sep='\t')

# pseudotime: same code as for Bone marrow analysis 

DimPlot(heartok)
heart_noMOMA <- subset(x = heartok, idents = c("0","1", "2", "3", "4", "5", "6", "8", "9", "10", "13", "14", "15", "16", "17", "18", "19"))
DimPlot(heart_noMOMA, label =TRUE)
new.grouping <- c("celltype")
heart_noMOMA[[new.grouping]] <- new.grouping
colnames(heart_noMOMA@meta.data)

heart_noMOMA$celltype[heart_noMOMA$seurat_clusters == "0"] <- "End"
heart_noMOMA$celltype[heart_noMOMA$seurat_clusters == "1"] <- "Fib"
heart_noMOMA$celltype[heart_noMOMA$seurat_clusters == "2"] <- "Fib"
heart_noMOMA$celltype[heart_noMOMA$seurat_clusters == "3"] <- "LymEnd"
heart_noMOMA$celltype[heart_noMOMA$seurat_clusters == "4"] <- "End"
heart_noMOMA$celltype[heart_noMOMA$seurat_clusters == "5"] <- "Fib"
heart_noMOMA$celltype[heart_noMOMA$seurat_clusters == "6"] <- "Bc"
heart_noMOMA$celltype[heart_noMOMA$seurat_clusters == "8"] <- "Tc"
heart_noMOMA$celltype[heart_noMOMA$seurat_clusters == "9"] <- "Per"
heart_noMOMA$celltype[heart_noMOMA$seurat_clusters == "10"] <- "Neu"
heart_noMOMA$celltype[heart_noMOMA$seurat_clusters == "13"] <- "NKT"
heart_noMOMA$celltype[heart_noMOMA$seurat_clusters == "14"] <- "Epi"
heart_noMOMA$celltype[heart_noMOMA$seurat_clusters == "15"] <- "ProEnd"
heart_noMOMA$celltype[heart_noMOMA$seurat_clusters == "16"] <- "Tc"
heart_noMOMA$celltype[heart_noMOMA$seurat_clusters == "17"] <- "Fib"
heart_noMOMA$celltype[heart_noMOMA$seurat_clusters == "18"] <- "SW"
heart_noMOMA$celltype[heart_noMOMA$seurat_clusters == "19"] <- "Fib"

head(pos_tdt_MOMA)
DimPlot(pos_tdt_MOMA, label =TRUE)
new.grouping <- c("celltype")
heartMOMA[[new.grouping]] <- new.grouping
colnames(heartMOMA@meta.data)

pos_tdt_MOMA$celltype[pos_tdt_MOMA$seurat_clusters == "0"] <- "MoL"
pos_tdt_MOMA$celltype[pos_tdt_MOMA$seurat_clusters == "1"] <- "TRM1"
pos_tdt_MOMA$celltype[pos_tdt_MOMA$seurat_clusters == "2"] <- "TRM2"
pos_tdt_MOMA$celltype[pos_tdt_MOMA$seurat_clusters == "3"] <- "MoMa"
pos_tdt_MOMA$celltype[pos_tdt_MOMA$seurat_clusters == "4"] <- "TRM3"
pos_tdt_MOMA$celltype[pos_tdt_MOMA$seurat_clusters == "6"] <- "MoH"
pos_tdt_MOMA$celltype[pos_tdt_MOMA$seurat_clusters == "7"] <- "MoMa2"

ccr2_merged <- merge(heart_noMOMA, pos_tdt_MOMA)
new.grouping <- c("celltype3")
ccr2_merged[[new.grouping]] <- new.grouping
colnames(ccr2_merged@meta.data)
tdt_expression <- GetAssayData(object = ccr2_merged, assay = "RNA")["tdtomato",]

ccr2_merged$celltype3[ccr2_merged$celltype == "End"] <- "End"
ccr2_merged$celltype3[ccr2_merged$celltype == "Fib"] <- "Fib"
ccr2_merged$celltype3[ccr2_merged$celltype == "LymEnd"] <- "LymEnd"
ccr2_merged$celltype3[ccr2_merged$celltype == "Bc"] <- "Bc"
ccr2_merged$celltype3[ccr2_merged$celltype == "Tc"] <- "Tc"
ccr2_merged$celltype3[ccr2_merged$celltype == "Per"] <- "Per"
ccr2_merged$celltype3[ccr2_merged$celltype == "Neu"] <- "Neu"
ccr2_merged$celltype3[ccr2_merged$celltype == "NKT"] <- "NKT"
ccr2_merged$celltype3[ccr2_merged$celltype == "Epi"] <- "Epi"
ccr2_merged$celltype3[ccr2_merged$celltype == "ProEnd"] <- "ProEnd"
ccr2_merged$celltype3[ccr2_merged$celltype == "MoL" & tdt_expression>0] <- "MoL_pos"
ccr2_merged$celltype3[ccr2_merged$celltype == "TRM1" & tdt_expression>0] <- "TRM1_pos"
ccr2_merged$celltype3[ccr2_merged$celltype == "TRM2" & tdt_expression>0] <- "TRM2_pos"
ccr2_merged$celltype3[ccr2_merged$celltype == "MoMa" & tdt_expression>0] <- "MoMa_pos"
ccr2_merged$celltype3[ccr2_merged$celltype == "TRM3" & tdt_expression>0] <- "TRM3_pos"
ccr2_merged$celltype3[ccr2_merged$celltype == "MoH" & tdt_expression>0] <- "MoH_pos"
ccr2_merged$celltype3[ccr2_merged$celltype == "MoMa2" & tdt_expression>0] <- "MoMa2_pos"

# only STROKE group #
Idents(object = ccr2_merged_new) <-  "group"
ccr2_merged_stroke <- subset(x = ccr2_merged_new, idents = "stroke")

# only CONTROL group #
Idents(object = ccr2_merged_new) <-  "group"
ccr2_merged_control <- subset(x = ccr2_merged_new, idents = "control")


# Code was obatined from CellChat vignette:
# https://htmlpreview.github.io/?https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.html
# Citation can be found in the manuscript

