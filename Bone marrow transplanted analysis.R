---
  title: "BM myeloid cells - BM transplant experiment _ Simats et al., 2023"

---
  
library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(umap)
library(plyr)
library(factoextra)      

matrix_dir = "X:/.../filtered_feature_bc_matrix/"
list.files(matrix_dir)

data <- Read10X(data.dir = matrix_dir)

BMtransplant.data <- CreateSeuratObject (counts = data, min.cells = 3, min.features = 200, project = "BM_transplant")

BMtransplant.data[['percent.mito']] <- percent.mito
VlnPlot(object = BMtransplant.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
FeatureScatter(object = BMtransplant.data, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object = BMtransplant.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

BMtransplant.data <- subset(x = BMtransplant.data, subset = nFeature_RNA > 500 & percent.mito < '0.07')
VlnPlot(object = BMtransplant.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
BMtransplant.data

BMtransplant.data <- NormalizeData(object = BMtransplant.data, normalization.method = "LogNormalize", scale.factor = 1e4)
BMtransplant.data <- FindVariableFeatures(object = BMtransplant.data, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
length(x = VariableFeatures(object = BMtransplant.data))

BMtransplant.data <- RunPCA(object = BMtransplant.data, features = VariableFeatures(object = BM.data), verbose = FALSE)

# Get batches based on cell names
samples_batches <- sapply(colnames(GetAssayData(object = BMtransplant.data, slot = "counts")),
                          FUN=function(x){substr(x,18,1)})
samples_batches <- as.numeric(as.character(samples_batches))
names(samples_batches) <- colnames(GetAssayData(object = BMtransplant.data, slot = "counts"))

complexity.per.cell <- apply(GetAssayData(object = BMtransplant.data, slot = "counts"),
                             2, function(x) sum(x>0))
plot(complexity.per.cell ~ jitter(samples_batches,2))      

sample.effect <- samples_batches
BMtransplant.data <- AddMetaData(BMtransplant.data, sample.effect, "sample.effect")

DimHeatmap(object = BMtransplant.data, dims = 1:18, cells = 500, balanced = TRUE)
BMtransplant.data <- JackStraw(object = BMtransplant.data, num.replicate = 100)
BMtransplant.data <- ScoreJackStraw(object = BMtransplant.data, dims = 1:20)

BMtransplant.data <- FindNeighbors(object = BMtransplant.data, dims = 1:9)
BMtransplant.data <- FindClusters(object = BMtransplant.data, resolution = 0.8)
BMtransplant.data <- RunTSNE(object = BMtransplant.data, dims = 1:9)

BMtransplant.data <- RunUMAP(object = BMtransplant.data, dims = 1:9)
DimPlot(object = BMtransplant.data, reduction = 'umap', label = TRUE)

# Identification of clusters
BM.markers <- FindAllMarkers(object = BMtransplant.data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
BM.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

# Remove T/B cells
FeaturePlot(BMtransplant.data, reduction = 'umap', features = c("Cd3e", "Cd79a"))
BM.ok <- WhichCells(object = BMtransplant.data, idents = c("0","1","2","4","5","8","10","11","12","13","14","15","17","18","20","22","23"))
BMtransplant.data.ok <- subset(x = BMtransplant.data, cells = BM.ok)
DimPlot(object = BMtransplant.data.ok, reduction = 'tsne', label = TRUE, pt.size = 0.5)


BMtransplant.data.ok <- NormalizeData(object = BMtransplant.data.ok, normalization.method = "LogNormalize", scale.factor = 1e4)

BMtransplant.data.ok <- FindVariableFeatures(object = BMtransplant.data.ok, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
length(x = VariableFeatures(object = BMtransplant.data.ok))
BMtransplant.data.ok <- ScaleData(object = BMtransplant.data.ok, features = rownames(x = BMtransplant.data.ok), vars.to.regress = c("nCount_RNA", "percent.mito"))
BMtransplant.data.ok <- RunPCA(object = BMtransplant.data.ok, features = VariableFeatures(object = BMtransplant.data.ok), verbose = FALSE)

ElbowPlot(object = BMtransplant.data.ok)
BMtransplant.data.ok <- FindNeighbors(object = BMtransplant.data.ok, dims = 1:9)
BMtransplant.data.ok <- FindClusters(object = BMtransplant.data.ok, resolution = 0.8)

BMtransplant.data.ok <- RunTSNE(object = BMtransplant.data.ok, dims = 1:9)
options(repr.plot.width=13, repr.plot.height=13)

BMtransplant.data.ok <- RunUMAP(object = BMtransplant.data.ok, dims = 1:9)
DimPlot(BMtransplant.data.ok, reduction = 'umap', label = FALSE)

# Remove more B cells
BMtransplant.ok <- WhichCells(object = BMtransplant.data.ok, idents = c("0","1","2","3","4","5","6","7","8","9","10","11","12","14","16","17","18","19","20","21","22"))
BMtransplant.data.ok <- subset(x = BMtransplant.data.ok, cells = BMtransplant.ok)

BMtransplant.data.ok <- NormalizeData(object = BMtransplant.data.ok, normalization.method = "LogNormalize", scale.factor = 1e4)

BMtransplant.data.ok <- FindVariableFeatures(object = BMtransplant.data.ok, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
length(x = VariableFeatures(object = BMtransplant.data.ok))
BMtransplant.data.ok <- ScaleData(object = BMtransplant.data.ok, features = rownames(x = BMtransplant.data.ok), vars.to.regress = c("nCount_RNA", "percent.mito"))
BMtransplant.data.ok <- RunPCA(object = BMtransplant.data.ok, features = VariableFeatures(object = BMtransplant.data.ok), verbose = FALSE)

BMtransplant.data.ok <- FindNeighbors(object = BMtransplant.data.ok, dims = 1:9)
BMtransplant.data.ok <- FindClusters(object = BMtransplant.data.ok, resolution = 0.8)
BMtransplant.data.ok <- RunTSNE(object = BMtransplant.data.ok, dims = 1:9)
options(repr.plot.width=13, repr.plot.height=13)
BMtransplant.data.ok <- RunUMAP(object = BMtransplant.data.ok, dims = 1:9)

# Identification of new clusters
BM.markers.ok <- FindAllMarkers(object = BMtransplant.data.ok, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
BM.markers.ok %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

# Feature plots
progenitors <- FeaturePlot(BMtransplant.data.ok, reduction = 'umap', features = c("Mpo", "Ms4a3", "Cd34", "Csf1r", "Flt3", "Itgb2l"))


# DEG between stroke and control
new.grouping <- c("group")
BMtransplant.data.ok[[new.grouping]] <- new.grouping
colnames(BMtransplant.data.ok@meta.data)

BMtransplant.data.ok$group[BMtransplant.data.ok$sample.effect == "1"] <- "Control"
BMtransplant.data.ok$group[BMtransplant.data.ok$sample.effect == "2"] <- "Stroke"
BMtransplant.data.ok$group[BMtransplant.data.ok$sample.effect == "3"] <- "Stroke"
BMtransplant.data.ok$group[BMtransplant.data.ok$sample.effect == "4"] <- "Control"

# Select only cells that are eGFP positive
egfp_expression <- GetAssayData(object = csf, assay = "RNA")["EGFP",]
head(egfp_expression)
positive_egfp_ids = names(which(egfp_expression>0))

BMtransplant.eGFP = subset(csf, cells=positive_egfp_ids)
saveRDS(BMtransplant.eGFP, file = "BMtransplant.eGFP.rds")

# Identification of new clusters
features.pos.eGFP <- rownames(BMtransplant.eGFP)
markers.remove <- grep(pattern = "^Rpl|^Rps", x = rownames(BMtransplant.eGFP), value = TRUE) # remove ribosomal genes from the list of DEG
features.pos.eGFP <- features.pos.eGFP[!(features.pos.eGFP%in%markers.remove)]
pos.egfp.markers.norib <- FindAllMarkers(object = BMtransplant.eGFP, features = features.pos.eGFP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pos.egfp.markers.norib %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

# Identify DEG per CLUSTER (only cluster 0 is showed)
Idents(object = BMtransplant.eGFP) <- "seurat_clusters"
egfp_0 <- FindMarkers(object = BMtransplant.eGFP, features = features.pos.eGFP, ident.1 = "Stroke", ident.2 = "Naive", group.by = "group", subset.ident = "0", only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(egfp_0, 'egfp.pos_0.cvs', sep='\t')


# Monocle3 pseudotime analysis --> same code as for Bone Marrow analysis
# Calculate PCA and Euclidean distances: same code as for Peripheral organs analysis


