---
title: "BM myeloid cells _ Simats et al., 2023"

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

BM.data <- CreateSeuratObject (counts = data, min.cells = 3, min.features = 200, project = "Macrophage_ALL_A")

BM.data[['percent.mito']] <- percent.mito
VlnPlot(object = BM.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
FeatureScatter(object = BM.data, feature1 = "nCount_RNA", feature2 = "percent.mito")
FeatureScatter(object = BM.data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

BM.data <- subset(x = BM.data, subset = nFeature_RNA > 500 & percent.mito < '0.07')
VlnPlot(object = BM.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
BM.data

BM.data <- NormalizeData(object = BM.data, normalization.method = "LogNormalize", scale.factor = 1e4)
BM.data <- FindVariableFeatures(object = BM.data, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
length(x = VariableFeatures(object = BM.data))

BM.data <- RunPCA(object = BM.data, features = VariableFeatures(object = BM.data), verbose = FALSE)

# Get batches based on cell names
samples_batches <- sapply(colnames(GetAssayData(object = BM.data, slot = "counts")),
                          FUN=function(x){substr(x,18,1)})
samples_batches <- as.numeric(as.character(samples_batches))
names(samples_batches) <- colnames(GetAssayData(object = BM.data, slot = "counts"))

complexity.per.cell <- apply(GetAssayData(object = BM.data, slot = "counts"),
                             2, function(x) sum(x>0))
plot(complexity.per.cell ~ jitter(samples_batches,2))      

sample.effect <- samples_batches
BM.data <- AddMetaData(BM.data, sample.effect, "sample.effect")

DimHeatmap(object = BM.data, dims = 1:18, cells = 500, balanced = TRUE)
BM.data <- JackStraw(object = BM.data, num.replicate = 100)
BM.data <- ScoreJackStraw(object = BM.data, dims = 1:20)

BM.data <- FindNeighbors(object = BM.data, dims = 1:12)
BM.data <- FindClusters(object = BM.data, resolution = 0.8)
BM.data <- RunTSNE(object = BM.data, dims = 1:12)

BM.data <- RunUMAP(object = BM.data, dims = 1:12)
DimPlot(object = BM.data, reduction = 'umap', label = TRUE)

# Identification of clusters
BM.markers <- FindAllMarkers(object = BM.data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
BM.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

# Remove T/B cells
FeaturePlot(BM.data, reduction = 'umap', features = c("Cd3e", "Cd79a"))
BM.ok <- WhichCells(object = BM.data, idents = c("0","1","2","5","6","7","8","9","10","11","12","13","14","16","18","20","21","22","23"))
BM.ok_cluster <- subset(x = BM.data, cells = BM.ok)
DimPlot(object = BM.ok_cluster, reduction = 'tsne', label = TRUE, pt.size = 0.5)

dataBM <- NormalizeData(object = BM.ok_cluster, normalization.method = "LogNormalize", scale.factor = 1e4)

dataBM <- FindVariableFeatures(object = dataBM, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
length(x = VariableFeatures(object = dataBM))
dataBM <- ScaleData(object = dataBM, features = rownames(x = dataBM), vars.to.regress = c("nCount_RNA", "percent.mito"))
dataBM <- RunPCA(object = dataBM, features = VariableFeatures(object = dataBM), verbose = FALSE)

ElbowPlot(object = dataBM)
dataBM <- FindNeighbors(object = dataBM, dims = 1:13)
dataBM <- FindClusters(object = dataBM, resolution = 0.8)

dataBM <- RunTSNE(object = dataBM, dims = 1:13)
options(repr.plot.width=13, repr.plot.height=13)

dataBM <- RunUMAP(object = dataBM, dims = 1:13)
DimPlot(dataBM, reduction = 'umap', label = FALSE)
        
# Identification of new clusters
BM.markers.ok <- FindAllMarkers(object = dataBM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
BM.markers.ok %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

# Feature plots
progenitors <- FeaturePlot(dataBM, reduction = 'umap', features = c("Mpo", "Ms4a3", "Cd34", "Csf1r", "Flt3", "Itgb2l"))


# DEG between stroke and control
new.grouping <- c("group")
dataBM[[new.grouping]] <- new.grouping
colnames(dataBM@meta.data)

dataBM$group[dataBM$sample.effect == "1"] <- "Naive"
dataBM$group[dataBM$sample.effect == "2"] <- "Stroke"
dataBM$group[dataBM$sample.effect == "3"] <- "Naive"
dataBM$group[dataBM$sample.effect == "4"] <- "Stroke"
dataBM$group[dataBM$sample.effect == "5"] <- "Naive"
dataBM$group[dataBM$sample.effect == "6"] <- "Stroke"
dataBM$group[dataBM$sample.effect == "7"] <- "Naive"
dataBM$group[dataBM$sample.effect == "8"] <- "Stroke"
head(dataBM@meta.data, 10000)

#For each cluster:
Idents(object = dataBM) <- "seurat_clusters"
BM.c0_ok <- FindMarkers(object = dataBM, ident.1 = "Naive", ident.2 = "Stroke", group.by = "group", subset.ident = "0", only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(BM.c0_ok, 'BM.c0.ok.tsv', sep='\t')

#For group of clusters: Monohigh, Monolow, Neutro, Baso, pDCs, Mek, HSPC
Idents(object = dataBM) <- "seurat_clusters"
BM.monohigh <- FindMarkers(object = dataBM, ident.1 = "Naive", ident.2 = "Stroke", group.by = "group", subset.ident = c("10", "0", "11", "6"), only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(BM.monohigh, 'BM.monohigh.tsv', sep='\t')

Idents(object = dataBM) <- "seurat_clusters"
BM.monolow <- FindMarkers(object = dataBM, ident.1 = "Naive", ident.2 = "Stroke", group.by = "group", subset.ident = c("16"), only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(BM.monolow, 'BM.monolow.tsv', sep='\t')

Idents(object = dataBM) <- "seurat_clusters"
BM.neutro <- FindMarkers(object = dataBM, ident.1 = "Naive", ident.2 = "Stroke", group.by = "group", subset.ident = c("2", "13", "1", "15"), only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(BM.neutro, 'BM.neutro.tsv', sep='\t')

Idents(object = dataBM) <- "seurat_clusters"
BM.baso <- FindMarkers(object = dataBM, ident.1 = "Naive", ident.2 = "Stroke", group.by = "group", subset.ident = c("17", "20", "14"), only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(BM.baso, 'BM.baso.tsv', sep='\t')

Idents(object = dataBM) <- "seurat_clusters"
BM.pDCs <- FindMarkers(object = dataBM, ident.1 = "Naive", ident.2 = "Stroke", group.by = "group", subset.ident = c("3", "19"), only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(BM.pDCs, 'BM.pDCs.tsv', sep='\t')

Idents(object = dataBM) <- "seurat_clusters"
BM.Mek <- FindMarkers(object = dataBM, ident.1 = "Naive", ident.2 = "Stroke", group.by = "group", subset.ident = c("12"), only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(BM.Mek, 'BM.Mek.tsv', sep='\t')

Idents(object = dataBM) <- "seurat_clusters"
BM.HSPC <- FindMarkers(object = dataBM, ident.1 = "Naive", ident.2 = "Stroke", group.by = "group", subset.ident = c("4", "9", "5", "18"), only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(BM.HSPC, 'BM.HSPC.tsv', sep='\t')



#Pseudotime trajectory analysis: Monocle3
library(monocle3)

#Remove NK cluster of cells 
dataBM
Idents(dataBM) <- "seurat_clusters"
dataBM <- subset(dataBM, idents = '7', invert = TRUE)
csf <- dataBM
unique(csf@meta.data$seurat_clusters)


#monocle v3, for all cells or stroke/control separately):
gene_annotation <- as.data.frame(rownames(csf@reductions[["pca"]]@feature.loadings),
                                 row.names = rownames(csf@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"


cell_metadata <- as.data.frame(csf@assays[["RNA"]]@counts@Dimnames[[2]],
                               row.names = csf@assays[["RNA"]]@counts@Dimnames[[2]])
colnames(cell_metadata) <- "barcode"


New_matrix <- csf@assays[["RNA"]]@counts
New_matrix <- New_matrix[rownames(csf@reductions[["pca"]]@feature.loadings), ]
expression_matrix <- New_matrix

cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata,
                                     gene_metadata = gene_annotation)


recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)

cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

list_cluster <- csf@active.ident
names(list_cluster) <- csf@assays[["RNA"]]@data@Dimnames[[2]]

cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <-csf@reductions[["umap"]]@cell.embeddings
cds_from_seurat@preprocess_aux$gene_loadings <- csf@reductions[["pca"]]@feature.loadings

cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = F)

plot_cells(cds_from_seurat, 
           color_cells_by = 'cluster',
           label_groups_by_cluster=TRUE,
           label_leaves=T,
           label_branch_points=TRUE,
           graph_label_size=4)

Idents(csf) <- "seurat_clusters"
cds_from_seurat <- order_cells(cds_from_seurat, root_cells = WhichCells(csf, idents = 9))

plot_cells(cds_from_seurat, 
           color_cells_by = 'pseudotime',
           label_cell_groups = FALSE,
           label_leaves = FALSE,
           label_branch_points = FALSE,
           cell_size = 0.7,
           graph_label_size=0)

csf <- AddMetaData(
        object = csf,
        metadata = cds_from_seurat@principal_graph_aux@listData$UMAP$pseudotime,
        col.name = "Pseudotime")


#Calculate PCA and Euclidean distances: same as for Peripheral organs


