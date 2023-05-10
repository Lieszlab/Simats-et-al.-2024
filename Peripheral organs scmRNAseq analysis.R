---
title: "Peripheral organs _ Simats et al., 2023 "

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

GEX <- mat[-31054:-31057,]
GEX_feature.names <- feature.names[-31054:-31057,]

HTO <- mat[31054:31057,]
HTO_feature.names <- feature.names[31054:31057,]
HTO_feature.names

colnames(GEX) = barcode.names$V1
colnames(HTO) = barcode.names$V1
rownames(GEX) = GEX_feature.names$V2
rownames(HTO) = HTO_feature.names$V1

rownames(GEX) <- make.unique(rownames(GEX))
rownames(HTO) <- paste0(rownames(HTO), "-HTO")

#Read Read HTO data for TotalSeq-A0305

#for all samples (D1-D8):
HTO_dir = "X:/.../umi_count/"
list.files(HTO_dir)

barcode.HTO.path <- paste0(HTO_dir, "barcodes.tsv.gz")
features.HTO.path <- paste0(HTO_dir, "features.tsv.gz")
matrix.HTO.path <- paste0(HTO_dir, "matrix.mtx.gz")
HTO.1 <- readMM(file = matrix.HTO.path)
HTO.1 <- HTO.1[1,]

barcode.HTO.1.names = read.delim(barcode.HTO.path, 
                                 header = FALSE,
                                 stringsAsFactors = FALSE)
feature.HTO.names = read.delim(features.HTO.path, 
                               header = FALSE,
                               stringsAsFactors = FALSE)


HTO.joined <- c(HTO.1, HTO.2, HTO.3, HTO.4, HTO.5, HTO.6, HTO.7, HTO.8)
feature.HTO.names <- feature.HTO.names[1,]
feature.HTO.names

barcode.HTO.names <- Reduce(function(x,y) merge(x,y,1,all=T,sort=F),list(barcode.HTO.1.names, barcode.HTO.2.names, barcode.HTO.3.names, barcode.HTO.4.names, barcode.HTO.5.names, barcode.HTO.6.names, barcode.HTO.7.names, barcode.HTO.8.names))

df <- as.matrix(HTO.joined)
df <- t(df)
df <- as.data.frame(df)
df

row.names(df) <- 1
df

mat <- as.sparse(df)
mat <- as(mat, "dgTMatrix")

colnames(mat) = barcode.HTO.names$V1
rownames(mat) = feature.HTO.names

# Create seurat Object #
peripheral.data <- CreateSeuratObject (counts = GEX, project = "Macrophage_ALL_C")
peripheral.data[["HTO"]] <- CreateAssayObject(counts = HTO)

peripheral.data <- NormalizeData(peripheral.data, assay = "HTO", normalization.method = "CLR")
peripheral.data <- HTODemux(peripheral.data, assay = "HTO", positive.quantile = 0.99)
table(peripheral.data$HTO_classification.global)


Idents(peripheral.data) <- "HTO_classification.global"
VlnPlot(peripheral.data, features = "nCount_RNA", pt.size = 0.1, log = TRUE)


peripheral.data <- subset(peripheral.data, idents = "Doublet", invert = TRUE)
peripheral.data

mito.features <- grep(pattern = "^mt-", x = rownames(x = peripheral.data), value = TRUE)
percent.mito <- Matrix::colSums(x = GetAssayData(object = peripheral.data, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = peripheral.data, slot = 'counts'))

peripheral.data[['percent.mito']] <- percent.mito
VlnPlot(object = peripheral.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)

peripheral.data <- subset(x = peripheral.data, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mito < "0.07")
VlnPlot(object = peripheral.data, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
peripheral.data

peripheral.data <- NormalizeData(object = peripheral.data, normalization.method = "LogNormalize", scale.factor = 1e4)

peripheral.data <- FindVariableFeatures(object = peripheral.data, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
length(x = VariableFeatures(object = peripheral.data))

peripheral.data <- ScaleData(object = peripheral.data, features = rownames(x = peripheral.data), vars.to.regress = c("nCount_RNA", "percent.mito"))


peripheral.data <- RunPCA(object = peripheral.data, features = VariableFeatures(object = peripheral.data), verbose = FALSE)
print(x = peripheral.data[['pca']], dims = 1:5, nfeatures = 5, projected = FALSE)

# Create batch effect
peripheral.data

samples_batches <- sapply(colnames(GetAssayData(object = peripheral.data, slot = "counts")),
                          FUN=function(x){substr(x,18,19)})

samples_batches <- as.numeric(as.factor(samples_batches))
names(samples_batches) <- colnames(GetAssayData(object = peripheral.data, slot = "counts"))

# Complexity per technical batch
complexity.per.cell <- apply(GetAssayData(object = peripheral.data, slot = "counts"),
                             2, function(x) sum(x>0))
plot(complexity.per.cell ~ jitter(samples_batches,2))

sample.effect <- samples_batches
peripheral.data <- AddMetaData(peripheral.data, sample.effect, "sample.effect")

peripheral.data <- ProjectDim(object = peripheral.data)
DimHeatmap(object = peripheral.data, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(object = peripheral.data)

peripheral.data <- FindNeighbors(object = peripheral.data, dims = 1:12)
peripheral.data <- FindClusters(object = peripheral.data, resolution = 0.8)

peripheral.data <- RunTSNE(object = peripheral.data, dims = 1:12)
Idents(object = peripheral.data) <- "RNA_snn_res.0.8"
DimPlot(object = peripheral.data, reduction = 'tsne', label = TRUE)

saveRDS(peripheral.data, file = "Macrophage_ALL_C.rds")

#Assignation of TotalSeq-A0305 - assay="SeqA"
rownames(x = peripheral.data[["SeqA"]])

md <- peripheral.data@meta.data
md$A0305.count <- GetAssayData(object = peripheral.data, assay = "SeqA", slot = "data")["A0305-CTTTGTCTTTGTGAG", ]
head(md)

aggregate(A0305.count ~ hash.ID, md, mean)
aggregate(A0305.count ~ hash.ID, md, sd)

for (i in 1:dim(md)[1]){
        if (md$hash.ID[i] == "Negative" & md$A0305.count[i] >= 2){
                md$HTO_final[i] <- "A0305-SeqA"
        } else {md$HTO_final[i] <- md$hash.ID[i]}
}

head(md)
unique(md$HTO_final)
count(md, "HTO_final")

peripheral.data <- AddMetaData(object=peripheral.data, metadata=md)
head(peripheral.data@meta.data)
DimPlot(object = peripheral.data, reduction = 'tsne', group.by= "HTO_final", label = TRUE)

saveRDS(peripheral.data, file = "Macrophage_ALL_C.rds")

______________

peripheral.data <- readRDS("Macrophage_ALL_C.rds")
colnames(peripheral.data@meta.data)

peripheral.data <- RunUMAP(object = peripheral.data, dims = 1:13)
DimPlot(object = peripheral.data, reduction = 'umap', group.by= "HTO_final", label = TRUE)

# Asign organ name
new.grouping <- c("organ")
peripheral.data[[new.grouping]] <- new.grouping
colnames(peripheral.data@meta.data)

peripheral.data$organ[peripheral.data$HTO_final == "A0305-SeqA"] <- "Spleen"
peripheral.data$organ[peripheral.data$HTO_final == "B0301-HTO"] <- "Heart"
peripheral.data$organ[peripheral.data$HTO_final == "B0302-HTO"] <- "Blood"
peripheral.data$organ[peripheral.data$HTO_final == "B0303-HTO"] <- "Liver"
peripheral.data$organ[peripheral.data$HTO_final == "B0304-HTO"] <- "Lung"
peripheral.data$organ[peripheral.data$HTO_final == "Negative"] <- "Negative"
head(peripheral.data@meta.data, 10000)

# Remove NEGATIVE cells
Idents(object = peripheral.data) <- "organ"
remove.data <- WhichCells(object = peripheral.data, idents = c("Negative"), invert = TRUE)
peripheral.data <- subset(x = peripheral.data, cells = remove.data)
DimPlot(object = peripheral.data, reduction = 'tsne', label = TRUE, pt.size = 1)

# Remove T/B cells
FeaturePlot(peripheral.data, reduction = 'tsne', features = c("Cd3e", "Cd19", "Cd79a"))
Idents(object = peripheral.data) <- "seurat_clusters"
peripheral.ok <- WhichCells(object = peripheral.data, idents = c("0", "1", "2", "3", "4", "5","6", "7", "8","10", "11", "13", "14","15", "16", "17", "18", "20", "21"))
peripheral.data <- subset(x = peripheral.data, cells = peripheral.ok)
DimPlot(object = peripheral.data, reduction = 'tsne', label = TRUE, pt.size = 0.5)


dataPERIPHERAL <- NormalizeData(object = peripheral.data, normalization.method = "LogNormalize", scale.factor = 1e4)

dataPERIPHERAL <- FindVariableFeatures(object = dataPERIPHERAL, selection.method = 'mean.var.plot', mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf))
length(x = VariableFeatures(object = dataPERIPHERAL))
dataPERIPHERAL <- ScaleData(object = dataPERIPHERAL, features = rownames(x = dataPERIPHERAL), vars.to.regress = c("nCount_RNA", "percent.mito"))
dataPERIPHERAL <- RunPCA(object = dataPERIPHERAL, features = VariableFeatures(object = dataPERIPHERAL), verbose = FALSE)

dataPERIPHERAL <- FindNeighbors(object = dataPERIPHERAL, dims = 1:12)
dataPERIPHERAL <- FindClusters(object = dataPERIPHERAL, resolution = 0.8)
dataPERIPHERAL <- RunTSNE(object = dataPERIPHERAL, dims = 1:12)
dataPERIPHERAL <- RunUMAP(object = dataPERIPHERAL, dims = 1:12)
DimPlot(dataPERIPHERAL, reduction = 'umap', label = FALSE, cols = col_vector)
write.table(dataPERIPHERAL@meta.data, 'meta.data.peripheral.tsv', sep='\t')

# Cluster identification
peripheral.markers <- FindAllMarkers(object = dataPERIPHERAL, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
peripheral.markers %>% group_by(cluster) %>% top_n(n = 2,wt = avg_log2FC)

# Downsample to 1000 for each specific organ
Idents(dataPERIPHERAL) <- "organ"
dataPERIPHERAL_heart <- WhichCells(object = dataPERIPHERAL, idents = c("Heart"))
dataPERIPHERAL_heart <- subset(x = dataPERIPHERAL, cells = dataPERIPHERAL_heart)
subset_blood <- subset(x = dataPERIPHERAL_heart, downsample = 1000)

# DEG between stroke and control conditions
featuresperipheralALL <- rownames(peripheral.markers)
markers.remove.2 <- grep(pattern = "^Rpl|^Rps|^mt", x = rownames(peripheral.markers), value = TRUE) ## remove ribosomal & mitocondrial genes from the list of genes
featuresperipheralALL <- featuresperipheralALL[!(featuresperipheralALL%in%markers.remove.2)]

# For each cluster:
Idents(object = dataPERIPHERAL) <- "RNA_snn_res.0.8"
Peripheral.c0_ok <- FindMarkers(object = dataPERIPHERAL, features = featuresperipheralALL, ident.1 = "Stroke", ident.2 = "Naive", group.by = "group", subset.ident = "0", only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(Peripheral.c0_ok, 'Peripheral.c0_ok.tsv', sep='\t')


# Isolation of the two clusters of monocytes/macrophages 
Peripheral.mono <- WhichCells(object = dataPERIPHERAL, idents = c("0", "12"))
Peripheral.mono <- subset(x = dataPERIPHERAL, cells = Peripheral.mono)
DimPlot(object = Peripheral.mono, reduction = 'umap', label = FALSE, pt.size = 1, split.by = "group")

Peripheral.mono0and12 <- FindMarkers(object = Peripheral.mono, ident.1 = "Stroke", ident.2 = "Naive", group.by = "group", only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.1)
write.table(Peripheral.mono0and12, 'Peripheral.mono0and12.tsv', sep='\t')

EnhancedVolcano(Peripheral.mono0and12,
                lab = rownames(Peripheral.mono0and12),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 0.20,
                pointSize = 3.0,
                labSize = 6.0,
                title = "Peripheral.mono0and12",
                xlim = c(-1,1),
                ylim = c(0, 30),
                col=c('gray', 'gray', 'gray', 'red3'),
                colAlpha = 1,
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                arrowheads = FALSE,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                max.overlaps = 10)


#Calculate PCA and Euclidean distances

head(dataPERIPHERAL@meta.data)

new.grouping <- c("condition")
dataPERIPHERAL[[new.grouping]] <- new.grouping
colnames(dataPERIPHERAL@meta.data)

dataPERIPHERAL$condition[dataPERIPHERAL$organ == "Heart" & dataPERIPHERAL$group == "Naive"] <- "Heart_Naive"
dataPERIPHERAL$condition[dataPERIPHERAL$organ == "Heart" & dataPERIPHERAL$group == "Stroke"] <- "Heart_Stroke"
dataPERIPHERAL$condition[dataPERIPHERAL$organ == "Lung" & dataPERIPHERAL$group == "Naive"] <- "Lung_Naive"
dataPERIPHERAL$condition[dataPERIPHERAL$organ == "Lung" & dataPERIPHERAL$group == "Stroke"] <- "Lung_Stroke"
dataPERIPHERAL$condition[dataPERIPHERAL$organ == "Liver" & dataPERIPHERAL$group == "Naive"] <- "Liver_Naive"
dataPERIPHERAL$condition[dataPERIPHERAL$organ == "Liver" & dataPERIPHERAL$group == "Stroke"] <- "Liver_Stroke"
dataPERIPHERAL$condition[dataPERIPHERAL$organ == "Spleen" & dataPERIPHERAL$group == "Naive"] <- "Spleen_Naive"
dataPERIPHERAL$condition[dataPERIPHERAL$organ == "Spleen" & dataPERIPHERAL$group == "Stroke"] <- "Spleen_Stroke"
dataPERIPHERAL$condition[dataPERIPHERAL$organ == "Blood" & dataPERIPHERAL$group == "Naive"] <- "Blood_Naive"
dataPERIPHERAL$condition[dataPERIPHERAL$organ == "Blood" & dataPERIPHERAL$group == "Stroke"] <- "Blood_Stroke"
head(dataPERIPHERAL@meta.data, 10000)

new.grouping <- c("sampleorgan")
dataPERIPHERAL[[new.grouping]] <- new.grouping
colnames(dataPERIPHERAL@meta.data)

dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Heart" & dataPERIPHERAL$sample.effect == "1"] <- "Heart_1"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Heart" & dataPERIPHERAL$sample.effect == "2"] <- "Heart_2"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Heart" & dataPERIPHERAL$sample.effect == "3"] <- "Heart_3"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Heart" & dataPERIPHERAL$sample.effect == "4"] <- "Heart_4"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Heart" & dataPERIPHERAL$sample.effect == "5"] <- "Heart_5"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Heart" & dataPERIPHERAL$sample.effect == "6"] <- "Heart_6"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Heart" & dataPERIPHERAL$sample.effect == "7"] <- "Heart_7"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Heart" & dataPERIPHERAL$sample.effect == "8"] <- "Heart_8"

dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Liver" & dataPERIPHERAL$sample.effect == "1"] <- "Liver_1"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Liver" & dataPERIPHERAL$sample.effect == "2"] <- "Liver_2"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Liver" & dataPERIPHERAL$sample.effect == "3"] <- "Liver_3"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Liver" & dataPERIPHERAL$sample.effect == "4"] <- "Liver_4"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Liver" & dataPERIPHERAL$sample.effect == "5"] <- "Liver_5"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Liver" & dataPERIPHERAL$sample.effect == "6"] <- "Liver_6"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Liver" & dataPERIPHERAL$sample.effect == "7"] <- "Liver_7"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Liver" & dataPERIPHERAL$sample.effect == "8"] <- "Liver_8"

dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Lung" & dataPERIPHERAL$sample.effect == "1"] <- "Lung_1"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Lung" & dataPERIPHERAL$sample.effect == "2"] <- "Lung_2"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Lung" & dataPERIPHERAL$sample.effect == "3"] <- "Lung_3"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Lung" & dataPERIPHERAL$sample.effect == "4"] <- "Lung_4"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Lung" & dataPERIPHERAL$sample.effect == "5"] <- "Lung_5"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Lung" & dataPERIPHERAL$sample.effect == "6"] <- "Lung_6"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Lung" & dataPERIPHERAL$sample.effect == "7"] <- "Lung_7"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Lung" & dataPERIPHERAL$sample.effect == "8"] <- "Lung_8"

dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Spleen" & dataPERIPHERAL$sample.effect == "1"] <- "Spleen_1"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Spleen" & dataPERIPHERAL$sample.effect == "2"] <- "Spleen_2"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Spleen" & dataPERIPHERAL$sample.effect == "3"] <- "Spleen_3"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Spleen" & dataPERIPHERAL$sample.effect == "4"] <- "Spleen_4"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Spleen" & dataPERIPHERAL$sample.effect == "5"] <- "Spleen_5"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Spleen" & dataPERIPHERAL$sample.effect == "6"] <- "Spleen_6"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Spleen" & dataPERIPHERAL$sample.effect == "7"] <- "Spleen_7"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Spleen" & dataPERIPHERAL$sample.effect == "8"] <- "Spleen_8"

dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Blood" & dataPERIPHERAL$sample.effect == "1"] <- "Blood_1"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Blood" & dataPERIPHERAL$sample.effect == "2"] <- "Blood_2"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Blood" & dataPERIPHERAL$sample.effect == "3"] <- "Blood_3"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Blood" & dataPERIPHERAL$sample.effect == "4"] <- "Blood_4"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Blood" & dataPERIPHERAL$sample.effect == "5"] <- "Blood_5"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Blood" & dataPERIPHERAL$sample.effect == "6"] <- "Blood_6"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Blood" & dataPERIPHERAL$sample.effect == "7"] <- "Blood_7"
dataPERIPHERAL$sampleorgan[dataPERIPHERAL$organ == "Blood" & dataPERIPHERAL$sample.effect == "8"] <- "Blood_8"

head(dataPERIPHERAL@meta.data, 10000)
tail(dataPERIPHERAL@meta.data, 10000)

GROUPaverages <- AverageExpression(dataPERIPHERAL, assay = "RNA", group.by = "sampleorgan") 
GROUPaverages <- as.data.frame(GROUPaverages)
head(GROUPaverages, n=100)

df2 <- data.frame(t(GROUPaverages))
df2$row_names <- row.names(df2)  
df3 <- df2[ , which(apply(df2, 2, var) != 0)]
ncol(df3)
ir.pca <- prcomp(df3[,1:20706], center = TRUE, scale. = TRUE)

ir.pca
summary(ir.pca)
loadings <- ir.pca$rotation
scores <- ir.pca$x
correlations <- t(loadings)*ir.pca$sdev
plot(ir.pca, type = "l",main="PCA")
ev <- ir.pca$sdev^2
evplot(ev)
df3$row_names <- row.names(df3)  

ggbiplot(ir.pca, choices=c(1,2), obs.scale = 1, var.scale = 1, varname.size = 0, var.axes = FALSE, groups = df3$row_names)

new.grouping <- c("group")
df3[[new.grouping]] <- new.grouping

df3$group[df3$row_names == "RNA.Heart_1"] <- "Heart_Naive"
df3$group[df3$row_names == "RNA.Heart_3"] <- "Heart_Naive"
df3$group[df3$row_names == "RNA.Heart_5"] <- "Heart_Naive"
df3$group[df3$row_names == "RNA.Heart_7"] <- "Heart_Naive"
df3$group[df3$row_names == "RNA.Heart_2"] <- "Heart_Stroke"
df3$group[df3$row_names == "RNA.Heart_4"] <- "Heart_Stroke"
df3$group[df3$row_names == "RNA.Heart_6"] <- "Heart_Stroke"
df3$group[df3$row_names == "RNA.Heart_8"] <- "Heart_Stroke"

df3$group[df3$row_names == "RNA.Liver_1"] <- "Liver_Naive"
df3$group[df3$row_names == "RNA.Liver_3"] <- "Liver_Naive"
df3$group[df3$row_names == "RNA.Liver_5"] <- "Liver_Naive"
df3$group[df3$row_names == "RNA.Liver_7"] <- "Liver_Naive"
df3$group[df3$row_names == "RNA.Liver_2"] <- "Liver_Stroke"
df3$group[df3$row_names == "RNA.Liver_4"] <- "Liver_Stroke"
df3$group[df3$row_names == "RNA.Liver_6"] <- "Liver_Stroke"
df3$group[df3$row_names == "RNA.Liver_8"] <- "Liver_Stroke"

df3$group[df3$row_names == "RNA.Lung_1"] <- "Lung_Naive"
df3$group[df3$row_names == "RNA.Lung_3"] <- "Lung_Naive"
df3$group[df3$row_names == "RNA.Lung_5"] <- "Lung_Naive"
df3$group[df3$row_names == "RNA.Lung_7"] <- "Lung_Naive"
df3$group[df3$row_names == "RNA.Lung_2"] <- "Lung_Stroke"
df3$group[df3$row_names == "RNA.Lung_4"] <- "Lung_Stroke"
df3$group[df3$row_names == "RNA.Lung_6"] <- "Lung_Stroke"
df3$group[df3$row_names == "RNA.Lung_8"] <- "Lung_Stroke"

df3$group[df3$row_names == "RNA.Spleen_1"] <- "Spleen_Naive"
df3$group[df3$row_names == "RNA.Spleen_3"] <- "Spleen_Naive"
df3$group[df3$row_names == "RNA.Spleen_5"] <- "Spleen_Naive"
df3$group[df3$row_names == "RNA.Spleen_7"] <- "Spleen_Naive"
df3$group[df3$row_names == "RNA.Spleen_2"] <- "Spleen_Stroke"
df3$group[df3$row_names == "RNA.Spleen_4"] <- "Spleen_Stroke"
df3$group[df3$row_names == "RNA.Spleen_6"] <- "Spleen_Stroke"
df3$group[df3$row_names == "RNA.Spleen_8"] <- "Spleen_Stroke"

df3$group[df3$row_names == "RNA.Blood_1"] <- "Blood_Naive"
df3$group[df3$row_names == "RNA.Blood_3"] <- "Blood_Naive"
df3$group[df3$row_names == "RNA.Blood_5"] <- "Blood_Naive"
df3$group[df3$row_names == "RNA.Blood_7"] <- "Blood_Naive"
df3$group[df3$row_names == "RNA.Blood_2"] <- "Blood_Stroke"
df3$group[df3$row_names == "RNA.Blood_4"] <- "Blood_Stroke"
df3$group[df3$row_names == "RNA.Blood_6"] <- "Blood_Stroke"
df3$group[df3$row_names == "RNA.Blood_8"] <- "Blood_Stroke"

ggbiplot(ir.pca, choices=c(1,2), obs.scale = 1, var.scale = 1, varname.size = 0, var.axes = FALSE, groups = df3$group, ellipse = TRUE, ellipse.prob = 0.5)

fviz_pca_ind(ir.pca, geom.ind = "point", pointshape = 21, 
             pointsize = 2, 
             col.ind = "black", 
             palette = c("red", "darkred", "blue", "darkblue", "orange", "darkorange", "green", "darkgreen", "pink", "purple"),
             #palette = "jco",
             fill.ind = df3$group,
             addEllipses = TRUE,
             ellipse.level=0.50,
             label = "var",
             col.var = "black",
             #repel = TRUE,
             legend.title = "DPCA ScSeq") +
        ggtitle("PCA") +
        theme(plot.title = element_text(hjust = 0.5))

ir.pca$x
pca <- ir.pca$x

Bnaive <- pca[c(1,3,5,7),]
Bstroke <- pca[c(2,4,6,8),]
Hnaive <- pca[c(9,11,13,15),]
Hstroke <- pca[c(10,12,14,16),]
Linaive <- pca[c(17,19,21,23),]
Listroke <- pca[c(18,20,22,24),]
Lunaive <- pca[c(25,27,29,31),]
Lustroke <- pca[c(26,28,30,32),]
Snaive <- pca[c(33,35,37,39),]
Sstroke <- pca[c(34,36,38,40),]

Heart_N <- apply(Hnaive, 2, mean)
Heart_S <- apply(Hstroke, 2, mean)
Liver_N <- apply(Linaive, 2, mean)
Liver_S <- apply(Listroke, 2, mean)
Lung_N <- apply(Lunaive, 2, mean)
Lung_S <- apply(Lustroke, 2, mean)
Spleen_N <- apply(Snaive, 2, mean)
Spleen_S <- apply(Sstroke, 2, mean)
Blood_N <- apply(Bnaive, 2, mean)
Blood_S <- apply(Bstroke, 2, mean)
res.pca <- rbind(Blood_N, Blood_S, Heart_N, Heart_S, Liver_N, Liver_S, Lung_N, Lung_S, Spleen_N, Spleen_S)

eu.dist <- dist(res.pca[,1:2],method = "euclidian")
ma.dist <- dist(res.pca[,1:2],method = "manhattan")

eu.dist <- as.matrix(eu.dist)
ma.dist <- as.matrix(ma.dist)





