library(Seurat)
library(ArchR)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
library(parallel)
library(dplyr)
library(Hmisc)
library(Cairo)

addArchRThreads(threads = 4)
addArchRGenome("mm10")

#####
#aux functions
#####
# 
featureToGR <- function(feature,pattern){
  featureSplit <- stringr::str_split(paste0(feature), pattern = pattern, n = 3, simplify = TRUE)
  gr <- GRanges(featureSplit[,1],IRanges(as.integer(featureSplit[,2]),as.integer(featureSplit[,3])))
  return(gr)
}

fastAnnoPeaks <- function(
    peaks = NULL,
    BSgenome = NULL,
    geneAnnotation = NULL,
    promoterRegion = c(2000, 100),
    logFile = NULL
){
  
  #Validate
  
  peakSummits <- resize(peaks,1,"center")
  BSgenome <- validBSgenome(BSgenome)
  
  #First Lets Get Distance to Nearest Gene Start
  distPeaks <- distanceToNearest(peakSummits, resize(geneAnnotation$genes, 1, "start"), ignore.strand = TRUE)
  mcols(peaks)$distToGeneStart <- mcols(distPeaks)$distance
  mcols(peaks)$nearestGene <- mcols(geneAnnotation$genes)$symbol[subjectHits(distPeaks)]
  promoters <- extendGR(resize(geneAnnotation$genes, 1, "start"), upstream = promoterRegion[1], downstream = promoterRegion[2])
  op <- overlapsAny(peakSummits, promoters, ignore.strand = TRUE)
  og <- overlapsAny(peakSummits, geneAnnotation$genes, ignore.strand = TRUE)
  oe <- overlapsAny(peakSummits, geneAnnotation$exons, ignore.strand = TRUE)
  type <- rep("Distal", length(peaks))
  type[which(og & oe)] <- "Exonic"
  type[which(og & !oe)] <- "Intronic"
  type[which(op)] <- "Promoter"
  mcols(peaks)$peakType <- type
  
  #First Lets Get Distance to Nearest TSS's
  distTSS <- distanceToNearest(peakSummits, resize(geneAnnotation$TSS, 1, "start"), ignore.strand = TRUE)
  mcols(peaks)$distToTSS <- mcols(distTSS)$distance
  if("symbol" %in% colnames(mcols(geneAnnotation$TSS))){
    mcols(peaks)$nearestTSS <- mcols(geneAnnotation$TSS)$symbol[subjectHits(distPeaks)]
  }else if("tx_name" %in% colnames(mcols(geneAnnotation$TSS))){
    mcols(peaks)$nearestTSS <- mcols(geneAnnotation$TSS)$tx_name[subjectHits(distPeaks)]
  }
  
  #Get NucleoTide Content
  nucFreq <- BSgenome::alphabetFrequency(getSeq(BSgenome, peaks))
  mcols(peaks)$GC <- round(rowSums(nucFreq[,c("G","C")]) / rowSums(nucFreq),4)
  mcols(peaks)$N <- round(nucFreq[,c("N")] / rowSums(nucFreq),4)
  peaks
  
}
#####
#Section 1 - Creation of ArchR object from scratch
#####

inputFiles <- getInputFiles("~/.../outs/")[1]
names(inputFiles) <- "Stroke"

ArrowFiles <- createArrowFiles(inputFiles,
                               minTSS = 0.1,
                               minFrags = 100,
                               maxFrags = 500000,
                               addGeneScoreMat = TRUE,
                               addTileMat = TRUE,
                               force = TRUE)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)
message('done doubScores')

proj <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = "stroke_full2",
  copyArrows = TRUE #This is recommened so that if you modify the Arrow files you have an original copy for later usage.
)

proj <- filterDoublets(ArchRProj = proj)
saveArchRProject(ArchRProj = proj, outputDirectory = ".",overwrite = T)
 
# #####
# # Subsetting cells from main ArchR project with cells from the coembed object
# #####
object <- readRDS("./coembed.rds")
subset_cells <- paste0('Stroke#', colnames(object@assays$peaks))
idxcode <- BiocGenerics::which(proj$cellNames %in% subset_cells)
cellscode <- proj$cellNames[idxcode]
subset_proj <- proj[cellscode, ]
saveArchRProject(ArchRProj = subset_proj, outputDirectory = ".",overwrite = TRUE)


# # #####
# # # Running commands for first pass placeholder
# # #####
subset_proj <- addIterativeLSI(
  ArchRProj = subset_proj,
  saveIterations = FALSE,
  useMatrix = "TileMatrix",
  depthCol = "nFrags",
  name = "LSI_ATAC"
)

subset_proj <- addClusters(subset_proj, reducedDims = "LSI_ATAC", name = "Clusters", resolution = 0.2, force = TRUE)
subset_proj <- addUMAP(subset_proj, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
# 
# #####
# # Adding the RNA object to redo integration, but with settings done previously to recapitulate
# # note that with the current settings, all the prediction scores table to >0.5, perhaps filter already present in internal setting
# #####
# 
rna_object <- readRDS("~/.../BM_scRNAseq.rds")
# 
subset_proj <- addImputeWeights(subset_proj, reducedDims = "LSI_ATAC")
subset_proj <- addGeneIntegrationMatrix(
  ArchRProj = subset_proj,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "LSI_ATAC",
  seRNA = rna_object,
  addToArrow = TRUE,
  force = TRUE,
  groupRNA = "group",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un",
  l2.norm=T,k.anchor = 20,k.filter=200,k.score = 30,max.features=500
)

#
# #####
# # Here reassigning the predicted clusters, conditions and UMAP calculated from the coembed object
# #####

cell_names <- paste0('Stroke#',colnames(object)) #get cell names for rna obj
#subset_proj <- subset_proj[cell_names]
archr_cellnames <- subset_proj$cellNames # get cell names for archr obj

seurat_umap <- Embeddings(object,'umap') # get embeddings from coembed object
row.names(seurat_umap) <- paste0('Stroke#',colnames(object)) #set rownames for embedding
seurat_umap <- seurat_umap[match(archr_cellnames,rownames(seurat_umap)),] #match rownames for umap
colnames(seurat_umap) <- colnames(subset_proj@embeddings$UMAP_ATAC$df) #set colnames for the umap df

seurat_clusters <- object$predicted.id # get predicted id from the coembed object
#names(seurat_clusters) <- paste0('Stroke#',cell_names)
names(seurat_clusters) <- cell_names
seurat_clusters <- seurat_clusters[match(archr_cellnames,cell_names)]

conditions <- data.frame(code=sapply(strsplit(names(object$orig.ident),"-"), `[`, 2))
conditions <- conditions %>% mutate(code_ID = ifelse(code == 1 | code == 5, 'naive',
                                                     ifelse(code == 4 | code == 7 | code == 8, 'stroke_treatment',
                                                            ifelse(code == 2 | code == 3 | code == 6, 'stroke', NA))))
seurat_cond <- conditions$code_ID
names(seurat_cond) <- paste0('Stroke#',cell_names)
seurat_cond <- seurat_cond[match(archr_cellnames,cell_names)]

subset_proj@embeddings$UMAP_ATAC$df <- seurat_umap
subset_proj$Clusters <- seurat_clusters
subset_proj$Cond <- seurat_cond
saveArchRProject(ArchRProj = subset_proj, outputDirectory = "./stroke_subset/",overwrite = T)

#
#####
# Adding the peak set manually and doing all the processing
#####

peak_set <-  featureToGR(row.names(object),'-')
chrs <- c('GL456233.1', 'GL456211.1', 'GL456212.1', 'JH584304.1', 'GL456216.1', 'JH584292.1 ')
#chrs <- c(paste("chr", 1:19, sep=''))
peak_set_final <- peak_set[grep('chr1|chr2|chr3', peak_set)]
peak_set_final <- peak_set[grep('GL456233.1|GL456211.1|GL456212.1|JH584304.1|GL456216.1|JH584292.1', peak_set, invert = TRUE)]
BSgenome <- ArchR::validBSgenome('mm10')

keepBSgenomeSequences <- function(genome, seqnames)
{
  stopifnot(all(seqnames %in% seqnames(genome)))
  genome@user_seqnames <- setNames(seqnames, seqnames)
  genome@seqinfo <- genome@seqinfo[seqnames]
  genome
}

sequences_to_keep <- paste0("chr", c(1:19, 'X', 'Y'))
BSgenome <- keepBSgenomeSequences(BSgenome, sequences_to_keep)

geneAnnotation <- subset_proj@geneAnnotation
#

peak_set_final <- fastAnnoPeaks(
  peak_set_final,
  BSgenome = BSgenome,
  geneAnnotation = geneAnnotation,
  promoterRegion = c(2000, 100)
)
subset_proj<- addPeakSet(ArchRProj = subset_proj, peakSet = peak_set_final, force = T)

subset_proj$ClusterCond <- paste0(subset_proj$Clusters,'_',subset_proj$Cond)
subset_proj <- addGroupCoverages(ArchRProj = subset_proj, groupBy = "Clusters",force=T)
subset_proj <- addGroupCoverages(ArchRProj = subset_proj, groupBy = "Cond",force=T)
subset_proj <- addGroupCoverages(ArchRProj = subset_proj, groupBy = "ClusterCond",force=T)
subset_proj <- addPeakMatrix(subset_proj)
subset_proj <- addMotifAnnotations(ArchRProj = subset_proj,motifSet = "JASPAR2020", annoName = "Motif",force=T)

subset_proj <- addBgdPeaks(subset_proj, method = 'ArchR')
message('done bgd peaks')
subset_proj <- addDeviationsMatrix(
  ArchRProj = subset_proj,
  peakAnnotation = "Motif",
  force = TRUE
)
# 
# 
# #################
# # Marker Peaks (following vignette)
# ###################
markersPeaks <- getMarkerFeatures(
  ArchRProj = subset_proj,
  useMatrix = "PeakMatrix",
  groupBy = "ClusterCond",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useGroups = "0_naive",
  bgdGroups = "0_stroke"
  
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1", returnGR = TRUE)

plotMarkerHeatmap(
  seMarker = markersPeaks,
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
# 
plotMarkers(seMarker = markersPeaks, name = "HSPC", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")

# 
# #################
# # Peak to gene links (following vignette)
# ###################
subset_proj <- addCoAccessibility(
  ArchRProj = subset_proj,
  reducedDims = "LSI_ATAC"
)

subset_proj <- addPeak2GeneLinks(
  ArchRProj = subset_proj,
  reducedDims = "LSI_ATAC"
)

p2g <- getPeak2GeneLinks(
  ArchRProj = subset_proj,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = TRUE
)
#

set.seed(123)
p <- plotPeak2GeneHeatmap(ArchRProj = subset_proj, groupBy = c("ClusterCond"), returnMatrices = TRUE, k = 25)
peak2genetable <- p@listData$Peak2GeneLinks
