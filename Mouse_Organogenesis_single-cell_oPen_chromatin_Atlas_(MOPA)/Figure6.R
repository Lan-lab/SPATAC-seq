setwd("/home/sunkeyong/MOPA_project/myocytes/")
#
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(TFBSTools)
library(ArchR)
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(pheatmap)
library(ComplexHeatmap)
library(factoextra)
library(matrixStats)
library("RColorBrewer")
library("circlize")
#
#
ArrowFiles <- c("/home/sunkeyong/MOPA_project/Major_cell_types_arrowfile/Myocytes.arrow")
#
MOPA_Myocytes <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  copyArrows = TRUE)
#
MOPA_Myocytes
load("/home/sunkeyong/MOPA_project/myocytes/myocytes.peakset.RData")
#
MOPA_Myocytes <- addPeakSet(MOPA_Myocytes,peakSet = mmyo.peakset)
MOPA_Myocytes <- addPeakMatrix(MOPA_Myocytes)
getAvailableMatrices(MOPA_Myocytes)
##
MOPA_Myocytes <- addIterativeLSI(
  ArchRProj = MOPA_Myocytes,
  useMatrix = "PeakMatrix", 
  name = "Peak_IterativeLSI", 
  iterations = 3, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(2,3), 
    sampleCells = 10000,
    maxClusters= 300, 
    n.start = 10
  ), 
  varFeatures = 100000, 
  dimsToUse = 1:50
)
#
MOPA_Myocytes <- addClusters(input = MOPA_Myocytes,reducedDims = "Peak_IterativeLSI",method = "Seurat",
                      name = "Peak_Clusters1", dimsToUse = 1:50,resolution = 0.5,maxClusters= 300,
                      force = T,sampleCells = NULL)
MOPA_Myocytes <- addClusters(input = MOPA_Myocytes,reducedDims = "Peak_IterativeLSI",method = "Seurat",
                      name = "Peak_Clusters2", dimsToUse = 1:50,resolution = 1,maxClusters= 300,
                      force = T,sampleCells = NULL)
MOPA_Myocytes <- addClusters(input = MOPA_Myocytes,reducedDims = "Peak_IterativeLSI",method = "Seurat",
                      name = "Peak_Clusters3", dimsToUse = 1:50,resolution = 2,maxClusters= 300,
                      force = T,sampleCells = NULL)
#
MOPA_Myocytes <- addUMAP(ArchRProj = MOPA_Myocytes, reducedDims = "Peak_IterativeLSI", name = "Peak_UMAP", 
                  nNeighbors = 50, dimsToUse = 1:50, minDist = 0.8, metric = "cosine",force = T)
##
MOPA_Myocytes <- addTSNE(ArchRProj = MOPA_Myocytes, reducedDims = "Peak_IterativeLSI", name = "Peak_TSNE", 
                  perplexity = 100, dimsToUse = 1:50,maxIterations = 1000,force = T)
MOPA_Myocytes <- addTSNE(ArchRProj = MOPA_Myocytes, reducedDims = "Peak_IterativeLSI", name = "Peak_TSNE2", 
                  perplexity = 200, dimsToUse = 1:30,maxIterations = 1500,force = T)
MOPA_Myocytes <- addTSNE(ArchRProj = MOPA_Myocytes, reducedDims = "Peak_IterativeLSI", name = "Peak_TSNE3", 
                  perplexity = 200, dimsToUse = 1:50,maxIterations = 1500,force = T)
########
MOPA_Myocytes <- addUMAP(ArchRProj = MOPA_Myocytes, reducedDims = "Peak_IterativeLSI", name = "Peak_UMAP2", 
                  nNeighbors = 100, dimsToUse = 1:50, minDist = 0.8, metric = "cosine",force = T)
MOPA_Myocytes <- addUMAP(ArchRProj = MOPA_Myocytes, reducedDims = "Peak_IterativeLSI", name = "Peak_UMAP3", 
                  nNeighbors = 200, dimsToUse = 1:50, minDist = 0.8, metric = "cosine",force = T)
MOPA_Myocytes <- addUMAP(ArchRProj = MOPA_Myocytes, reducedDims = "Peak_IterativeLSI", name = "Peak_UMAP4", 
                  nNeighbors = 100, dimsToUse = 1:50, minDist = 0.5, metric = "cosine",force = T)
MOPA_Myocytes <- addUMAP(ArchRProj = MOPA_Myocytes, reducedDims = "Peak_IterativeLSI", name = "Peak_UMAP5", 
                  nNeighbors = 100, dimsToUse = 1:50, minDist = 0.3, metric = "cosine",force = T)
#
MOPA_Myocytes <- addUMAP(ArchRProj = MOPA_Myocytes, reducedDims = "Peak_IterativeLSI", name = "Peak_UMAP6", 
                  nNeighbors = 100, dimsToUse = 2:50, minDist = 0.3, metric = "cosine",force = T)
#
###
##
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "Peak_Clusters1", embedding = "Peak_UMAP")
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "Peak_Clusters1", embedding = "Peak_UMAP2")
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "Peak_Clusters1", embedding = "Peak_UMAP3")
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "Peak_Clusters1", embedding = "Peak_UMAP4")
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "Peak_Clusters1", embedding = "Peak_UMAP5")
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "Peak_Clusters3", embedding = "Peak_UMAP5")
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "Peak_Clusters3", embedding = "Peak_TSNE")
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "Peak_Clusters3", embedding = "Peak_TSNE2")
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "Peak_Clusters3", embedding = "Peak_TSNE3")
#
p1 <- plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "Peak_Clusters1", embedding = "Peak_UMAP")
p2 <- plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "Peak_Clusters3", embedding = "Peak_TSNE")
ggAlignPlots(p1, p2, type = "h")
#
p1 <- plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "Sample", embedding = "TSNE")
ggAlignPlots(p1, p2, type = "h")
#
MOPA_Myocytes$embryo_day <- myo.metadata$embryo_day
p1 <- plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "embryo_day", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "embryo_day", embedding = "TSNE")
ggAlignPlots(p1, p2, type = "h")
####
##
getAvailableMatrices(MOPA_Myocytes)
MOPA_Myocytes <- addImputeWeights(MOPA_Myocytes,reducedDims = "Peak_IterativeLSI", dimsToUse = 1:50)
######
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Myf5"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Pax3"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Pax5"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Pax7"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Myh3"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Myod1"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Myog"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Lhx2"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
##
####
MOPA_Myocytes <- addMotifAnnotations(ArchRProj = MOPA_Myocytes, motifSet = "cisbp", name = "Motif")
MOPA_Myocytes <- addBgdPeaks(MOPA_Myocytes,force = TRUE)
MOPA_Myocytes <- addDeviationsMatrix(
  ArchRProj = MOPA_Myocytes, 
  peakAnnotation = "Motif",
  force = TRUE
)
#
#
#
load("TOME_myocyte.scRNA.RData")
DimPlot(TOME_myocyte)
FeaturePlot(TOME_myocyte,features = c("Pax3"))
FeaturePlot(TOME_myocyte,features = c("Pax5"))
FeaturePlot(TOME_myocyte,features = c("Pax7"))
FeaturePlot(TOME_myocyte,features = c("Myh3"))
FeaturePlot(TOME_myocyte,features = c("Myf5"))
FeaturePlot(TOME_myocyte,features = c("Myod1"))
FeaturePlot(TOME_myocyte,features = c("Myog"))
FeaturePlot(TOME_myocyte,features = c("Lhx2"))
#
MOPA_Myocytes <- addGeneIntegrationMatrix(
  ArchRProj = MOPA_Myocytes, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "Peak_IterativeLSI",
  seRNA = TOME_myocyte,
  addToArrow = T,
  force = TRUE,
  groupRNA = "sub_cluster_id",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un"
)
####
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "predictedGroup_Un", embedding = "Peak_UMAP5")
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "predictedScore_Un", embedding = "Peak_UMAP5")
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "DoubletScore", embedding = "Peak_UMAP5")
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "Timepoint", embedding = "Peak_UMAP5")
####
getAvailableMatrices(MOPA_Myocytes)
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Lhx2"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Lhx2"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Pax3"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Pax3"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Pax7"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Pax7"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Myf5"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Myf5"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Myod1"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Myod1"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Myog"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Myh3"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Myog"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Lhx2"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
#
motifs <- c("Mef")
markerMotifs <- getFeatures(MOPA_Myocytes, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
markerMotifs
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "MotifMatrix", name = c("z:Myf5_778"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "MotifMatrix", name = c("z:Myog_45"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "MotifMatrix", name = c("z:Myf6_63"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "MotifMatrix", name = c("z:Myod1_24"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "MotifMatrix", name = c("z:Pax3_613"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "MotifMatrix", name = c("z:Pax7_467"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "MotifMatrix", name = c("z:Trp63_853"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "MotifMatrix", name = c("z:Mef2d_842"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "MotifMatrix", name = c("z:Mef2c_638"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "MotifMatrix", name = c("z:Trp63_853"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
###########
getAvailableMatrices(MOPA_Myocytes)
markersGS <- getMarkerFeatures(
  ArchRProj = MOPA_Myocytes, 
  useMatrix = "GeneIntegrationMatrix", 
  groupBy = "Peak_Clusters1",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  #labelMarkers = markerGenes,
  transpose = TRUE
)
#
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
#
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$C11
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Isl1"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Isl1"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Pcolce2"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Pcolce2"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Mark1"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Mark1"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Htr1a"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Htr1a"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Ntn1"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Ntn1"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Robo1"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Robo1"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Robo2"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Robo2"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
##
#
markerList$C3$name[1:12]
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Pax7"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Pax7"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Trp63"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Trp63"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Meox1"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Meox1"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Mef2c"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Mef2c"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Mef2d"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Mef2d"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Meox1"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Meox1"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Lhx9"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Lhx9"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
#
FeaturePlot(TOME_myocyte,features = c("Trp63"))
FeaturePlot(TOME_myocyte,features = c("Meox1"))
FeaturePlot(TOME_myocyte,features = c("Gli2"))
FeaturePlot(TOME_myocyte,features = c("Hoxa9"))
FeaturePlot(TOME_myocyte,features = c("Hoxa3"))
FeaturePlot(TOME_myocyte,features = c("Lhx2"))
FeaturePlot(TOME_myocyte,features = c("Pax3"))
FeaturePlot(TOME_myocyte,features = c("Pax5"))
FeaturePlot(TOME_myocyte,features = c("Pax7"))
FeaturePlot(TOME_myocyte,features = c("Myf5"))
FeaturePlot(TOME_myocyte,features = c("Gdnf"))
FeaturePlot(TOME_myocyte,features = c("Sim2"))
FeaturePlot(TOME_myocyte,features = c("Gabra2"))
FeaturePlot(TOME_myocyte,features = c("Hoxb9"))
FeaturePlot(TOME_myocyte,features = c("Taf4b"))
FeaturePlot(TOME_myocyte,features = c("Fmr1nb"))
FeaturePlot(TOME_myocyte,features = c("Meox2"))
FeaturePlot(TOME_myocyte,features = c("Cnr2"))
FeaturePlot(TOME_myocyte,features = c("Ddx39"))
FeaturePlot(TOME_myocyte,features = c("Mef2c"))
FeaturePlot(TOME_myocyte,features = c("Mef2d"))
FeaturePlot(TOME_myocyte,features = c("Trp63"))
FeaturePlot(TOME_myocyte,features = c("Myod1"))
FeaturePlot(TOME_myocyte,features = c("Sp5"))
FeaturePlot(TOME_myocyte,features = c("Pbx3"))
FeaturePlot(TOME_myocyte,features = c("Etv1"))
FeaturePlot(TOME_myocyte,features = c("Dmrt2"))
FeaturePlot(TOME_myocyte,features = c("Msc"))
FeaturePlot(TOME_myocyte,features = c("Robo2"))
FeaturePlot(TOME_myocyte,features = c("Robo1"))
FeaturePlot(TOME_myocyte,features = c("Msc"))
FeaturePlot(TOME_myocyte,features = c("Msc"))
FeaturePlot(TOME_myocyte,features = c("Dnm3"))
FeaturePlot(TOME_myocyte,features = c("Meis2"))
FeaturePlot(TOME_myocyte,features = c("Apobec2"))
#
p <- plotBrowserTrack(ArchRProj = MOPA_Myocytes, groupBy = "Peak_Clusters1", geneSymbol = "Hoxa9", 
                      upstream = 50000,downstream = 50000,ylim = c(0.1,1))
grid::grid.newpage()
grid::grid.draw(p$Hoxa9)
p <- plotBrowserTrack(ArchRProj = MOPA_Myocytes, groupBy = "Peak_Clusters1", geneSymbol = "Myog", 
                      upstream = 70000,downstream = 70000,ylim = c(0.1,1))
grid::grid.newpage()
grid::grid.draw(p$Myog)
p <- plotBrowserTrack(ArchRProj = MOPA_Myocytes, groupBy = "Peak_Clusters1", geneSymbol = "Pax3", 
                      upstream = 70000,downstream = 70000,ylim = c(0.1,1))
grid::grid.newpage()
grid::grid.draw(p$Pax3)
p <- plotBrowserTrack(ArchRProj = MOPA_Myocytes, groupBy = "Peak_Clusters1", geneSymbol = "Trp63", 
                      upstream = 70000,downstream = 300000,ylim = c(0.1,1))
grid::grid.newpage()
grid::grid.draw(p$Trp63)
p <- plotBrowserTrack(ArchRProj = MOPA_Myocytes, groupBy = "Peak_Clusters1", geneSymbol = "Meox1", 
                      upstream = 30000,downstream = 30000,ylim = c(0.01,1))
grid::grid.newpage()
grid::grid.draw(p$Meox1)
p <- plotBrowserTrack(ArchRProj = MOPA_Myocytes, groupBy = "Peak_Clusters1", geneSymbol = "Mef2c", 
                      upstream = 50000,downstream = 220000,ylim = c(0.01,1))
grid::grid.newpage()
grid::grid.draw(p$Mef2c)
p <- plotBrowserTrack(ArchRProj = MOPA_Myocytes, groupBy = "Peak_Clusters1", geneSymbol = "Mef2d", 
                      upstream = 20000,downstream = 50000,ylim = c(0.01,1))
grid::grid.newpage()
grid::grid.draw(p$Mef2d)
p <- plotBrowserTrack(ArchRProj = MOPA_Myocytes, groupBy = "Peak_Clusters1", geneSymbol = "Bcl11b", 
                      upstream = 20000,downstream = 50000,ylim = c(0.01,1))
grid::grid.newpage()
grid::grid.draw(p$Bcl11b)
p <- plotBrowserTrack(ArchRProj = MOPA_Myocytes, groupBy = "Peak_Clusters1", geneSymbol = "Apobec2", 
                      upstream = 30000,downstream =10000,ylim = c(0.01,1))
grid::grid.newpage()
grid::grid.draw(p$Apobec2)
p <- plotBrowserTrack(ArchRProj = MOPA_Myocytes, groupBy = "Peak_Clusters1", geneSymbol = "Meis3", 
                      upstream = 20000,downstream =20000,ylim = c(0.01,1),loops = getPeak2GeneLinks(MOPA_Myocytes))
grid::grid.newpage()
grid::grid.draw(p$Meis3)
p <- plotBrowserTrack(ArchRProj = MOPA_Myocytes, groupBy = "Peak_Clusters1", geneSymbol = "Twist1", 
                      upstream = 20000,downstream =20000,ylim = c(0.01,1),loops = getPeak2GeneLinks(MOPA_Myocytes))
grid::grid.newpage()
grid::grid.draw(p$Twist1)
##
MOPA_Myocytes <- addCoAccessibility(
  ArchRProj = MOPA_Myocytes,
  dimsToUse = 1:30,
  scaleDims = NULL,
  corCutOff = 0.75,
  reducedDims = "Peak_IterativeLSI"
)
#
#####
p <- plotBrowserTrack(ArchRProj = MOPA_Myocytes, groupBy = "Peak_Clusters1", geneSymbol = "Myog", 
                      upstream = 80000,downstream = 80000,ylim = c(0.01,1),
                      loops = getCoAccessibility(MOPA_Myocytes))
grid::grid.newpage()
grid::grid.draw(p$Myog)
####
MOPA_Myocytes <- addPeak2GeneLinks(
  ArchRProj = MOPA_Myocytes,
  dimsToUse = 1:50,
  corCutOff = 0.6,
  reducedDims = "Peak_IterativeLSI"
)
#
p2g <- getPeak2GeneLinks(
  ArchRProj = MOPA_Myocytes,
  corCutOff = 0.4,
  resolution = 1,
  returnLoops = FALSE
)
#####
# myog
p <- plotBrowserTrack(ArchRProj = MOPA_Myocytes, groupBy = "Peak_Clusters1", geneSymbol = "Myog",
                      upstream = 80000,downstream = 80000,ylim = c(0.1,1),
                      loops = getPeak2GeneLinks(MOPA_Myocytes))
grid::grid.newpage()
grid::grid.draw(p$Myog)
###
###
###
p <- plotPeak2GeneHeatmap(ArchRProj = MOPA_Myocytes, groupBy = "Peak_Clusters1")
p
###
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "Peak_Clusters1", embedding = "Peak_UMAP5")
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "Peak_Clusters1", embedding = "Peak_UMAP5")
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "Peak_Clusters3", embedding = "Peak_UMAP5")
trajectory <- c("C15","C22","C23","C21","C14","C13","C9","C10","C12","C11","C7","C8","C6","C3","C1","C2")
trajectory
MOPA_Myocytes <- addTrajectory(
  ArchRProj = MOPA_Myocytes,
  reducedDims = "Peak_IterativeLSI",
  name = "MBMT", 
  groupBy = "Peak_Clusters3",
  trajectory = trajectory, 
  embedding = "Peak_UMAP5", 
  force = TRUE
)
p <- plotTrajectory(MOPA_Myocytes, embedding = "Peak_UMAP5",trajectory = "MBMT", colorBy = "cellColData", name = "MBMT")
p[[1]]
#
#
###
trajMM  <- getTrajectory(ArchRProj = MOPA_Myocytes, name = "MBMT", useMatrix = "MotifMatrix", log2Norm = FALSE)
trajGSM <- getTrajectory(ArchRProj = MOPA_Myocytes, name = "MBMT", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
trajGIM <- getTrajectory(ArchRProj = MOPA_Myocytes, name = "MBMT", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
trajPM  <- getTrajectory(ArchRProj = MOPA_Myocytes, name = "MBMT", useMatrix = "PeakMatrix", log2Norm = TRUE)
####
p4 <- plotTrajectoryHeatmap(trajGIM, pal = paletteContinuous(set = "solarExtra"))
p4
##
###
corGSM_MM <- correlateTrajectories(trajGSM, trajMM)
corGSM_MM[[1]]
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]
trajCombined <- trajGSM2
assay(trajCombined) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "horizonExtra"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
ht1 + ht2
##
corGIM_MM <- correlateTrajectories(trajGIM, trajMM,corCutOff = 0.1)
corGIM_MM[[1]]
trajGIM2 <- trajGIM[corGIM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGIM_MM[[1]]$name2, ]
trajCombined <- trajGIM2
assay(trajCombined) <- t(apply(assay(trajGIM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGIM2))
ht1 <- plotTrajectoryHeatmap(trajGIM2,  pal = paletteContinuous(set = "blueYellow"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "solarExtra"), varCutOff = 0, rowOrder = rowOrder)
ht1 + ht2
###
load("../Figure1/MOPA.seurat_type")
#
MOPA.seurat_type$cells <- unlist(strsplit(rownames(MOPA.seurat_type@meta.data),"#"))[seq(2,2*357171,2)]
MOPA.seurat_type.metadata <- as.data.frame(MOPA.seurat_type@meta.data)
rownames(MOPA.seurat_type.metadata) <- MOPA.seurat_type.metadata$cells
MOPA.seurat_type.metadata.tmp <- MOPA.seurat_type.metadata[unlist(strsplit(MOPA_Myocytes$cellNames,"#"))[seq(2,2*10837,2)],]
#
MOPA_Myocytes$DoubletScore <-  as.numeric(MOPA.seurat_type.metadata.tmp$DoubletScore)
MOPA_Myocytes$DoubletEnrich <-  as.numeric(MOPA.seurat_type.metadata.tmp$DoubletEnrich)
MOPA_Myocytes$Timepoint <- as.character(MOPA.seurat_type.metadata.tmp$Timepoint)
#
#
#
x.num <- c(rep(1,800000),rep(0,800000),rep(2,800000))
head(x.num)
x.num.sample <- cbind(sample(x.num,length(MOPA_Myocytes$cellNames)),
                      sample(x.num,length(MOPA_Myocytes$cellNames)),
                      sample(x.num,length(MOPA_Myocytes$cellNames)),
                      sample(x.num,length(MOPA_Myocytes$cellNames)),
                      sample(x.num,length(MOPA_Myocytes$cellNames)))
head(x.num.sample)
colnames(x.num.sample) <- c("F1","F2","F3","F4","F5")
rownames(x.num.sample) <- MOPA_Myocytes$cellNames
#
MOPA_Myocytes_seurat <- CreateSeuratObject(counts = t(x.num.sample), assay = "ATAC", project = "CardiacMuscle", min.cells = 0, min.features = 0)
MOPA_Myocytes_seurat
MOPA_Myocytes_seurat <- NormalizeData(MOPA_Myocytes_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
MOPA_Myocytes_seurat <- FindVariableFeatures(MOPA_Myocytes_seurat, selection.method = "vst", nfeatures = 4)
MOPA_Myocytes_seurat <- ScaleData(MOPA_Myocytes_seurat, features = rownames(MOPA_Myocytes_seurat))
variablegene <- VariableFeatures(object = MOPA_Myocytes_seurat)
MOPA_Myocytes_seurat <- RunPCA(MOPA_Myocytes_seurat, features = variablegene,npcs =2)
DimPlot(MOPA_Myocytes_seurat)
#
x <- MOPA_Myocytes@embeddings$Peak_UMAP5
y <- x$df
s.UMAP <- cbind(y$`Peak_IterativeLSI#UMAP_Dimension_1`,y$`Peak_IterativeLSI#UMAP_Dimension_2`)
#
dim(s.UMAP)
head(s.UMAP)
rownames(s.UMAP) <- rownames(y)
head(s.UMAP)
colnames(s.UMAP) <- c("PC_1","PC_2")
head(s.UMAP)
class(s.UMAP)
#
MOPA_Myocytes_seurat@reductions$pca@cell.embeddings <- s.UMAP
table(MOPA_Myocytes_seurat@active.ident)
DimPlot(MOPA_Myocytes_seurat)
DimPlot(MOPA_Myocytes_seurat,split.by = "batch")
DimPlot(MOPA_Myocytes_seurat,group.by = "sample")
#
MOPA_Myocytes_seurat$DoubletEnrich <- MOPA_Myocytes$DoubletEnrich
MOPA_Myocytes_seurat$DoubletScore <- MOPA_Myocytes$DoubletScore
MOPA_Myocytes_seurat$Timepoint <- MOPA_Myocytes$Timepoint
MOPA_Myocytes_seurat$predictedScore_Un <- MOPA_Myocytes$predictedScore_Un
MOPA_Myocytes_seurat$Peak_Clusters3 <- MOPA_Myocytes$Peak_Clusters3
MOPA_Myocytes_seurat$Peak_Clusters1 <- MOPA_Myocytes$Peak_Clusters1
#
FeaturePlot(MOPA_Myocytes_seurat,features = c("DoubletScore"))
FeaturePlot(MOPA_Myocytes_seurat,features = c("DoubletEnrich"))
FeaturePlot(MOPA_Myocytes_seurat,features = c("predictedScore_Un"))
#
DimPlot(MOPA_Myocytes_seurat,group.by = "Timepoint")
DimPlot(MOPA_Myocytes_seurat,group.by = "Peak_Clusters3",label = T,repel = T)
DimPlot(MOPA_Myocytes_seurat,group.by = "Peak_Clusters1",label = T,repel = T)
DimPlot(MOPA_Myocytes_seurat,group.by = "Timepoint")

#
#
#
MOPA_Myocytes <- addHarmony(
  ArchRProj = MOPA_Myocytes,
  reducedDims = "Peak_IterativeLSI",
  name = "PeakMatrix_Harmony250",
  groupBy = "Timepoint",
  dimsToUse = 2:50,
  force = T
)
MOPA_Myocytes <- addClusters(input = MOPA_Myocytes,reducedDims = "PeakMatrix_Harmony250",method = "Seurat",
                      name = "Peak_H_C1", dimsToUse = 1:50,resolution = 2,maxClusters= 300,
                      force = T,sampleCells = NULL)
#
MOPA_Myocytes <- addUMAP(ArchRProj = MOPA_Myocytes, reducedDims = "PeakMatrix_Harmony250", name = "Peak_UMAP_H1", 
                  nNeighbors = 50, dimsToUse = 1:50, minDist = 0.8, metric = "cosine",force = T)
##
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "Timepoint", embedding = "Peak_UMAP_H1")
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Lhx2"), embedding = "Peak_UMAP_H1",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Lhx2"), embedding = "Peak_UMAP_H1",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Pax3"), embedding = "Peak_UMAP_H1",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Pax3"), embedding = "Peak_UMAP_H1",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Myog"), embedding = "Peak_UMAP_H1",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Myog"), embedding = "Peak_UMAP_H1",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Myod1"), embedding = "Peak_UMAP_H1",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Myod1"), embedding = "Peak_UMAP_H1",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Lhx9"), embedding = "Peak_UMAP_H1",imputeWeights = getImputeWeights(MOPA_Myocytes))
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneIntegrationMatrix", name = c("Lhx9"), embedding = "Peak_UMAP_H1",imputeWeights = getImputeWeights(MOPA_Myocytes))
#
MOPA_Myocytes <- addUMAP(ArchRProj = MOPA_Myocytes, reducedDims = "Peak_IterativeLSI", name = "Peak_UMAP6", 
                  nNeighbors = 100, dimsToUse = 2:50, minDist = 0.3, metric = "cosine",force = T)
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "Timepoint", embedding = "Peak_UMAP6")

#
#
#
M6genes <- c("Arfgef1","Ncoa2","Tmem131","Map4k4","Pard3b","Xrcc5","Agap1","Enah","Tgfb2",
             "Plxna2","Acvr2a","Cdh4","Zfp704","Dpyd","Pitx2","Ppp3ca","Bmpr1b","Hs2st1","Chd7","Tgfbr1",
             "Pde4b","Macf1","Camta1","Cacna2d1","Cnot6l","Fras1","Hip1",
             "Col26a1","Actb","Pdzrn3","Tmcc1","Mef2a","St8sia2","Akap13","Hbb-y","H19","Kcnq1ot1",
             "Palld","Large1","Gse1","Pard3","Jam3","Pknox2","Ncam1","Herc1","Mir6236","Rnf217","Marcks","Micu1","Grip1",
             "Ehbp1","Rtn4","Sptbn1","Msi2","Tanc2","Rptor","Rcor1","2010111I01Rik","Zswim6","Parp8","Enox1","Ghr","Zfr",
             "Fam49b","Rbfox2","Nectin3","Tulp4","Qk","Airn","Hsp90ab1",
             "Arhgap28","Svil","Zeb1","Nedd4l","Tmem2","Dock11","Eda","Diaph2","Mid1","mt-Nd5")
#
ArchRgenes <- getGenes(MOPA_Myocytes)
ArchRgenes <- as.data.frame(ArchRgenes$symbol)
head(ArchRgenes)
colnames(ArchRgenes) <- "genename"
#
M6genes2 <- intersect(ArchRgenes$genename,M6genes)
M6genes3 <- M6genes2[sample(78,20)]
#
features <- list(M6genes2Score = M6genes2,
                 M6genes3Score = M6genes3)
#
MOPA_Myocytes <- addModuleScore(ArchRProj = MOPA_Myocytes,useMatrix = "GeneScoreMatrix", features =features )
#
plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.M6genes2Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()
plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.M6genes3Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()
#
modulescore <- read.csv("scimodule.csv")
head(modulescore)
#
mM1 <- intersect(ArchRgenes$genename, subset(modulescore,gene_module=="1")$gene_short_name)
mM2 <- intersect(ArchRgenes$genename, subset(modulescore,gene_module=="2")$gene_short_name)
mM3 <- intersect(ArchRgenes$genename, subset(modulescore,gene_module=="3")$gene_short_name)
mM4 <- intersect(ArchRgenes$genename, subset(modulescore,gene_module=="4")$gene_short_name)
mM5 <- intersect(ArchRgenes$genename, subset(modulescore,gene_module=="5")$gene_short_name)
mM6 <- intersect(ArchRgenes$genename, subset(modulescore,gene_module=="6")$gene_short_name)
mM7 <- intersect(ArchRgenes$genename, subset(modulescore,gene_module=="7")$gene_short_name)
mM8 <- intersect(ArchRgenes$genename, subset(modulescore,gene_module=="8")$gene_short_name)
mM9 <- intersect(ArchRgenes$genename, subset(modulescore,gene_module=="9")$gene_short_name)[1:500]
mM10 <- intersect(ArchRgenes$genename, subset(modulescore,gene_module=="10")$gene_short_name)
mM11 <- intersect(ArchRgenes$genename, subset(modulescore,gene_module=="11")$gene_short_name)
mM12 <- intersect(ArchRgenes$genename, subset(modulescore,gene_module=="12")$gene_short_name)
mM13 <- intersect(ArchRgenes$genename, subset(modulescore,gene_module=="13")$gene_short_name)
mM14 <- intersect(ArchRgenes$genename, subset(modulescore,gene_module=="14")$gene_short_name)
#
features <- list(mM1Score = mM1,
                 mM2Score = mM2,
                 mM3Score = mM3,
                 mM4Score = mM4,
                 mM5Score = mM5,
                 mM6Score = mM6,
                 mM7Score = mM7,
                 mM8Score = mM8,
                 mM9Score = mM9,
                 mM10Score = mM10,                 
                 mM11Score = mM11,
                 mM12Score = mM12,
                 mM13Score = mM13,
                 mM14Score = mM14)
#
MOPA_Myocytes <- addModuleScore(ArchRProj = MOPA_Myocytes,useMatrix = "GeneScoreMatrix", features =features )
#
plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM14Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()
#
mM9.1 <- intersect(ArchRgenes$genename, subset(modulescore,gene_module=="9")$gene_short_name)[1:100]
mM9.2 <- intersect(ArchRgenes$genename, subset(modulescore,gene_module=="9")$gene_short_name)[101:200]
mM9.3 <- intersect(ArchRgenes$genename, subset(modulescore,gene_module=="9")$gene_short_name)[201:300]
mM9.4 <- intersect(ArchRgenes$genename, subset(modulescore,gene_module=="9")$gene_short_name)[301:400]
mM9.5 <- intersect(ArchRgenes$genename, subset(modulescore,gene_module=="9")$gene_short_name)[401:500]
#
features <- list(mM9.1Score = mM9.1,
                 mM9.2Score = mM9.2)
features <- list(mM9.3Score = mM9.3,
                 mM9.4Score = mM9.4,
                 mM9.5Score = mM9.5)
features <- list(mM8Score = mM8,
                 mM9Score = mM9)
#
plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM1Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()
plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM2Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()
plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM3Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()
plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM4Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()
plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM5Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()
plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM6Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()
plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM7Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()
plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM8Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()
plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM9Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()
plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM10Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()
plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM11Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()
plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM12Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()
plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM13Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()
plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM14Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()
#
plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM9.1Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "Timepoint", embedding = "Peak_UMAP5")
#
MOPA_Myocytes_seurat@active.ident <- as.factor(MOPA_Myocytes_seurat$Peak_Clusters3)
MOPA_Myocytes_seurat <- RenameIdents(MOPA_Myocytes_seurat,"C2"="C2357","C3"="C2357","C5"="C2357","C7"="C2357")
MOPA_Myocytes_seurat <- RenameIdents(MOPA_Myocytes_seurat,"C8"="C89","C9"="C89")
#MOPA_Myocytes_seurat <- RenameIdents(MOPA_Myocytes_seurat,"C22"="C2223","C23"="C2223")
MOPA_Myocytes_seurat <- RenameIdents(MOPA_Myocytes_seurat,"C14"="C1314","C13"="C1314")
DimPlot(MOPA_Myocytes_seurat,label = T,repel = T)
MOPA_Myocytes$Peak_Clusters3.2 <- as.character(MOPA_Myocytes_seurat@active.ident)
#
#"C15","C22","C23","C21",
trajectory <- c("C14","C13","C11","C12","C9","C8","C5","C7","C3","C2","C18","C19","C24","C23","C22")
trajectory
MOPA_Myocytes <- addTrajectory(
  ArchRProj = MOPA_Myocytes,
  reducedDims = "Peak_IterativeLSI",
  name = "MBMT", 
  groupBy = "Peak_Clusters3",
  trajectory = trajectory, 
  embedding = "Peak_UMAP5", 
  force = TRUE
)
p <- plotTrajectory(MOPA_Myocytes, embedding = "Peak_UMAP5",trajectory = "MBMT", colorBy = "cellColData", name = "MBMT")
p[[1]]
#
trajectory <- c("C1314","C11","C12","C89","C2357","C18","C19","C24","C23","C22")
trajectory
MOPA_Myocytes <- addTrajectory(
  ArchRProj = MOPA_Myocytes,
  reducedDims = "Peak_IterativeLSI",
  name = "MBMT.2", 
  groupBy = "Peak_Clusters3.2",
  trajectory = trajectory, 
  embedding = "Peak_UMAP5", 
  force = TRUE
)
p <- plotTrajectory(MOPA_Myocytes, embedding = "Peak_UMAP5",trajectory = "MBMT.2", colorBy = "cellColData", name = "MBMT.2")
p[[1]]
#
#######
MOPA_Myocytes.pseudo <- as.data.frame(cbind(MOPA_Myocytes$MBMT,MOPA_Myocytes$cellNames))
head(MOPA_Myocytes.pseudo)
x.rank <- na.omit(MOPA_Myocytes.pseudo)
head(x.rank)
dim(x.rank)
hist(as.numeric(x.rank$V1),breaks = 100)
#
t = x.rank[order(x.rank$V1),]
head(t)
#
t$T100 <- paste("C",c(rep(1:100,each=70),rep(100,8)),sep = "")
head(t)
rownames(t) <- t$V2
#
MOPA_MyocytesT100 <- MOPA_Myocytes[t$V2,]
t.tmp <- t[MOPA_MyocytesT100$cellNames,]
MOPA_MyocytesT100$T100 <- t.tmp$T100
#
plotEmbedding(ArchRProj = MOPA_MyocytesT100, colorBy = "cellColData", name = "T100", embedding = "Peak_UMAP5")+NoLegend()
#
MOPA_MyocytesT100 <- addPeak2GeneLinks(
  ArchRProj = MOPA_MyocytesT100,
  reducedDims = "Peak_IterativeLSI",
  useMatrix = "GeneScoreMatrix"#GeneIntegrationMatrix
)
p2g <- getPeak2GeneLinks(
  ArchRProj = MOPA_MyocytesT100,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)
p2g
#
p <- plotPeak2GeneHeatmap(ArchRProj = MOPA_MyocytesT100, groupBy = "T100",k = 25)
p
#
p.P2GT100m <- plotPeak2GeneHeatmap(ArchRProj = MOPA_MyocytesT100, groupBy = "T100",returnMatrices=T)
p.P2GT100m
#
#
#
trajPM  <- getTrajectory(ArchRProj = MOPA_MyocytesT100, name = "MBMT", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"),returnMatrix=F)
p4
####
p4.m <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"),returnMatrix=T)
#
p4.m
dim(p4.m)
pheatmap(p4.m,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F)
#
p1 <- plotTrajectory(MOPA_Myocytes, trajectory = "MBMT", colorBy = "GeneScoreMatrix", name = "Myod1", continuousSet = "horizonExtra", embedding = "Peak_UMAP5")
p2 <- plotTrajectory(MOPA_Myocytes, trajectory = "MBMT", colorBy = "GeneIntegrationMatrix", name = "Myod1", continuousSet = "blueYellow", embedding = "Peak_UMAP5")
ggAlignPlots(p1[[1]], p2[[1]], type = "h")
ggAlignPlots(p1[[2]], p2[[2]], type = "h")
#
#
#
#
#
MOPA_Myocytes <- addPeak2GeneLinks(
  ArchRProj = MOPA_Myocytes,
  reducedDims = "Peak_IterativeLSI",
  useMatrix = "GeneScoreMatrix"#GeneIntegrationMatrix
)
p2g <- getPeak2GeneLinks(
  ArchRProj = MOPA_Myocytes,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)
p2g
#
p <- plotPeak2GeneHeatmap(ArchRProj = MOPA_Myocytes, groupBy = "T100",k = 25)
p
#
p.p2g.total.m <- plotPeak2GeneHeatmap(ArchRProj = MOPA_Myocytes,nPlot = 70000, groupBy = "Peak_Clusters3",k = 25,returnMatrices =T)
#p.p2g.total.m
p.p2g.total.m$ATAC
dim(p.p2g.total.m$ATAC$matrix)
dim(p.p2g.total.m$RNA$matrix)
max(p.p2g.total.m$RNA$matrix)
min(p.p2g.total.m$RNA$matrix)
#
pheatmap(p.p2g.total.m$RNA$matrix[1:1000,],cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F)
#
#
##########
trajPM  <- getTrajectory(ArchRProj = MOPA_Myocytes, name = "MBMT", useMatrix = "PeakMatrix", log2Norm = TRUE)
p4 <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"),returnMatrix=F)
p4
####
p4.peak.m <- plotTrajectoryHeatmap(trajPM, pal = paletteContinuous(set = "solarExtra"),returnMatrix=T)
p4.peak.m
dim(p4.peak.m);max(p4.peak.m);min(p4.peak.m)
head(p4.peak.m)
pheatmap(p4.peak.m,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F)
#################
trajGSM <- getTrajectory(ArchRProj = MOPA_Myocytes, name = "MBMT", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p2 <- trajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
p2
#
p2.GSM.m <- trajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"),returnMatrix=T)
dim(p2.GSM.m);max(p2.GSM.m);min(p2.GSM.m)
head(p2.GSM.m)
pheatmap(p2.GSM.m,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F)
###############
trajGIM <- getTrajectory(ArchRProj = MOPA_Myocytes, name = "MBMT", useMatrix = "GeneIntegrationMatrix", log2Norm = FALSE)
p3 <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"))
p3 
#
p3.GIM.m <- plotTrajectoryHeatmap(trajGIM,  pal = paletteContinuous(set = "blueYellow"),returnMatrix=T)
dim(p3.GIM.m);max(p3.GIM.m);min(p3.GIM.m)
pheatmap(p3.GIM.m,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F)
#
#
#
p2geneDF <- metadata(MOPA_Myocytes@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
p2geneDF
p2geneDF.tmp <- as.data.frame(p2geneDF)
p2geneDF.tmp$peakName
#
p2geneDF.tmp$peak <- paste(paste(unlist(strsplit(p2geneDF.tmp$peakName,"_"))[seq(1,3*nrow(p2geneDF.tmp),3)],
                                 unlist(strsplit(p2geneDF.tmp$peakName,"_"))[seq(2,3*nrow(p2geneDF.tmp),3)],sep = ":"),
                           unlist(strsplit(p2geneDF.tmp$peakName,"_"))[seq(3,3*nrow(p2geneDF.tmp),3)],sep = "_")
p2geneDF.tmp$gene <- paste(unlist(strsplit(p2geneDF.tmp$peakName,"_"))[seq(1,3*nrow(p2geneDF.tmp),3)],p2geneDF.tmp$geneName,sep = ":")
#
dim(p2geneDF.tmp)
p2geneDF.tmp <- subset(p2geneDF.tmp,Correlation > 0.45)
dim(p2geneDF.tmp)
p2geneDF.tmp <- subset(p2geneDF.tmp,VarQATAC > 0.2)
p2geneDF.tmp <- subset(p2geneDF.tmp,VarQRNA > 0.2)
dim(p2geneDF.tmp)
head(p2geneDF.tmp)
#
length(unique(p2geneDF.tmp$peak))
p2gene.peak2gene <- as.data.frame(table(p2geneDF.tmp$peak))
dim(p2gene.peak2gene)
hist(table(p2geneDF.tmp$peak),breaks = 100)
#
head(p2gene.peak2gene)
p2gene.peak2gene <- subset(p2gene.peak2gene,Freq < 4)
dim(p2gene.peak2gene)
#
geneindex1 <- c()
geneindex2 <- c()
for (i in 1:nrow(p2gene.peak2gene)) {
  geneindex1 <- which(p2geneDF.tmp$peak==p2gene.peak2gene$Var1[i])
  geneindex2 <- c(geneindex2,geneindex1)
}
dim(p2geneDF.tmp)
p2geneDF.tmp <- p2geneDF.tmp[geneindex2,]
dim(p2geneDF.tmp)
#
length(rownames(p4.peak.m))
length(unique(rownames(p4.peak.m)))
#
p4.peak.m.1 <- p4.peak.m
dim(p4.peak.m.1)
p4.peak.m.2 <- p4.peak.m.1[intersect(rownames(p4.peak.m.1),p2geneDF.tmp$peak),]
dim(p4.peak.m.2)
pheatmap(p4.peak.m.2,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F)
#
head(p2.GSM.m)
p2.GSM.m.1 <- p2.GSM.m
dim(p2.GSM.m.1)
p2.GSM.m.2 <- p2.GSM.m.1[intersect(rownames(p2.GSM.m.1),p2geneDF.tmp$gene),]
dim(p2.GSM.m.2)
head((p2.GSM.m.2))
pheatmap(p2.GSM.m.2,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F)
#
#
#### 
geneindex1 <- c()
geneindex2 <- c()
dim(p2geneDF.tmp)
p2geneDF.tmp.peakrm <- setdiff(unique(p2geneDF.tmp$peak),rownames(p4.peak.m.2))
length(p2geneDF.tmp.peakrm)
for (i in 1:length(p2geneDF.tmp.peakrm)) {
  geneindex1 <- which(p2geneDF.tmp$peak==p2geneDF.tmp.peakrm[i])
  geneindex2 <- c(geneindex2,geneindex1)
}
length(geneindex2)
head(geneindex2)
p2geneDF.tmp1 <- p2geneDF.tmp[-geneindex2,]
#
geneindex1 <- c()
geneindex2 <- c()
dim(p2geneDF.tmp1)
head(p2geneDF.tmp1)
p2geneDF.tmp.generm <- setdiff(unique(p2geneDF.tmp1$gene),rownames(p2.GSM.m.2))
length(p2geneDF.tmp.generm)
for (i in 1:length(p2geneDF.tmp.generm)) {
  geneindex1 <- which(p2geneDF.tmp1$gene==p2geneDF.tmp.generm[i])
  geneindex2 <- c(geneindex2,geneindex1)
}
length(geneindex2)
head(geneindex2)
p2geneDF.tmp2 <- p2geneDF.tmp1[-geneindex2,]
dim(p2geneDF.tmp2)
head(p2geneDF.tmp2)
##
#
#
geneindex1 <- c()
geneindex2 <- c()
geneindex3 <- c()
nrow(p4.peak.m.2)
for (i in 1:nrow(p4.peak.m.2)) {
  geneindex1 <- which(p2geneDF.tmp2$peak==rownames(p4.peak.m.2)[i])
  geneindex2 <- c(geneindex2,geneindex1)
  geneindex3 <- c(geneindex3,length(geneindex1))
}
#
length(geneindex2)
length(geneindex3)
head(geneindex2)
head(geneindex3)
#
p2.GSM.m.3 <- p2.GSM.m.2[p2geneDF.tmp2[geneindex2,]$gene,]
dim(p2.GSM.m.3)
head(p2.GSM.m.3)
tail(rownames(p2.GSM.m.3))
length(unique(rownames(p2.GSM.m.3)))
length(rownames(p2.GSM.m.3))
pheatmap(p2.GSM.m.3,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,
         color = paletteContinuous(set = "solarExtra"))
pheatmap(p4.peak.m.2,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F)
#
dim(p4.peak.m.2)
p4.peak.m.3 <- p4.peak.m.2[intersect(rownames(p4.peak.m.2),p2geneDF.tmp2$peak),]
dim(p4.peak.m.3)
pheatmap(p4.peak.m.3,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F)
#
geneindex4 <- as.data.frame(geneindex3)
head(geneindex4)
geneindex4 <- geneindex4[-which(geneindex4$geneindex3==0),]
head(geneindex4)
sum(geneindex4)
head(rep(1:nrow(p4.peak.m.3),geneindex4))
p4.peak.m.4 <- p4.peak.m.3[rep(1:nrow(p4.peak.m.3),geneindex4),]
dim(p4.peak.m.4)
pheatmap(p4.peak.m.4,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F)
pheatmap(p4.peak.m.4,cluster_rows = T,cluster_cols = F,show_rownames = F,show_colnames = F)
#
pheatmap(p4.peak.m.4,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,
         color = paletteContinuous(set = "horizonExtra"))
#
#horizonExtra paletteContinuous(set = "solarExtra")
p4.peak.m.5 <- p4.peak.m.4
head(p4.peak.m.5)
rownames(p4.peak.m.5) <- paste(1:7267,rownames(p4.peak.m.5),sep = "#")
km_result1 <- kmeans(p4.peak.m.5,4,nstart = 24)
head(km_result1$cluster)
pheatmap(cbind(km_result1$cluster,km_result1$cluster),cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F)
#
km_result1.tmp <- as.data.frame(km_result1$cluster)
km_result1.tmp$rank <- 1:7267
dim(km_result1.tmp)
head(km_result1.tmp)
km_result1.tmp.m1 <- subset(km_result1.tmp,km_result1$cluster=="2")
km_result1.tmp.m2 <- subset(km_result1.tmp,km_result1$cluster=="3")
km_result1.tmp.m3 <- subset(km_result1.tmp,km_result1$cluster=="4")
km_result1.tmp.m4 <- subset(km_result1.tmp,km_result1$cluster=="1")
#
#
#
pheatmap(p2.GSM.m.3[km_result1.tmp.m1$rank,],use_raster=T,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,
         color = paletteContinuous(set = "solarExtra"))
pheatmap(p2.GSM.m.3[km_result1.tmp.m2$rank,],use_raster=T,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,
         color = paletteContinuous(set = "solarExtra"))
pheatmap(p2.GSM.m.3[km_result1.tmp.m3$rank,],use_raster=T,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,
         color = paletteContinuous(set = "solarExtra"))
pheatmap(p2.GSM.m.3[km_result1.tmp.m4$rank,],use_raster=T,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,
         color = paletteContinuous(set = "solarExtra"))
#
pheatmap(p4.peak.m.4[km_result1.tmp.m1$rank,],use_raster=T,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,
         color = paletteContinuous(set = "horizonExtra"))
pheatmap(p4.peak.m.4[km_result1.tmp.m2$rank,],use_raster=T,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,
         color = paletteContinuous(set = "horizonExtra"))
pheatmap(p4.peak.m.4[km_result1.tmp.m3$rank,],use_raster=T,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,
         color = paletteContinuous(set = "horizonExtra"))
pheatmap(p4.peak.m.4[km_result1.tmp.m4$rank,],use_raster=T,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,
         color = paletteContinuous(set = "horizonExtra"))
#
km1.peaks <- as.character(unique(rownames(p4.peak.m.4[km_result1.tmp.m1$rank,])))
km2.peaks <- as.character(unique(rownames(p4.peak.m.4[km_result1.tmp.m2$rank,])))
km3.peaks <- as.character(unique(rownames(p4.peak.m.4[km_result1.tmp.m3$rank,])))
km4.peaks <- as.character(unique(rownames(p4.peak.m.4[km_result1.tmp.m4$rank,])))
#
km1.peaks.peak <- cbind(unlist(strsplit(km1.peaks,":"))[seq(1,2*length(km1.peaks),2)],
                        unlist(strsplit(unlist(strsplit(km1.peaks,":"))[seq(2,2*length(km1.peaks),2)],"_"))[seq(1,2*length(km1.peaks),2)],
                        unlist(strsplit(unlist(strsplit(km1.peaks,":"))[seq(2,2*length(km1.peaks),2)],"_"))[seq(1,2*length(km1.peaks),2)])
km2.peaks.peak <- cbind(unlist(strsplit(km2.peaks,":"))[seq(1,2*length(km2.peaks),2)],
                        unlist(strsplit(unlist(strsplit(km2.peaks,":"))[seq(2,2*length(km2.peaks),2)],"_"))[seq(1,2*length(km2.peaks),2)],
                        unlist(strsplit(unlist(strsplit(km2.peaks,":"))[seq(2,2*length(km2.peaks),2)],"_"))[seq(1,2*length(km2.peaks),2)])
km3.peaks.peak <- cbind(unlist(strsplit(km3.peaks,":"))[seq(1,2*length(km3.peaks),2)],
                        unlist(strsplit(unlist(strsplit(km3.peaks,":"))[seq(2,2*length(km3.peaks),2)],"_"))[seq(1,2*length(km3.peaks),2)],
                        unlist(strsplit(unlist(strsplit(km3.peaks,":"))[seq(2,2*length(km3.peaks),2)],"_"))[seq(1,2*length(km3.peaks),2)])
km4.peaks.peak <- cbind(unlist(strsplit(km4.peaks,":"))[seq(1,2*length(km4.peaks),2)],
                        unlist(strsplit(unlist(strsplit(km4.peaks,":"))[seq(2,2*length(km4.peaks),2)],"_"))[seq(1,2*length(km4.peaks),2)],
                        unlist(strsplit(unlist(strsplit(km4.peaks,":"))[seq(2,2*length(km4.peaks),2)],"_"))[seq(1,2*length(km4.peaks),2)])
dim(km1.peaks.peak)
dim(km2.peaks.peak)
dim(km3.peaks.peak)
dim(km4.peaks.peak)
#
write.table(km1.peaks.peak,file = paste("p2g.module/km1.peaks.peak.csv",sep = ""),quote=F, sep = "\t",row.names = F,col.names = F)
write.table(km2.peaks.peak,file = paste("p2g.module/km2.peaks.peak.csv",sep = ""),quote=F, sep = "\t",row.names = F,col.names = F)
write.table(km3.peaks.peak,file = paste("p2g.module/km3.peaks.peak.csv",sep = ""),quote=F, sep = "\t",row.names = F,col.names = F)
write.table(km4.peaks.peak,file = paste("p2g.module/km4.peaks.peak.csv",sep = ""),quote=F, sep = "\t",row.names = F,col.names = F)
#
#as.character(unique(rownames(p2.GSM.m.3[km_result1.tmp.m1$rank,])))
km1.genes <- as.character(unique(rownames(p2.GSM.m.3[km_result1.tmp.m1$rank,])))
km2.genes <- as.character(unique(rownames(p2.GSM.m.3[km_result1.tmp.m2$rank,])))
km3.genes <- as.character(unique(rownames(p2.GSM.m.3[km_result1.tmp.m3$rank,])))
km4.genes <- as.character(unique(rownames(p2.GSM.m.3[km_result1.tmp.m4$rank,])))
#
km1.genes.gene <- unlist(strsplit(km1.genes,":"))[seq(2,2*length(km1.genes),2)]
km2.genes.gene <- unlist(strsplit(km2.genes,":"))[seq(2,2*length(km2.genes),2)]
km3.genes.gene <- unlist(strsplit(km3.genes,":"))[seq(2,2*length(km3.genes),2)]
km4.genes.gene <- unlist(strsplit(km4.genes,":"))[seq(2,2*length(km4.genes),2)]
#
write.table(km1.genes.gene,file = paste("p2g.module/km1.genes.gene.csv",sep = ""),quote=F, sep = "\t",row.names = F,col.names = F)
write.table(km2.genes.gene,file = paste("p2g.module/km2.genes.gene.csv",sep = ""),quote=F, sep = "\t",row.names = F,col.names = F)
write.table(km3.genes.gene,file = paste("p2g.module/km3.genes.gene.csv",sep = ""),quote=F, sep = "\t",row.names = F,col.names = F)
write.table(km4.genes.gene,file = paste("p2g.module/km4.genes.gene.csv",sep = ""),quote=F, sep = "\t",row.names = F,col.names = F)
#
#
write.csv(p2.GSM.m.3[km_result1.tmp.m1$rank,],file = "./p2g.module/p2g.matrix/P2G.heatmap.gene.module1.csv")
write.csv(p2.GSM.m.3[km_result1.tmp.m2$rank,],file = "./p2g.module/p2g.matrix/P2G.heatmap.gene.module2.csv")
write.csv(p2.GSM.m.3[km_result1.tmp.m3$rank,],file = "./p2g.module/p2g.matrix/P2G.heatmap.gene.module3.csv")
write.csv(p2.GSM.m.3[km_result1.tmp.m4$rank,],file = "./p2g.module/p2g.matrix/P2G.heatmap.gene.module4.csv")
#
#
write.csv(p4.peak.m.4[km_result1.tmp.m1$rank,],file = "./p2g.module/p2g.matrix/P2G.heatmap.peak.module1.csv")
write.csv(p4.peak.m.4[km_result1.tmp.m2$rank,],file = "./p2g.module/p2g.matrix/P2G.heatmap.peak.module2.csv")
write.csv(p4.peak.m.4[km_result1.tmp.m3$rank,],file = "./p2g.module/p2g.matrix/P2G.heatmap.peak.module3.csv")
write.csv(p4.peak.m.4[km_result1.tmp.m4$rank,],file = "./p2g.module/p2g.matrix/P2G.heatmap.peak.module4.csv")
#
#
trajGSM <- getTrajectory(ArchRProj = MOPA_Myocytes, name = "MBMT", useMatrix = "GeneScoreMatrix", log2Norm = TRUE)
p1.GSM <- plotTrajectoryHeatmap(trajGSM,  pal = paletteContinuous(set = "horizonExtra"))
p1.GSM
#
trajMM  <- getTrajectory(ArchRProj = MOPA_Myocytes, name = "MBMT", useMatrix = "MotifMatrix", log2Norm = FALSE)
p2.MM <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"))
p2.MM
p2.MM.m  <- plotTrajectoryHeatmap(trajMM, pal = paletteContinuous(set = "solarExtra"),returnMatrix =T)
dim(p2.MM.m)
#
corGSM_MM <- correlateTrajectories(trajGSM, trajMM,corCutOff = 0.55)
corGSM_MM[[1]]
trajGSM2 <- trajGSM[corGSM_MM[[1]]$name1, ]
trajMM2 <- trajMM[corGSM_MM[[1]]$name2, ]
trajCombined <- trajGSM2
assay(trajCombined) <- t(apply(assay(trajGSM2), 1, scale)) + t(apply(assay(trajMM2), 1, scale))
combinedMat <- plotTrajectoryHeatmap(trajCombined, returnMat = TRUE, varCutOff = 0)
rowOrder <- match(rownames(combinedMat), rownames(trajGSM2))
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "solarExtra"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "blueYellow"), varCutOff = 0, rowOrder = rowOrder)
ht1 + ht2
#
#
#
#
#
plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "cellColData", name = "Timepoint", embedding = "Peak_UMAP5")
MOPA_Myocytes_seurat
MOPA_Myocytes_seurat@active.ident <- as.factor(MOPA_Myocytes_seurat$Timepoint)
DimPlot(MOPA_Myocytes_seurat)
#
MOPA_Myocytes_seurat.metadata.myocytes <- subset(MOPA_Myocytes_seurat.metadata,emID.F2=="Myocytes")
dim(MOPA_Myocytes_seurat.metadata.myocytes)
head(MOPA_Myocytes_seurat.metadata.myocytes)
MOPA_Myocytes_seurat.metadata.myocytes$cells2 <- paste("Myocytes",MOPA_Myocytes_seurat.metadata.myocytes$cells,sep = "#")
head(MOPA_Myocytes_seurat)
MOPA_Myocytes_seurat.filter <- subset(MOPA_Myocytes_seurat,cells=MOPA_Myocytes_seurat.metadata.myocytes$cells2)
MOPA_Myocytes_seurat.filter
#
DimPlot(MOPA_Myocytes_seurat.filter)
Timepoint.color <- c("#98304E","#DBCCE0","#48D1CC","#6861A9")
DimPlot(MOPA_Myocytes_seurat,label = T,repel = T,cols = Timepoint.color )
p <- DimPlot(MOPA_Myocytes_seurat.filter,label =F,repel = T,cols =Timepoint.color,label.size = 6,pt.size = 0.0001)+NoLegend()+
  #ggtitle("hpf72h_Batch.ID")+
  theme(plot.title = element_text(color="black", size=14))+
  theme(plot.title = element_text(hjust = 0.5))+
  #labs(x = "tSNE_1", y = "tSNE_2")+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("Myocytes.timepoint.png",plot=p,device="png",dpi=400,units = "cm",width = 13,height = 13)
p <- DimPlot(MOPA_Myocytes_seurat.filter,label =F,repel = T,cols =c("gray","gray","gray","gray"),label.size = 6,pt.size = 0.0001)+NoLegend()+
  #ggtitle("hpf72h_Batch.ID")+
  theme(plot.title = element_text(color="black", size=14))+
  theme(plot.title = element_text(hjust = 0.5))+
  #labs(x = "tSNE_1", y = "tSNE_2")+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("Myocytes.gray2.png",plot=p,device="png",dpi=400,units = "cm",width = 13,height = 13)
#
#
p <- plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Myf5"),quantCut = c(0.35, 0.99), pal = paletteContinuous(set = "solarExtra"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./feature/Myf5.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
#
p <- plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Pax3"), pal = paletteContinuous(set = "solarExtra"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./feature/Pax3.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
#
p <- plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Pax5"),quantCut = c(0.1, 0.99), pal = paletteContinuous(set = "solarExtra"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./feature/Pax5.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
#
p <- plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Pax7"),quantCut = c(0.1, 0.99), pal = paletteContinuous(set = "solarExtra"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./feature/Pax7.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
#
p <- plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Myh3"), pal = paletteContinuous(set = "solarExtra"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./feature/Myh3.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
#
p <- plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Myod1"), pal = paletteContinuous(set = "solarExtra"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./feature/Myod1.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
#
p <- plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Myog"), pal = paletteContinuous(set = "solarExtra"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./feature/Myog.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
#
p <- plotEmbedding(ArchRProj = MOPA_Myocytes, colorBy = "GeneScoreMatrix", name = c("Lhx2"), pal = paletteContinuous(set = "solarExtra"), embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./feature/Lhx2.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
#
##
p <- plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM1Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./feature/mM1.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
p <- plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM2Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./feature/mM2.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
p <- plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM3Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./feature/mM3.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
p <- plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM4Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./feature/mM4.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
p <- plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM5Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./feature/mM5.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
p <- plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM6Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./feature/mM6.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
p <- plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM7Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./feature/mM7.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
p <- plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM8Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./feature/mM8.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
p <- plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM9Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./feature/mM9.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
p <- plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM10Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./feature/mM10.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
p <- plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM11Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./feature/mM11.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
p <- plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM12Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./feature/mM12.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
p <- plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM13Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./feature/mM13.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
p <- plotEmbedding(ArchRProj = MOPA_Myocytes,colorBy = "cellColData",size = 0.0001, name = "Module.mM14Score", embedding = "Peak_UMAP5",imputeWeights = getImputeWeights(MOPA_Myocytes))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("./feature/mM14.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
#
trajectory <- c("C14","C10","C21","C23","C22")
trajectory
MOPA_Myocytes <- addTrajectory(
  ArchRProj = MOPA_Myocytes,
  reducedDims = "Peak_IterativeLSI",
  name = "MBMT.Myf5", 
  groupBy = "Peak_Clusters3",
  trajectory = trajectory, 
  embedding = "Peak_UMAP5", 
  force = TRUE
)
p <- plotTrajectory(MOPA_Myocytes, embedding = "Peak_UMAP5",trajectory = "MBMT.Myf5", colorBy = "cellColData", name = "MBMT.Myf5")
p[[1]]
#
trajectory <- c("C17","C4","C7","C3","C2","C5","C18","C19","C24","C23","C22")
trajectory
MOPA_Myocytes <- addTrajectory(
  ArchRProj = MOPA_Myocytes,
  reducedDims = "Peak_IterativeLSI",
  name = "MBMT.Lhx2", 
  groupBy = "Peak_Clusters3",
  trajectory = trajectory, 
  embedding = "Peak_UMAP5", 
  force = TRUE
)
p <- plotTrajectory(MOPA_Myocytes, embedding = "Peak_UMAP5",trajectory = "MBMT.Lhx2", colorBy = "cellColData", name = "MBMT.Lhx2")
p[[1]]
#
trajectory <- c("C15","C6","C7","C3","C2","C5","C18","C19","C24","C23","C22")
trajectory
MOPA_Myocytes <- addTrajectory(
  ArchRProj = MOPA_Myocytes,
  reducedDims = "Peak_IterativeLSI",
  name = "MBMT.der", 
  groupBy = "Peak_Clusters3",
  trajectory = trajectory, 
  embedding = "Peak_UMAP5", 
  force = TRUE
)
p <- plotTrajectory(MOPA_Myocytes, embedding = "Peak_UMAP5",trajectory = "MBMT.der", colorBy = "cellColData", name = "MBMT.der")
p[[1]]
#
ht1 <- plotTrajectoryHeatmap(trajGSM2,  pal = paletteContinuous(set = "solarExtra"),  varCutOff = 0, rowOrder = rowOrder)
ht2 <- plotTrajectoryHeatmap(trajMM2, pal = paletteContinuous(set = "blueYellow"), varCutOff = 0, rowOrder = rowOrder)
ht1 + ht2
#
pheatmap(ht1@ht_list$GeneScoreMatrix@matrix[-2,],use_raster=T,color = paletteContinuous(set = "solarExtra"),cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F)
pheatmap(ht2@ht_list$MotifMatrix@matrix[-2,],use_raster=T,color = paletteContinuous(set = "blueYellow"),cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F)
#
head(ht1@ht_list$GeneScoreMatrix@matrix[-2,])
x <- ht1@ht_list$GeneScoreMatrix@matrix[-2,]
dim(x)
colnames(x) <- paste("T",1:100,sep = "")
write.csv(x,file = "myocyte.pseudotime.motifmatrix.csv")
#
x <- ht2@ht_list$MotifMatrix@matrix[-2,]
dim(x)
head(x)
colnames(x) <- paste("T",1:100,sep = "")
write.csv(x,file = "myocyte.pseudotime.genescorematrix.csv")
#
#
library(ComplexHeatmap)
Heatmap(cbind(c(1:100),c(1:100)),use_raster=T,color = c("#D6EACB","#4699C8"),cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F)
Heatmap(cbind(c(1:100),c(1:100)),use_raster=T,col = c("#D6EACB","#4699C8"),cluster_rows =F,cluster_columns =F)
Heatmap(c(1:100),use_raster=T,col = c("#D6EACB","#4699C8"),cluster_rows =F,cluster_columns =F)
#
#)
#
#
#
MOPA_Myocytes$MBMT
MOPA_MyocytesT100$MBMT
#
MOPA_Myocytes.MBMT <- as.data.frame(cbind(MOPA_Myocytes$cellNames,MOPA_Myocytes$MBMT))
colnames(MOPA_Myocytes.MBMT) <- c("cells","pseudovalues")
dim(MOPA_Myocytes.MBMT)
head(MOPA_Myocytes.MBMT)
MOPA_Myocytes.MBMT <- na.omit(MOPA_Myocytes.MBMT)
dim(MOPA_Myocytes.MBMT)
head(MOPA_Myocytes.MBMT)
#
MOPA_Myocytes.MBMT$pseudo2 <- MOPA_Myocytes.MBMT$pseudovalues
head(MOPA_Myocytes.MBMT)
MOPA_Myocytes.MBMT[MOPA_Myocytes.MBMT$pseudovalues > 80 ,3] <- "T5"
MOPA_Myocytes.MBMT[MOPA_Myocytes.MBMT$pseudovalues > 60 & MOPA_Myocytes.MBMT$pseudovalues <= 80,3] <- "T4"
MOPA_Myocytes.MBMT[MOPA_Myocytes.MBMT$pseudovalues > 40 & MOPA_Myocytes.MBMT$pseudovalues <= 60,3] <- "T3"
MOPA_Myocytes.MBMT[MOPA_Myocytes.MBMT$pseudovalues > 20 & MOPA_Myocytes.MBMT$pseudovalues <= 40,3] <- "T2"
MOPA_Myocytes.MBMT[MOPA_Myocytes.MBMT$pseudovalues >= 0 & MOPA_Myocytes.MBMT$pseudovalues <= 20,3] <- "T1"
head(MOPA_Myocytes.MBMT)
#
MOPA_Myocytes.MBMT.ArchR <- MOPA_Myocytes[MOPA_Myocytes.MBMT$V1,]
MOPA_Myocytes.MBMT.ArchR
#
MOPA_Myocytes.MBMT.ArchR$pseudoT15 <- as.character(projEu.MBMT.meatadata$V3)
#
addArchRThreads(threads = 1) 
MOPA_Myocytes.MBMT.ArchR <- addGroupCoverages(ArchRProj = MOPA_Myocytes.MBMT.ArchR, groupBy = "pseudoT15")
#
p <- plotBrowserTrack(ArchRProj = MOPA_Myocytes.MBMT.ArchR, groupBy = "pseudoT15", geneSymbol = "Myog", 
                      upstream = 70000,downstream = 70000,ylim = c(0.01,1), loops = getPeak2GeneLinks(MOPA_Myocytes.MBMT.ArchR))
grid::grid.newpage()
grid::grid.draw(p$Myog)
#
MOPA_Myocytes.MBMT.ArchR <- addMotifAnnotations(ArchRProj = MOPA_Myocytes.MBMT.ArchR, motifSet = "cisbp", name = "Motif")
MOPA_Myocytes.MBMT.ArchR <- addBgdPeaks(MOPA_Myocytes.MBMT.ArchR)
motifPositions <- getPositions(MOPA_Myocytes.MBMT.ArchR)
#
motifs <- c("1","2","3","4","5","6","7","8","9")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs
markerMotifs <- unique(markerMotifs)
length(markerMotifs)
markerMotifs
#
seFoot <- getFootprints(
  ArchRProj = MOPA_Myocytes.MBMT.ArchR, 
  positions = motifPositions[markerMotifs[751:884]], 
  groupBy = "pseudoT15"
)
#
plotFootprints(
  seFoot = seFoot,
  ArchRProj = MOPA_Myocytes.MBMT.ArchR, 
  normMethod = "Subtract",
  plotName = "Footprints-6",
  addDOC = FALSE,
  smoothWindow = 5
)
#
save.image("MOPA.Myocytes.RData")
#