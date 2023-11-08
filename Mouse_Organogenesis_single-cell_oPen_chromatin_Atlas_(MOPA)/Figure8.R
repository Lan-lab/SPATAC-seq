# For the analysis in the Figure 8
setwd("/home/sunkeyong/MOPA_project/Fetal_adult")
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
library(tidyverse)
library(ggExtra)
library(ComplexHeatmap)
library("RColorBrewer")
library("circlize")
require(org.Mm.eg.db)
#
addArchRThreads(threads = 30) 
#
load("Mouse.Fetal.Adult.RData")
#
table(Fetal.Adult.seuratOBO$MAE.ID)
Fetal.Adult.seuratOBO@active.ident <- as.factor(Fetal.Adult.seuratOBO$MAE.ID)
Fetal.Adult.seuratOBO.Cardio <- subset(Fetal.Adult.seuratOBO,idents=c("Adult.Cardiomyocytes","Fetal.Cardiac muscle lineages"))
Fetal.Adult.seuratOBO.Cardio
#
Fetal.Adult.ArchROBO.tmp <- Fetal.Adult.ArchROBO[Cells(Fetal.Adult.seuratOBO.Cardio),]
Fetal.Adult.ArchROBO.tmp
table(Fetal.Adult.ArchROBO.tmp$Sample)
#
saveArchRProject(ArchRProj = Fetal.Adult.ArchROBO.tmp, outputDirectory = "Save-ProjHeme_Cardiomyocytes", load = T,dropCells = T)
#
table(Fetal.Adult.seuratOBO$MAE.ID)
Fetal.Adult.seuratOBO@active.ident <- as.factor(Fetal.Adult.seuratOBO$MAE.ID)
Fetal.Adult.seuratOBO.ExIn <- subset(Fetal.Adult.seuratOBO,idents=c("Adult.Ex. neurons SCPN","Adult.Ex or in neurons","Adult.Inhibitory neurons",
                                                        "Adult.Ex or in neurons","Fetal.Inhibitory neurons","Fetal.Inhibitory neuron progenitors","Fetal.Excitatory neurons"))
Fetal.Adult.seuratOBO.ExIn
#
Fetal.Adult.ArchROBO.tmp <- Fetal.Adult.ArchROBO[Cells(Fetal.Adult.seuratOBO.ExIn),]
Fetal.Adult.ArchROBO.tmp
table(Fetal.Adult.ArchROBO.tmp$Sample)
#
saveArchRProject(ArchRProj = Fetal.Adult.ArchROBO.tmp, outputDirectory = "Save-ProjHeme_neuron", load = T,dropCells = T)
#
table(Fetal.Adult.seuratOBO$MAE.ID)
Fetal.Adult.seuratOBO@active.ident <- as.factor(Fetal.Adult.seuratOBO$MAE.ID)
Fetal.Adult.seuratOBO.Hepatocyte <- subset(Fetal.Adult.seuratOBO,idents=c("Adult.Hepatocytes","Fetal.Hepatocytes"))
Fetal.Adult.seuratOBO.Hepatocyte
#
#
Fetal.Adult.ArchROBO.tmp <- Fetal.Adult.ArchROBO[Cells(Fetal.Adult.seuratOBO.Hepatocyte),]
Fetal.Adult.ArchROBO.tmp
table(Fetal.Adult.ArchROBO.tmp$Sample)
saveArchRProject(ArchRProj = Fetal.Adult.ArchROBO.tmp, outputDirectory = "Save-ProjHeme_Hepatocytes", load = T,dropCells = T)
#
########
table(Fetal.Adult.seuratOBO$MAE.ID)
Fetal.Adult.seuratOBO@active.ident <- as.factor(Fetal.Adult.seuratOBO$MAE.ID)
Fetal.Adult.seuratOBO.Oligo <- subset(Fetal.Adult.seuratOBO,idents=c("Adult.Oligodendrocytes","Fetal.Oligodendrocyte Progenitors","Fetal.Premature oligodendrocyte"))
Fetal.Adult.seuratOBO.Oligo
#
Fetal.Adult.ArchROBO.tmp <- Fetal.Adult.ArchROBO[Cells(Fetal.Adult.seuratOBO.Oligo),]
Fetal.Adult.ArchROBO.tmp
table(Fetal.Adult.ArchROBO.tmp$Sample)
#
saveArchRProject(ArchRProj = Fetal.Adult.ArchROBO.tmp, outputDirectory = "Save-ProjHeme_Oligodendrocytes", load = T,dropCells = T)
#
table(Fetal.Adult.seuratOBO$MAE.ID)
Fetal.Adult.seuratOBO@active.ident <- as.factor(Fetal.Adult.seuratOBO$MAE.ID)
Fetal.Adult.seuratOBO.granule <- subset(Fetal.Adult.seuratOBO,idents=c("Adult.Cerebellar granule cells","Fetal.Granule neurons"))
Fetal.Adult.seuratOBO.granule
#
#
Fetal.Adult.ArchROBO.tmp <- Fetal.Adult.ArchROBO[Cells(Fetal.Adult.seuratOBO.granule),]
Fetal.Adult.ArchROBO.tmp
table(Fetal.Adult.ArchROBO.tmp$Sample)
#
saveArchRProject(ArchRProj = Fetal.Adult.ArchROBO.tmp, outputDirectory = "Save-ProjHeme_adult_granule", load = T,dropCells = T)
#
# Taking Erythroid cells as example
setwd("/home/sunkeyong/MOPA_project/Fetal_adult/Erythroid")
#
addArchRThreads(threads = 13) 
#
ArrowFiles <- c("BoneMarrow_62216.arrow",
                "BoneMarrow_62016.arrow",
                "Spleen_62016.arrow",
                "E8All.arrow",
                "E8Large.arrow",
                "E8Small.arrow",
                "Definitive.erythroid.lineage.arrow",
                "Primitive.erythroid.lineage.arrow")
Fetal.Adult.Erythroid.ArchR <- ArchRProject(ArrowFiles = ArrowFiles,
                          copyArrows = TRUE)
#
x.num <- c(rep(1,800000),rep(0,800000),rep(2,800000))
head(x.num)
x.num.sample <- cbind(sample(x.num,length(Fetal.Adult.Erythroid.ArchR$cellNames)),
                      sample(x.num,length(Fetal.Adult.Erythroid.ArchR$cellNames)),
                      sample(x.num,length(Fetal.Adult.Erythroid.ArchR$cellNames)),
                      sample(x.num,length(Fetal.Adult.Erythroid.ArchR$cellNames)),
                      sample(x.num,length(Fetal.Adult.Erythroid.ArchR$cellNames)))
head(x.num.sample)
colnames(x.num.sample) <- c("F1","F2","F3","F4","F5")
rownames(x.num.sample) <- Fetal.Adult.Erythroid.ArchR$cellNames
#
Fetal.Adult.Erythroid.Seurat <- CreateSeuratObject(counts = t(x.num.sample), assay = "ATAC", project = "Fetal.Adult.Erythroid", min.cells = 0, min.features = 0)
Fetal.Adult.Erythroid.Seurat
Fetal.Adult.Erythroid.Seurat <- NormalizeData(Fetal.Adult.Erythroid.Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
Fetal.Adult.Erythroid.Seurat <- FindVariableFeatures(Fetal.Adult.Erythroid.Seurat, selection.method = "vst", nfeatures = 4)
Fetal.Adult.Erythroid.Seurat <- ScaleData(Fetal.Adult.Erythroid.Seurat, features = rownames(Fetal.Adult.Erythroid.Seurat))
variablegene <- VariableFeatures(object = Fetal.Adult.Erythroid.Seurat)
Fetal.Adult.Erythroid.Seurat <- RunPCA(Fetal.Adult.Erythroid.Seurat, features = variablegene,npcs =2)
DimPlot(Fetal.Adult.Erythroid.Seurat)
#
Fetal.Adult.Erythroid.Seurat$sample <- Fetal.Adult.Erythroid.ArchR$Sample
Fetal.Adult.Erythroid.Seurat@active.ident <- as.factor(Fetal.Adult.Erythroid.Seurat$sample)
table(Fetal.Adult.Erythroid.Seurat@active.ident )
Fetal.Adult.Erythroid.Seurat <- RenameIdents(Fetal.Adult.Erythroid.Seurat, 'BoneMarrow_62016' = 'Adult', 'BoneMarrow_62216' = 'Adult',"Spleen_62016"="Adult")
Fetal.Adult.Erythroid.Seurat <- RenameIdents(Fetal.Adult.Erythroid.Seurat, 'E8All' = 'E8', 'E8Large' = 'E8',"E8Small"="E8")
Fetal.Adult.Erythroid.Seurat <- RenameIdents(Fetal.Adult.Erythroid.Seurat, 'Primitive.erythroid.lineage' = 'our', 'Definitive.erythroid.lineage' = 'our')
#
Fetal.Adult.Erythroid.Seurat$batch <- Fetal.Adult.Erythroid.Seurat@active.ident
table(Fetal.Adult.Erythroid.Seurat$batch)
Fetal.Adult.Erythroid.ArchR$batch <- as.character(Fetal.Adult.Erythroid.Seurat$batch)
#
Fetal.Adult.Erythroid.ArchR <- addIterativeLSI(
  ArchRProj = Fetal.Adult.Erythroid.ArchR,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI_Tile", 
  iterations = 3, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(3,4), 
    sampleCells = 10000, 
    n.start = 10,
    maxClusters = 100
  ), 
  varFeatures = 100000, 
  dimsToUse = 1:100,
  force=T
)
#
Fetal.Adult.Erythroid.ArchR <- addHarmony(
  ArchRProj = Fetal.Adult.Erythroid.ArchR,
  reducedDims = "IterativeLSI_Tile",
  name = "tile_Harmony250",
  groupBy = "batch",
  dimsToUse = 2:50,
  force = T
)
#
Fetal.Adult.Erythroid.ArchR <- addClusters(
  input = Fetal.Adult.Erythroid.ArchR,
  reducedDims = "tile_Harmony250",
  method = "Seurat",
  name = "tile_Harmony250_R5",
  resolution = 5, 
  dimsToUse = 1:50,
  maxClusters= 300,
  force = T,
  sampleCells = NULL
)
Fetal.Adult.Erythroid.ArchR <- addUMAP(
  ArchRProj = Fetal.Adult.Erythroid.ArchR, 
  reducedDims = "tile_Harmony250", 
  name = "tile_Harmony250_UMAP_1", 
  nNeighbors = 100, 
  dimsToUse = 1:50, 
  minDist = 0.4, 
  metric = "cosine",force = T
)
Fetal.Adult.Erythroid.ArchR <- addTSNE(
  ArchRProj = Fetal.Adult.Erythroid.ArchR, 
  reducedDims = "tile_Harmony250", 
  name = "tile_Harmony250_TSNE_1", 
  perplexity = 50, 
  dimsToUse = 1:50,
  maxIterations = 1500,
  force = T
)
#
addArchRThreads(threads = 1) 
####
Fetal.Adult.Erythroid.ArchR <- addGroupCoverages(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "tile_Harmony250_R5")
#
addArchRThreads(threads = 1) 
Fetal.Adult.Erythroid.ArchR <- addReproduciblePeakSet(
  ArchRProj = Fetal.Adult.Erythroid.ArchR, 
  groupBy = "tile_Harmony250_R5",
  pathToMacs2 = "/home/sunkeyong/.local/bin/macs2", force=T
)
#
addArchRThreads(threads = 25)
Fetal.Adult.Erythroid.ArchR <- addPeakMatrix(Fetal.Adult.Erythroid.ArchR,force = T)
#
getAvailableMatrices(Fetal.Adult.Erythroid.ArchR)
Fetal.Adult.Erythroid.Peakset <- getPeakSet(Fetal.Adult.Erythroid.ArchR)
#
Fetal.Adult.Erythroid.ArchR <- addIterativeLSI(
  ArchRProj = Fetal.Adult.Erythroid.ArchR,
  useMatrix = "PeakMatrix", 
  name = "PeakMatrix_LSI", 
  iterations = 3, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(3,4), 
    sampleCells = 10000, 
    n.start = 10,
    maxClusters = 100
  ), 
  varFeatures = 100000, 
  dimsToUse = 1:100,
  force=T
)
#
Fetal.Adult.Erythroid.ArchR <- addHarmony(
  ArchRProj = Fetal.Adult.Erythroid.ArchR,
  reducedDims = "PeakMatrix_LSI",
  name = "PeakMatrix_Harmony250",
  groupBy = "batch",
  dimsToUse = 2:50,
  force = T
)
##
Fetal.Adult.Erythroid.ArchR <- addClusters(
  input = Fetal.Adult.Erythroid.ArchR,
  reducedDims = "PeakMatrix_Harmony250",
  method = "Seurat",
  name = "PeakMatrix_Harmony250_R1",
  resolution = 3, 
  dimsToUse = 1:100,
  maxClusters= 300,
  force = T,
  sampleCells = NULL
)
Fetal.Adult.Erythroid.ArchR <- addClusters(
  input = Fetal.Adult.Erythroid.ArchR,
  reducedDims = "PeakMatrix_Harmony250",
  method = "Seurat",
  name = "PeakMatrix_Harmony250_R3",
  resolution = 3, 
  dimsToUse = 1:100,
  maxClusters= 300,
  force = T,
  sampleCells = NULL
)
#
Fetal.Adult.Erythroid.ArchR <- addUMAP(
  ArchRProj = Fetal.Adult.Erythroid.ArchR, 
  reducedDims = "PeakMatrix_Harmony250", 
  name = "PeakMatrix_Harmony250_UMAP1", 
  nNeighbors = 100, 
  dimsToUse = 1:100, 
  minDist = 0.4, 
  metric = "cosine",force = T
)
##
Fetal.Adult.Erythroid.ArchR <- addUMAP(
  ArchRProj = Fetal.Adult.Erythroid.ArchR, 
  reducedDims = "PeakMatrix_Harmony250", 
  name = "PeakMatrix_Harmony250_UMAP2", 
  nNeighbors = 50, 
  dimsToUse = 1:100, 
  minDist = 0.4, 
  metric = "cosine",force = T
)
##
Fetal.Adult.Erythroid.ArchR <- addTSNE(
  ArchRProj = Fetal.Adult.Erythroid.ArchR, 
  reducedDims = "PeakMatrix_Harmony250", 
  name = "PeakMatrix_Harmony250_TSNE1", 
  perplexity = 50, 
  dimsToUse = 1:100,
  maxIterations = 1500,
  force = T
)
#
##
Fetal.Adult.Erythroid.ArchR <- addTSNE(
  ArchRProj = Fetal.Adult.Erythroid.ArchR, 
  reducedDims = "PeakMatrix_Harmony250", 
  name = "PeakMatrix_Harmony250_TSNE2", 
  perplexity = 100, 
  dimsToUse = 1:100,
  maxIterations = 1500,
  force = T
)
#
#
#
plotEmbedding(ArchRProj = Fetal.Adult.Erythroid.ArchR, colorBy = "cellColData", name = "Sample", embedding = "PeakMatrix_Harmony250_UMAP1")
plotEmbedding(ArchRProj = Fetal.Adult.Erythroid.ArchR, colorBy = "cellColData", name = "batch", embedding = "PeakMatrix_Harmony250_UMAP1")
#
plotEmbedding(ArchRProj = Fetal.Adult.Erythroid.ArchR, colorBy = "cellColData", name = "Sample", embedding = "PeakMatrix_Harmony250_UMAP2")
plotEmbedding(ArchRProj = Fetal.Adult.Erythroid.ArchR, colorBy = "cellColData", name = "batch", embedding = "PeakMatrix_Harmony250_UMAP2")
#
#
plotEmbedding(ArchRProj = Fetal.Adult.Erythroid.ArchR, colorBy = "cellColData", name = "Sample", embedding = "PeakMatrix_Harmony250_TSNE1")
plotEmbedding(ArchRProj = Fetal.Adult.Erythroid.ArchR, colorBy = "cellColData", name = "batch", embedding = "PeakMatrix_Harmony250_TSNE1")
#
plotEmbedding(ArchRProj = Fetal.Adult.Erythroid.ArchR, colorBy = "cellColData", name = "Sample", embedding = "PeakMatrix_Harmony250_TSNE2")
plotEmbedding(ArchRProj = Fetal.Adult.Erythroid.ArchR, colorBy = "cellColData", name = "batch", embedding = "PeakMatrix_Harmony250_TSNE2")
#
#
#
x <- Fetal.Adult.Erythroid.ArchR@embeddings$PeakMatrix_Harmony250_UMAP1
y <- x$df
s.tSNE <- cbind(y$`PeakMatrix_Harmony250#UMAP_Dimension_1`,y$`PeakMatrix_Harmony250#UMAP_Dimension_2`)
#
dim(s.tSNE)
head(s.tSNE)
rownames(s.tSNE) <- rownames(y)
head(s.tSNE)
colnames(s.tSNE) <- c("PC_1","PC_2")
head(s.tSNE)
class(s.tSNE)
#
Fetal.Adult.Erythroid.Seurat@reductions$pca@cell.embeddings <- s.tSNE
table(Fetal.Adult.Erythroid.Seurat@active.ident)
DimPlot(Fetal.Adult.Erythroid.Seurat)
DimPlot(Fetal.Adult.Erythroid.Seurat,split.by = "batch")
#
plot <- DimPlot(subset(Fetal.Adult.Erythroid.Seurat,idents=c("E8")))
plot 
select.cells <- CellSelector(plot = plot)
head(select.cells)
Idents(Fetal.Adult.Erythroid.Seurat, cells = select.cells) <- "E8.def"
DimPlot(Fetal.Adult.Erythroid.Seurat)
table(Fetal.Adult.Erythroid.Seurat@active.ident)
Fetal.Adult.Erythroid.Seurat$batch3 <- Fetal.Adult.Erythroid.Seurat@active.ident
#
#
Fetal.Adult.Erythroid.Seurat@active.ident <- as.factor(Fetal.Adult.Erythroid.Seurat$sample)
table(Fetal.Adult.Erythroid.Seurat@active.ident)
Fetal.Adult.Erythroid.Seurat <- RenameIdents(Fetal.Adult.Erythroid.Seurat, 'BoneMarrow_62016' = 'Adult', 'BoneMarrow_62216' = 'Adult',"Spleen_62016"="Adult")
Fetal.Adult.Erythroid.Seurat <- RenameIdents(Fetal.Adult.Erythroid.Seurat, 'E8All' = 'E8', 'E8Large' = 'E8',"E8Small"="E8")
table(Fetal.Adult.Erythroid.Seurat@active.ident)
Fetal.Adult.Erythroid.Seurat$batch2 <- Fetal.Adult.Erythroid.Seurat@active.ident
#
Fetal.Adult.Erythroid.ArchR$batch2 <- as.character(Fetal.Adult.Erythroid.Seurat$batch2)
Fetal.Adult.Erythroid.ArchR$batch3 <- as.character(Fetal.Adult.Erythroid.Seurat$batch3)
#
#
DimPlot(Fetal.Adult.Erythroid.Seurat,label = T,repel = T,group.by = "batch2")+NoLegend()
DimPlot(Fetal.Adult.Erythroid.Seurat,label = T,repel = T,group.by = "batch3")+NoLegend()
#
Fetal.Adult.Erythroid.Seurat2 <- Fetal.Adult.Erythroid.Seurat
Fetal.Adult.Erythroid.Seurat2@active.ident <- as.factor(Fetal.Adult.Erythroid.Seurat2$batch2)
DimPlot(Fetal.Adult.Erythroid.Seurat2,label = T,repel = T)+NoLegend()
table(Fetal.Adult.Erythroid.Seurat2@active.ident)
Fetal.Adult.Erythroid.Seurat2 <- RenameIdents(Fetal.Adult.Erythroid.Seurat2, 'E8' = 'A.E8', 'Adult' = 'D.Adult',
                      "Definitive.erythroid.lineage"="C.Definitive","Primitive.erythroid.lineage"="B.Primitive")
DimPlot(Fetal.Adult.Erythroid.Seurat2,label = F,repel = T)+NoLegend()
levels(Fetal.Adult.Erythroid.Seurat2) <- c("A.E8","B.Primitive","C.Definitive","D.Adult")
p <- DimPlot(Fetal.Adult.Erythroid.Seurat2,label = F,repel = T,cols =c("#00C185","#99CCCC","#EE82EE","#F8766D"),pt.size = 0.1)+NoLegend()+
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
ggsave("MAEer.celltypes.png",plot=p,device="png",dpi=400,units = "cm",width = 7,height = 7)
###
#
Fetal.Adult.Erythroid.Seurat2$bacth4 <- Fetal.Adult.Erythroid.Seurat2@active.ident
Fetal.Adult.Erythroid.ArchR$batch4 <- as.character(Fetal.Adult.Erythroid.Seurat2$bacth4)
#
#
Fetal.Adult.Erythroid.Seurat$UMAP1 <- s.tSNE[,1]
Fetal.Adult.Erythroid.Seurat$UMAP2 <- s.tSNE[,2]
#
Fetal.Adult.Erythroid.Seurat.metadata <- as.data.frame(Fetal.Adult.Erythroid.Seurat@meta.data)
#
#
#
table(Fetal.Adult.Erythroid.ArchR$batch3)
markersPeaks <- getMarkerFeatures(
  ArchRProj = Fetal.Adult.Erythroid.ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "batch2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList
markerList$Adult
markerList$E8
markerList$Primitive.erythroid.lineage
markerList$Definitive.erythroid.lineage
#
x.p2 <- as.data.frame(cbind(as.character(markerList$Definitive.erythroid.lineage$seqnames),
                            as.character(markerList$Definitive.erythroid.lineage$start),as.character(markerList$Definitive.erythroid.lineage$end)))
write.table(x.p2,file = "rGREAT/Definitive.erythroid.lineage.peak.csv",quote=F, sep = "\t",row.names = F,col.names = F)
x.p2 <- as.data.frame(cbind(as.character(markerList$Primitive.erythroid.lineage$seqnames),
                            as.character(markerList$Primitive.erythroid.lineage$start),as.character(markerList$Primitive.erythroid.lineage$end)))
write.table(x.p2,file = "rGREAT/Primitive.erythroid.lineage.peak.csv",quote=F, sep = "\t",row.names = F,col.names = F)
x.p2 <- as.data.frame(cbind(as.character(markerList$E8$seqnames),as.character(markerList$E8$start),as.character(markerList$E8$end)))
write.table(x.p2,file = "rGREAT/E8.peak.csv",quote=F, sep = "\t",row.names = F,col.names = F)
x.p2 <- as.data.frame(cbind(as.character(markerList$Adult$seqnames),as.character(markerList$Adult$start),as.character(markerList$Adult$end)))
write.table(x.p2,file = "rGREAT/Adult.peak.csv",quote=F, sep = "\t",row.names = F,col.names = F)
#
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
#
Fetal.Adult.Erythroid.ArchR <- addImputeWeights(Fetal.Adult.Erythroid.ArchR, reducedDims = "PeakMatrix_Harmony250")
#
plotEmbedding(ArchRProj = Fetal.Adult.Erythroid.ArchR, colorBy = "GeneScoreMatrix", 
              name = "Hba-a1", embedding = "PeakMatrix_Harmony250_UMAP1",imputeWeights = getImputeWeights(Fetal.Adult.Erythroid.ArchR))
plotEmbedding(ArchRProj = Fetal.Adult.Erythroid.ArchR, colorBy = "GeneScoreMatrix", 
              name = "Myh6", embedding = "PeakMatrix_Harmony250_UMAP1",imputeWeights = getImputeWeights(Fetal.Adult.Erythroid.ArchR))
plotEmbedding(ArchRProj = Fetal.Adult.Erythroid.ArchR, colorBy = "GeneScoreMatrix", 
              name = "Myl7", embedding = "PeakMatrix_Harmony250_UMAP1",imputeWeights = getImputeWeights(Fetal.Adult.Erythroid.ArchR))
plotEmbedding(ArchRProj = Fetal.Adult.Erythroid.ArchR, colorBy = "GeneScoreMatrix", 
              name = "Snx9", embedding = "PeakMatrix_Harmony250_UMAP1",imputeWeights = getImputeWeights(Fetal.Adult.Erythroid.ArchR))
plotEmbedding(ArchRProj = Fetal.Adult.Erythroid.ArchR, colorBy = "GeneScoreMatrix", 
              name = "Myl4", embedding = "PeakMatrix_Harmony250_UMAP1",imputeWeights = getImputeWeights(Fetal.Adult.Erythroid.ArchR))
plotEmbedding(ArchRProj = Fetal.Adult.Erythroid.ArchR, colorBy = "GeneScoreMatrix", 
              name = "Pax6", embedding = "PeakMatrix_Harmony250_UMAP1",imputeWeights = getImputeWeights(Fetal.Adult.Erythroid.ArchR))
plotEmbedding(ArchRProj = Fetal.Adult.Erythroid.ArchR, colorBy = "GeneScoreMatrix", 
              name = "Entpd4b", embedding = "PeakMatrix_Harmony250_UMAP1",imputeWeights = getImputeWeights(Fetal.Adult.Erythroid.ArchR))
#
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Car2"), 
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),
                      upstream = 100000,downstream = 100000,ylim = c(0.01,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Car2)
#
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Slc4a1"), 
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),
                      upstream = 100000,downstream = 100000,ylim = c(0.01,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Slc4a1)
#
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Ldb1"), 
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),
                      upstream = 90000,downstream = 45000,ylim = c(0.01,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Ldb1)
#
#
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Hba-a1"), 
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),
                      upstream = 90000,downstream = 45000,ylim = c(0.01,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$`Hba-a1`)
#
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Hbb-y"), 
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 100,
                      upstream = 200000,downstream = 200000,ylim = c(0.2,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$`Hbb-y`)
#
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Hbb-y"), 
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 100,
                      upstream = 100000,downstream = 100000,ylim = c(0.2,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$`Hbb-y`)
#
#
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Hba-a1"), pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),
                      upstream = 60000,downstream = 30000,ylim = c(0.01,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$`Hba-a1`)
#
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Hbb-y"), pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 100,
                      upstream = 100000,downstream = 100000,ylim = c(0.001,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$`Hbb-y`)
#
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Tfcp2"), 
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 100,
                      upstream = 200000,downstream = 200000,ylim = c(0.2,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Tfcp2)
#
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Klf1"), pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 250,
                      upstream = 50000,downstream = 50000,ylim = c(0.2,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Klf1)
#
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Sox6"), pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 250,
                      upstream = 600000,downstream = 100000,ylim = c(0.2,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Sox6)
#
#
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Myb"), 
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 250,
                      upstream = 100000,downstream = 200000,ylim = c(0.2,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Myb)
#
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Bcl11a"), pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 250,
                      upstream = 300000,downstream = 300000,ylim = c(0.2,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Bcl11a)
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Bcl11b"), pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 250,
                      upstream = 300000,downstream = 300000,ylim = c(0.2,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Bcl11b)
#
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Lmo2"), pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 250,
                      upstream = 80000,downstream = 55000,ylim = c(0.2,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Lmo2)
#
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Tal1"), pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 250,
                      upstream = 100000,downstream = 100000,ylim = c(0.2,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Tal1)  ##
#
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Gata1"), pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 250,
                      upstream = 40000,downstream = 30000,ylim = c(0.2,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Gata1)
#
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Hbb-b2"),pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 250,
                      upstream = 300000,downstream = 300000,ylim = c(0.01,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$`Hbb-b2`)
#
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Fgf18"),pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 250,
                      upstream = 30000,downstream = 50000,ylim = c(0.01,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Fgf18)
#
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Fgf8"),pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 250,
                      upstream = 130000,downstream = 150000,ylim = c(0.01,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Fgf8)
#
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Etv6"),pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 250,
                      upstream = 130000,downstream = 150000,ylim = c(0.01,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Etv6)
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Etv1"),pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 250,
                      upstream = 130000,downstream = 150000,ylim = c(0.01,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Etv1)
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Nr2f2"),pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 250,
                      upstream = 130000,downstream = 150000,ylim = c(0.1,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Nr2f2)
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Cyp46a1"),pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 250,
                      upstream = 130000,downstream = 150000,ylim = c(0.01,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Cyp46a1)
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Sox2"),pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 250,
                      upstream = 130000,downstream = 150000,ylim = c(0.01,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Sox2)
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Sox4"),pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 250,
                      upstream = 130000,downstream = 150000,ylim = c(0.01,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Sox4)
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Nr5a1"),pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 250,
                      upstream = 110000,downstream = 250000,ylim = c(0.01,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Nr5a1)
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("E2f4"),pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 250,
                      upstream = 30000,downstream = 50000,ylim = c(0.01,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$E2f4)
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Egr3"),pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 250,
                      upstream = 50000,downstream = 50000,ylim = c(0.01,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Egr3)
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Egr2"),pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 3),tileSize = 250,
                      upstream = 10000,downstream = 20000,ylim = c(0.01,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$Egr2)
#
markersGS <- getMarkerFeatures(
  ArchRProj = Fetal.Adult.Erythroid.ArchR, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "batch2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.001 & Log2FC >= 1.5")
markerList$Adult$name
markerList$E8$name
markerList$de$name
markerList$pr$name
#
heatmapGS <- markerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.001 & Log2FC >= 1.5", 
  #labelMarkers = markerGenes,
  transpose = TRUE
)
ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
###
###
##
###
table(Fetal.Adult.Erythroid.ArchR$batch4)
Fetal.Adult.Erythroid.ArchR <- addCoAccessibility(
  ArchRProj = Fetal.Adult.Erythroid.ArchR,
  cellsToUse=Cells(subset(Fetal.Adult.Erythroid.Seurat2,idents=c("A.E8"))),
  corCutOff =0.75,
  maxDist = 2e+05,
  dimsToUse = 1:100,
  overlapCutoff = 0.8,
  reducedDims = "PeakMatrix_Harmony250"
)
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Hba-a1"),pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      tileSize = 250,useGroups = c("A.E8"),
                      upstream = 60000,downstream = 30000,ylim = c(0.01,1),log2Norm = T,
                      loops=getCoAccessibility(Fetal.Adult.Erythroid.ArchR))
grid::grid.newpage()
grid::grid.draw(p$`Hba-a1`)
plotPDF(plotList = p, name = "A.E8.hba.Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", ArchRProj = Fetal.Adult.Erythroid.ArchR, addDOC = FALSE, width = 5, height = 5)
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Hbb-y"), pal = c("#00C185","#99CCCC","#EE82EE","#F8766D"),
                      tileSize = 100,useGroups = c("A.E8"),
                      upstream = 100000,downstream = 100000,ylim = c(0.001,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$`Hbb-y`)
plotPDF(plotList = p, name = "A.E8.hbb.Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", ArchRProj = Fetal.Adult.Erythroid.ArchR, addDOC = FALSE, width = 5, height = 5)
##
#
Fetal.Adult.Erythroid.ArchR <- addCoAccessibility(
  ArchRProj = Fetal.Adult.Erythroid.ArchR,
  cellsToUse=Cells(subset(Fetal.Adult.Erythroid.Seurat2,idents=c("B.Primitive"))),
  corCutOff =0.75,
  maxDist = 2e+05,
  dimsToUse = 1:100,
  overlapCutoff = 0.8,
  reducedDims = "PeakMatrix_Harmony250"
)
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Hba-a1"),pal = c("#99CCCC"),
                      tileSize = 250,useGroups = c("B.Primitive"),
                      upstream = 60000,downstream = 30000,ylim = c(0.01,1),log2Norm = T,
                      loops=getCoAccessibility(Fetal.Adult.Erythroid.ArchR))
grid::grid.newpage()
grid::grid.draw(p$`Hba-a1`)
plotPDF(plotList = p, name = "B.Primitive.hba.Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", ArchRProj = Fetal.Adult.Erythroid.ArchR, addDOC = FALSE, width = 5, height = 5)

p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Hbb-y"), pal = c("#99CCCC"),
                      tileSize = 100,useGroups = c("B.Primitive"),
                      upstream = 100000,downstream = 100000,ylim = c(0.001,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$`Hbb-y`)
plotPDF(plotList = p, name = "B.Primitive.hbb.Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", ArchRProj = Fetal.Adult.Erythroid.ArchR, addDOC = FALSE, width = 5, height = 5)
##
#
#
#
#
table(Fetal.Adult.Erythroid.ArchR$batch4)
Fetal.Adult.Erythroid.ArchR <- addCoAccessibility(
  ArchRProj = Fetal.Adult.Erythroid.ArchR,
  cellsToUse=Cells(subset(Fetal.Adult.Erythroid.Seurat2,idents=c("C.Definitive"))),
  corCutOff =0.75,
  maxDist = 2e+05,
  dimsToUse = 1:100,
  overlapCutoff = 0.8,
  reducedDims = "PeakMatrix_Harmony250"
)
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Hba-a1"),pal = c("#EE82EE"),
                      tileSize = 250,useGroups = c("C.Definitive"),
                      upstream = 60000,downstream = 30000,ylim = c(0.01,1),log2Norm = T,
                      loops=getCoAccessibility(Fetal.Adult.Erythroid.ArchR))
grid::grid.newpage()
grid::grid.draw(p$`Hba-a1`)
plotPDF(plotList = p, name = "C.Definitive.hba.Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", ArchRProj = Fetal.Adult.Erythroid.ArchR, addDOC = FALSE, width = 5, height = 5)
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Hbb-y"), pal = c("#EE82EE"),
                      tileSize = 100,useGroups = c("C.Definitive"),
                      upstream = 100000,downstream = 100000,ylim = c(0.001,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$`Hbb-y`)
plotPDF(plotList = p, name = "C.Definitive.hbb.Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", ArchRProj = Fetal.Adult.Erythroid.ArchR, addDOC = FALSE, width = 5, height = 5)
#
table(Fetal.Adult.Erythroid.ArchR$batch4)
Fetal.Adult.Erythroid.ArchR <- addCoAccessibility(
  ArchRProj = Fetal.Adult.Erythroid.ArchR,
  cellsToUse=Cells(subset(Fetal.Adult.Erythroid.Seurat2,idents=c("D.Adult"))),
  corCutOff =0.75,
  maxDist = 2e+05,
  dimsToUse = 1:100,
  overlapCutoff = 0.8,
  reducedDims = "PeakMatrix_Harmony250"
)
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Hba-a1"),pal = c("#F8766D"),
                      tileSize = 250,useGroups = c("D.Adult"),
                      upstream = 60000,downstream = 30000,ylim = c(0.01,1),log2Norm = T,
                      loops=getCoAccessibility(Fetal.Adult.Erythroid.ArchR))
grid::grid.newpage()
grid::grid.draw(p$`Hba-a1`)
plotPDF(plotList = p, name = "D.Adult.hba.Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", ArchRProj = Fetal.Adult.Erythroid.ArchR, addDOC = FALSE, width = 5, height = 5)
p <- plotBrowserTrack(ArchRProj = Fetal.Adult.Erythroid.ArchR, groupBy = "batch4", geneSymbol = c("Hbb-y"), pal = c("#F8766D"),
                      tileSize = 100,useGroups = c("D.Adult"),
                      upstream = 100000,downstream = 100000,ylim = c(0.001,1),log2Norm = T)
grid::grid.newpage()
grid::grid.draw(p$`Hbb-y`)
plotPDF(plotList = p, name = "D.Adult.hbb.Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", ArchRProj = Fetal.Adult.Erythroid.ArchR, addDOC = FALSE, width = 5, height = 5)
##
##
#
# 28-May-2022 09:33
load("pro1.RData")
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList
markerList$Adult
markerList$E8
markerList$Primitive.erythroid.lineage
markerList$Definitive.erythroid.lineage
###
#
MAEendoPeakset.metadata <- as.data.frame(cbind(as.character(MAEendoPeakset@seqnames),MAEendoPeakset@ranges@start,MAEendoPeakset@ranges@start+500,MAEendoPeakset$peakType))
head(MAEendoPeakset.metadata)
colnames(MAEendoPeakset.metadata) <- c("chr","start","end","peaktype")
MAEendoPeakset.metadata
rownames(MAEendoPeakset.metadata) <- paste(MAEendoPeakset.metadata$chr,MAEendoPeakset.metadata$start,MAEendoPeakset.metadata$end,sep = "-")
head(MAEendoPeakset.metadata)
##
Pri.peak <- MAEendoPeakset.metadata[paste(markerList$Primitive.erythroid.lineage$seqnames,markerList$Primitive.erythroid.lineage$start,markerList$Primitive.erythroid.lineage$end,sep="-"),]
Pri.peak <- subset(Pri.peak, peaktype == c("Distal","Intronic"))
table(Pri.peak$peaktype)
#
Def.peak <- MAEendoPeakset.metadata[paste(markerList$Definitive.erythroid.lineage$seqnames,markerList$Definitive.erythroid.lineage$start,markerList$Definitive.erythroid.lineage$end,sep="-"),]
Def.peak <- subset(Def.peak, peaktype == c("Distal","Intronic"))
table(Def.peak$peaktype)
#
Adult.peak <- MAEendoPeakset.metadata[paste(markerList$Adult$seqnames,markerList$Adult$start,markerList$Adult$end,sep="-"),]
Adult.peak <- subset(Adult.peak, peaktype == c("Distal","Intronic"))
table(Adult.peak$peaktype)
#
E8.peak <- MAEendoPeakset.metadata[paste(markerList$E8$seqnames,markerList$E8$start,markerList$E8$end,sep="-"),]
E8.peak <- subset(E8.peak, peaktype == c("Distal","Intronic"))
table(E8.peak$peaktype)
#
Pri.peak.F <- as.data.frame(cbind(as.character(Pri.peak$chr),as.numeric(Pri.peak$start),as.numeric(Pri.peak$start)+1))
Def.peak.F <- as.data.frame(cbind(as.character(Def.peak$chr),as.numeric(Def.peak$start),as.numeric(Def.peak$start)+1))
Adult.peak.F <- as.data.frame(cbind(as.character(Adult.peak$chr),as.numeric(Adult.peak$start),as.numeric(Adult.peak$start)+1))
E8.peak.F <- as.data.frame(cbind(as.character(E8.peak$chr),as.numeric(E8.peak$start),as.numeric(E8.peak$start)+1))
#
write.table(Pri.peak.F,file = "conservation/Pri.peak.F.csv",sep = "\t",quote = F,col.names = F,row.names =F)
write.table(Def.peak.F,file = "conservation/Def.peak.F.csv",sep = "\t",quote = F,col.names = F,row.names =F)
write.table(Adult.peak.F,file = "conservation/Adult.peak.F.csv",sep = "\t",quote = F,col.names = F,row.names =F)
write.table(E8.peak.F,file = "conservation/E8.peak.F.csv",sep = "\t",quote = F,col.names = F,row.names =F)
#
Pri.peak.F.phastCons <- read.table("conservation/Pri.Bin500.phastCons.center.2")
Pri.peak.F.phyloP60way <- read.table("conservation/Pri.Bin500.phyloP60way.center.2")
Def.peak.F.phastCons <- read.table("conservation/Def.Bin500.phastCons.center.2")
Def.peak.F.phyloP60way <- read.table("conservation/Def.Bin500.phyloP60way.center.2")
Adult.peak.F.phastCons <- read.table("conservation/Adult.Bin500.phastCons.center.2")
Adult.peak.F.phyloP60way <- read.table("conservation/Adult.Bin500.phyloP60way.center.2")
E8.peak.F.phastCons <- read.table("conservation/E8.Bin500.phastCons.center.2")
E8.peak.F.phyloP60way <- read.table("conservation/E8.Bin500.phyloP60way.center.2")
#
boxplot(E8.peak.F.phastCons$V7,Pri.peak.F.phastCons$V7,Def.peak.F.phastCons$V7,Adult.peak.F.phastCons$V7)
boxplot(E8.peak.F.phyloP60way$V7,Pri.peak.F.phyloP60way$V7,Def.peak.F.phyloP60way$V7,Adult.peak.F.phyloP60way$V7)
#
boxplot(E8.peak.F.phastCons$V7,c(Pri.peak.F.phastCons$V7,Def.peak.F.phastCons$V7),Adult.peak.F.phastCons$V7)
boxplot(E8.peak.F.phyloP60way$V7,c(Pri.peak.F.phyloP60way$V7,Def.peak.F.phyloP60way$V7),Adult.peak.F.phyloP60way$V7)
#
boxplot(c(E8.peak.F.phastCons$V7,c(Pri.peak.F.phastCons$V7,Def.peak.F.phastCons$V7)),Adult.peak.F.phastCons$V7)
boxplot(c(E8.peak.F.phyloP60way$V7,c(Pri.peak.F.phyloP60way$V7,Def.peak.F.phyloP60way$V7)),Adult.peak.F.phyloP60way$V7)
#
#######
# loading scRNA-seq data
###
####
setwd("/home/sunkeyong/MAE/MAE_Erythroid/")
#
NCB.count <- read.table("scRNA_download/NCB/normalisedCounts.tsv")
head(NCB.count)
dim(NCB.count)
NCB.count[1:2,1:2]
colnames(NCB.count)
#
genenames <- as.data.frame(mapIds(org.Mm.eg.db,keys = rownames(NCB.count), column = 'SYMBOL',keytype = 'ENSEMBL'))
head(genenames)
colnames(genenames) <- "geneid"
head(genenames)
genenames2 <- genenames %>% drop_na(geneid)
head(genenames2)
genenames2$a <- "a"
genenames3 <- genenames2[!duplicated(genenames2[,"geneid"]),]
#
NCB.count <- NCB.count[rownames(genenames3),]
NCB.count[1:3,1:3]
NCB.count <- NCB.count
rownames(NCB.count) <- genenames3$geneid
#
NCB <- CreateSeuratObject(NCB.count,project="E8",min.cells = 0,min.features = 0)
NCB
VlnPlot(NCB, features = c("nFeature_RNA", "nCount_RNA"),pt.size = 0)
#NCB <- subset(NCB, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
NCB <- NormalizeData(NCB, normalization.method = "LogNormalize", scale.factor = 10000)
NCB <- FindVariableFeatures(NCB, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(NCB)
NCB <- ScaleData(NCB, features = all.genes)
#Perform linear dimensional reduction
NCB <- RunPCA(NCB, features = VariableFeatures(object = NCB))
ElbowPlot(NCB)
NCB <- FindNeighbors(NCB, dims = 1:20)
NCB <- FindClusters(NCB, resolution = 0.5)
NCB <- RunUMAP(NCB, dims = 1:20)
NCB <- RunTSNE(NCB, dims = 1:20)
#
DimPlot(NCB,repel = T, label = T)
DimPlot(NCB,repel = T, label = T,reduction = "umap")
#
FeaturePlot(NCB,features = c("Hbb-y"))
FeaturePlot(NCB,features = c("Hbb-bs"))
FeaturePlot(NCB,features = c("Hbb-bt"))
FeaturePlot(NCB,features = c("Hbb-bh1"),min.cutoff = 0)
FeaturePlot(NCB,features = c("Hbb-bh2"))
#
NCB.er <- subset(NCB,idents=c("0","18"))
#
load("MOCA/Fetal.Adult.Erythroid.Seurat.RData")
Fetal.Adult.Erythroid.Seurat@active.ident <- as.factor(Fetal.Adult.Erythroid.Seurat$Main_cell_type)
table(Fetal.Adult.Erythroid.Seurat@active.ident)
sci.er <- subset(Fetal.Adult.Erythroid.Seurat,idents = c("Primitive erythroid lineage","Definitive erythroid lineage"))
table(sci.er@active.ident)
sci.er@active.ident <- as.factor(sci.er$Main_cell_type.Stage)
sci.er <- subset(sci.er,downsample= 500)
table(sci.er@active.ident)
sci.er@active.ident <- as.factor(sci.er$Main_cell_type)
#
#
load("scRNA_download/DroSQ.er.RData")
load("scRNA_download/GGJs.er.RData")
#
DroSQ.er$stage <- "Adult.SQ"
GGJs.er$stage <- "Adult.GGJ"
sci.er$stage <- "E10.5.E13.5.SD"
NCB.er$stage <- "E8.25"
#
Fetal.Adult.Erythroid.Seurat <- merge(DroSQ.er,y=c(GGJs.er,sci.er,NCB.er))
Fetal.Adult.Erythroid.Seurat
Fetal.Adult.Erythroid.Seurat@active.ident <- as.factor(Fetal.Adult.Erythroid.Seurat$stage)
#
VlnPlot(Fetal.Adult.Erythroid.Seurat, features = c("nFeature_RNA", "nCount_RNA"),pt.size = 0)
#Fetal.Adult.Erythroid.Seurat <- subset(Fetal.Adult.Erythroid.Seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
Fetal.Adult.Erythroid.Seurat <- NormalizeData(Fetal.Adult.Erythroid.Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
Fetal.Adult.Erythroid.Seurat <- FindVariableFeatures(Fetal.Adult.Erythroid.Seurat, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Fetal.Adult.Erythroid.Seurat)
Fetal.Adult.Erythroid.Seurat <- ScaleData(Fetal.Adult.Erythroid.Seurat, features = all.genes)
#Perform linear dimensional reduction
Fetal.Adult.Erythroid.Seurat <- RunPCA(Fetal.Adult.Erythroid.Seurat, features = VariableFeatures(object = Fetal.Adult.Erythroid.Seurat))
ElbowPlot(Fetal.Adult.Erythroid.Seurat,ndims = 50)
Fetal.Adult.Erythroid.Seurat <- FindNeighbors(Fetal.Adult.Erythroid.Seurat, dims = 1:20)
#Fetal.Adult.Erythroid.Seurat <- FindClusters(Fetal.Adult.Erythroid.Seurat, resolution = 0.5)
Fetal.Adult.Erythroid.Seurat <- RunUMAP(Fetal.Adult.Erythroid.Seurat, dims = 1:20)
Fetal.Adult.Erythroid.Seurat <- RunTSNE(Fetal.Adult.Erythroid.Seurat, dims = 1:20)
#
DimPlot(Fetal.Adult.Erythroid.Seurat,repel=T,label=T)+NoLegend()
#
Idents(Fetal.Adult.Erythroid.Seurat,cells=Cells(subset(sci.er,idents="Definitive erythroid lineage"))) <- "SD.Def"
Idents(Fetal.Adult.Erythroid.Seurat,cells=Cells(subset(sci.er,idents="Primitive erythroid lineage"))) <- "SD.Pri"
#
VlnPlot(Fetal.Adult.Erythroid.Seurat,features = c("Hbb-y","Hbb-bh1","Hbb-bh2","Hbb-bt","Hbb-bs"),pt.size = 0)+NoLegend()
VlnPlot(Fetal.Adult.Erythroid.Seurat,features = c("Hba-x","Hba-a1","Hba-a2"),pt.size = 0)+NoLegend()
###
levels(Fetal.Adult.Erythroid.Seurat) <- c("E8.25","SD.Pri","SD.Def","Adult.GGJ","Adult.SQ")
VlnPlot(Fetal.Adult.Erythroid.Seurat,features = c("Hbb-y","Hbb-bh1","Hbb-bh2","Hbb-bt","Hbb-bs"),pt.size = 0)+NoLegend()
VlnPlot(Fetal.Adult.Erythroid.Seurat,features = c("Hbb-b1","Hbb-b2"),pt.size = 0)+NoLegend()
VlnPlot(Fetal.Adult.Erythroid.Seurat,features = c("Hba-x","Hba-a1","Hba-a2"),pt.size = 0)+NoLegend()
VlnPlot(Fetal.Adult.Erythroid.Seurat,features = c("Hbb-y","Hbb-bh1","Hbb-bh2","Hbb-bt","Hbb-bs","Hbb-b1","Hbb-b2","Hba-x","Hba-a1","Hba-a2"),
        pt.size = 0,cols = c("#00C185","#99CCCC","#EE82EE","#AE00C1","#C10061"))+NoLegend()
###
StackedVlnPlot(Fetal.Adult.Erythroid.Seurat,features = c("Hbb-y","Hbb-bh1","Hbb-bh2","Hbb-bt","Hbb-bs","Hbb-b1","Hbb-b2","Hba-x","Hba-a1","Hba-a2"),
               pt.size = 0,cols = c("#00C185","#99CCCC","#EE82EE","#AE00C1","#C10061"))+NoLegend()
###
###
load("MOCA/Fetal.Adult.Erythroid.Seurat.Main_cell_type.Stage.pseudo.RData")
VlnPlot(Fetal.Adult.Erythroid.Seurat.Main_cell_type.Stage.pseudo.seurat,
        features = c("Hbb-y","Hbb-bh1","Hbb-bh2","Hbb-bt","Hbb-bs"),pt.size = 0)+NoLegend()
Fetal.Adult.Erythroid.Seurat.Main_cell_type.Stage.pseudo.seurat.er <- subset(Fetal.Adult.Erythroid.Seurat.Main_cell_type.Stage.pseudo.seurat,idents=(c("Primitive erythroid lineage","Definitive erythroid lineage")))
DimPlot(Fetal.Adult.Erythroid.Seurat.Main_cell_type.Stage.pseudo.seurat.er)
#
Fetal.Adult.Erythroid.Seurat.Main_cell_type.Stage.pseudo.seurat.er$stage <- "E10.5.E13.5.SD"
Fetal.Adult.Erythroid.Seurat2 <- merge(GGJs.er,y=c(Fetal.Adult.Erythroid.Seurat.Main_cell_type.Stage.pseudo.seurat.er,NCB.er))
Fetal.Adult.Erythroid.Seurat2
Fetal.Adult.Erythroid.Seurat2@active.ident <- as.factor(Fetal.Adult.Erythroid.Seurat2$stage)
#
VlnPlot(Fetal.Adult.Erythroid.Seurat2, features = c("nFeature_RNA", "nCount_RNA"),pt.size = 0)
#Fetal.Adult.Erythroid.Seurat2 <- subset(Fetal.Adult.Erythroid.Seurat2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
Fetal.Adult.Erythroid.Seurat2 <- NormalizeData(Fetal.Adult.Erythroid.Seurat2, normalization.method = "LogNormalize", scale.factor = 10000)
Fetal.Adult.Erythroid.Seurat2 <- FindVariableFeatures(Fetal.Adult.Erythroid.Seurat2, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(Fetal.Adult.Erythroid.Seurat2)
Fetal.Adult.Erythroid.Seurat2 <- ScaleData(Fetal.Adult.Erythroid.Seurat2, features = all.genes)
#Perform linear dimensional reduction
Fetal.Adult.Erythroid.Seurat2 <- RunPCA(Fetal.Adult.Erythroid.Seurat2, features = VariableFeatures(object = Fetal.Adult.Erythroid.Seurat2))
ElbowPlot(Fetal.Adult.Erythroid.Seurat2,ndims = 50)
Fetal.Adult.Erythroid.Seurat2 <- FindNeighbors(Fetal.Adult.Erythroid.Seurat2, dims = 1:20)
#Fetal.Adult.Erythroid.Seurat2 <- FindClusters(Fetal.Adult.Erythroid.Seurat2, resolution = 0.5)
Fetal.Adult.Erythroid.Seurat2 <- RunUMAP(Fetal.Adult.Erythroid.Seurat2, dims = 1:20)
Fetal.Adult.Erythroid.Seurat2 <- RunTSNE(Fetal.Adult.Erythroid.Seurat2, dims = 1:20)
#
DimPlot(Fetal.Adult.Erythroid.Seurat2,repel=T,label=T)+NoLegend()
#
Idents(Fetal.Adult.Erythroid.Seurat2,cells=Cells(subset(Fetal.Adult.Erythroid.Seurat.Main_cell_type.Stage.pseudo.seurat,idents="Definitive erythroid lineage"))) <- "SD.Def"
Idents(Fetal.Adult.Erythroid.Seurat2,cells=Cells(subset(Fetal.Adult.Erythroid.Seurat.Main_cell_type.Stage.pseudo.seurat,idents="Primitive erythroid lineage"))) <- "SD.Pri"
#
VlnPlot(Fetal.Adult.Erythroid.Seurat2,features = c("Hbb-y","Hbb-bh1","Hbb-bh2","Hbb-bt","Hbb-bs"),pt.size = 0)+NoLegend()
VlnPlot(Fetal.Adult.Erythroid.Seurat2,features = c("Hba-x","Hba-a1","Hba-a2"),pt.size = 0)+NoLegend()
###
save.image("Mouse.Fetal.Adult_Erythroid.RData")
#