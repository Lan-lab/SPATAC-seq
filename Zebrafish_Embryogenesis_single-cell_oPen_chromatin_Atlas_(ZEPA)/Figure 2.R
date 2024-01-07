# Figure 2
# cell type annotation 
# taking 10 hpf as an example 
setwd("/home/sunkeyong/ZEPA/OBO/hpf10")
#
###Zebrafish_
library(ArchR)
library(Seurat)
library("BSgenome.Drerio.UCSC.danRer11")
library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(pheatmap)
library(ComplexHeatmap)
library("RColorBrewer")
library("circlize")
#
# loading annotation files
load("/home/sunkeyong/ZEPA/GRCz11_anno/Lawson/storage/GRCz11.geneAnnotation.RData")
load("/home/sunkeyong/ZEPA/GRCz11_anno/Lawson/storage/genomeAnnotation.RData")
chr_rm <- read.table("/home/sunkeyong/ZEPA/GRCz11_anno/GRCz11_toplevel.chromosome.V2.txt",text = )
##
addArchRThreads(threads = 12) 
getArchRChrPrefix()
# load bed files
BedFiles <- c("/home/sunkeyong/ZEPA/All_arrow/10hpf.SPATACseq.batch1.bed.gz")
names(BedFiles) = c("10hpf.SPATACseq.batch1")
# Create ArchR arrow file for downstream analysis
ArrowFiles <- createArrowFiles(
  inputFiles = BedFiles,
  sampleNames = names(BedFiles),
  minTSS = 4,
  minFrags = 500,
  offsetPlus = 0,
  offsetMinus = 0,
  geneAnnotation = GRCz11.geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  addTileMat = T,
  addGeneScoreMat = T,
  excludeChr =c("MT",paste("chr",chr_rm$V1,sep = "")))
# calculate Doublet score 
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)
#
zhpf10 <- ArchRProject(ArrowFiles = ArrowFiles,
                       geneAnnotation = GRCz11.geneAnnotation,
                       genomeAnnotation = genomeAnnotation,
                       copyArrows = T)
#
zhpf10
#
df <- getCellColData(zhpf10, select = c("log10(nFrags)", "TSSEnrichment"))
df
p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 5, lty = "dashed") + geom_vline(xintercept = log10(1000), lty = "dashed")
p
#
idxPass <- which(zhpf10$TSSEnrichment >= 5 & zhpf10$nFrags >= 1000)
length(idxPass)
cellsPass <- zhpf10$cellNames[idxPass]
zhpf10 <- zhpf10[cellsPass, ]
#
#
#
df <- getCellColData(zhpf10, select = c("log10(nFrags)", "TSSEnrichment"))
df
p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = log10(1000), lty = "dashed")
p
#
zhpf10 <- addIterativeLSI(ArchRProj = zhpf10,
                          useMatrix = "TileMatrix", 
                          name = "IterativeLSI_Tile_LM1", 
                          iterations = 3, 
                          clusterParams = list( resolution = c(2,3), sampleCells = 10000, n.start = 10,maxClusters = 45), 
                          corCutOff = 0.5,
                          varFeatures = 50000, 
                          dimsToUse = 1:100,
                          LSIMethod = 1,
                          force=T)
zhpf10 <- addIterativeLSI(ArchRProj = zhpf10,
                          useMatrix = "TileMatrix", 
                          name = "IterativeLSI_Tile_LM1x", 
                          iterations = 3, 
                          clusterParams = list( resolution = c(2,3), sampleCells = 10000, n.start = 10,maxClusters = 45), 
                          corCutOff = 0.75,
                          varFeatures = 50000, 
                          dimsToUse = 1:100,
                          LSIMethod = 1,
                          force=T)
save.image("zhpf10.RData")
barplot(zhpf10@reducedDims$IterativeLSI_Tile_LM1$corToDepth$scaled)
barplot(zhpf10@reducedDims$IterativeLSI_Tile_LM1x$corToDepth$scaled)
#
zhpf10 <- addTSNE(ArchRProj = zhpf10, reducedDims = "IterativeLSI_Tile_LM1", name = "Tile_LM1_TSNE1", 
                  perplexity = 100, dimsToUse = 1:100,corCutOff = 0.4,maxIterations = 1000,force = T)
zhpf10 <- addTSNE(ArchRProj = zhpf10, reducedDims = "IterativeLSI_Tile_LM1", name = "Tile_LM1_TSNE2", 
                  perplexity = 100, dimsToUse = 1:100,corCutOff = 0.5,maxIterations = 1000,force = T)
zhpf10 <- addTSNE(ArchRProj = zhpf10, reducedDims = "IterativeLSI_Tile_LM1", name = "Tile_LM1_TSNE3", 
                  perplexity = 100, dimsToUse = 1:100,corCutOff = 0.8,maxIterations = 1000,force = T)
save.image("zhpf10.RData")
#
zhpf10 <- addUMAP(ArchRProj = zhpf10, reducedDims = "IterativeLSI_Tile_LM1", name = "Tile_LM1_UMAP1", 
                  nNeighbors = 50, dimsToUse = 1:100, minDist = 0.3, metric = "cosine",force = T,corCutOff = 0.4)
zhpf10 <- addUMAP(ArchRProj = zhpf10, reducedDims = "IterativeLSI_Tile_LM1", name = "Tile_LM1_UMAP2", 
                  nNeighbors = 50, dimsToUse = 1:100, minDist = 0.3, metric = "cosine",force = T,corCutOff = 0.5)
zhpf10 <- addUMAP(ArchRProj = zhpf10, reducedDims = "IterativeLSI_Tile_LM1", name = "Tile_LM1_UMAP3", 
                  nNeighbors = 50, dimsToUse = 1:100, minDist = 0.3, metric = "cosine",force = T,corCutOff = 0.8)
#
save.image("zhpf10.RData")
zhpf10 <- addTSNE(ArchRProj = zhpf10, reducedDims = "IterativeLSI_Tile_LM1x", name = "Tile_LM1_TSNE1x", 
                  perplexity = 100, dimsToUse = 1:100,corCutOff = 0.4,maxIterations = 1000,force = T)
zhpf10 <- addTSNE(ArchRProj = zhpf10, reducedDims = "IterativeLSI_Tile_LM1x", name = "Tile_LM1_TSNE2x", 
                  perplexity = 100, dimsToUse = 1:100,corCutOff = 0.5,maxIterations = 1000,force = T)
zhpf10 <- addTSNE(ArchRProj = zhpf10, reducedDims = "IterativeLSI_Tile_LM1x", name = "Tile_LM1_TSNE3x", 
                  perplexity = 100, dimsToUse = 1:100,corCutOff = 0.8,maxIterations = 1000,force = T)
save.image("zhpf10.RData")
#
zhpf10 <- addUMAP(ArchRProj = zhpf10, reducedDims = "IterativeLSI_Tile_LM1x", name = "Tile_LM1_UMAP1x", 
                  nNeighbors = 50, dimsToUse = 1:100, minDist = 0.3, metric = "cosine",force = T,corCutOff = 0.4)
zhpf10 <- addUMAP(ArchRProj = zhpf10, reducedDims = "IterativeLSI_Tile_LM1x", name = "Tile_LM1_UMAP2x", 
                  nNeighbors = 50, dimsToUse = 1:100, minDist = 0.3, metric = "cosine",force = T,corCutOff = 0.5)
zhpf10 <- addUMAP(ArchRProj = zhpf10, reducedDims = "IterativeLSI_Tile_LM1x", name = "Tile_LM1_UMAP3x", 
                  nNeighbors = 50, dimsToUse = 1:100, minDist = 0.3, metric = "cosine",force = T,corCutOff = 0.8)
#
save.image("zhpf10.RData")
#
plotEmbedding(ArchRProj = zhpf10, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_TSNE1")
plotEmbedding(ArchRProj = zhpf10, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_TSNE2")
plotEmbedding(ArchRProj = zhpf10, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_TSNE3")
plotEmbedding(ArchRProj = zhpf10, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_UMAP1")
plotEmbedding(ArchRProj = zhpf10, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_UMAP2")
plotEmbedding(ArchRProj = zhpf10, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_UMAP3")
plotEmbedding(ArchRProj = zhpf10, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_TSNE1x")
plotEmbedding(ArchRProj = zhpf10, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_TSNE2x")
plotEmbedding(ArchRProj = zhpf10, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_TSNE3x")
plotEmbedding(ArchRProj = zhpf10, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_UMAP1x")
plotEmbedding(ArchRProj = zhpf10, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_UMAP2x")
plotEmbedding(ArchRProj = zhpf10, colorBy = "cellColData", name = "Sample", embedding = "Tile_LM1_UMAP3x")
#
#
#
zhpf10 <- addClusters(input = zhpf10,name = "Tile_LM1.C3",reducedDims = "IterativeLSI_Tile_LM1",corCutOff = 0.5,
                      resolution = 3,method = "Seurat", dimsToUse = 1:100,maxClusters= 300,force = T)
zhpf10 <- addClusters(input = zhpf10,name = "Tile_LM1.C3x",reducedDims = "IterativeLSI_Tile_LM1x",corCutOff = 0.5,
                      resolution = 3,method = "Seurat", dimsToUse = 1:100,maxClusters= 300,force = T)
zhpf10 <- addClusters(input = zhpf10,name = "Tile_LM1.C5",reducedDims = "IterativeLSI_Tile_LM1",corCutOff = 0.5,
                      resolution = 5,method = "Seurat", dimsToUse = 1:100,maxClusters= 300,force = T)
zhpf10 <- addClusters(input = zhpf10,name = "Tile_LM1.C5x",reducedDims = "IterativeLSI_Tile_LM1x",corCutOff = 0.5,
                      resolution = 5,method = "Seurat", dimsToUse = 1:100,maxClusters= 300,force = T)
#
#
plotEmbedding(ArchRProj = zhpf10, colorBy = "cellColData", name = "Tile_LM1.C3", embedding = "Tile_LM1_TSNE2")
plotEmbedding(ArchRProj = zhpf10, colorBy = "cellColData", name = "Tile_LM1.C3x", embedding = "Tile_LM1_TSNE2x")
plotEmbedding(ArchRProj = zhpf10, colorBy = "cellColData", name = "Tile_LM1.C5", embedding = "Tile_LM1_TSNE2")
plotEmbedding(ArchRProj = zhpf10, colorBy = "cellColData", name = "Tile_LM1.C5x", embedding = "Tile_LM1_TSNE2x")
#
#
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf10.F.RData")
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf10.W.RData")
hpf10.F <- NormalizeData(hpf10.F, normalization.method = "LogNormalize", scale.factor = 10000)
hpf10.F <- FindVariableFeatures(hpf10.F, selection.method = "vst", nfeatures = 3000)
hpf10.F <- ScaleData(hpf10.F, features = rownames(hpf10.F))
hpf10.F <- RunPCA(hpf10.F, features = VariableFeatures(object = hpf10.F))
hpf10.F <- FindNeighbors(hpf10.F, dims = 1:50)
hpf10.F <- RunUMAP(hpf10.F, dims = 1:50)
hpf10.F <- RunTSNE(hpf10.F, dims = 1:50)
DimPlot(hpf10.F, reduction = "umap",group.by = "cell_type",label=T,repel = T)+NoLegend()
#
hpf10.W <- NormalizeData(hpf10.W, normalization.method = "LogNormalize", scale.factor = 10000)
hpf10.W <- FindVariableFeatures(hpf10.W, selection.method = "vst", nfeatures = 3000)
hpf10.W <- ScaleData(hpf10.W, features = rownames(hpf10.W))
hpf10.W <- RunPCA(hpf10.W, features = VariableFeatures(object = hpf10.W))
hpf10.W <- FindNeighbors(hpf10.W, dims = 1:50)
hpf10.W <- RunUMAP(hpf10.W, dims = 1:50)
hpf10.W <- RunTSNE(hpf10.W, dims = 1:50)
DimPlot(hpf10.W, reduction = "umap",group.by = "cell_type",label=T,repel = T)+NoLegend()
#
DimPlot(hpf10.F, reduction = "tsne",group.by = "cell_type",label=T,repel = T)+NoLegend()
DimPlot(hpf10.W, reduction = "tsne",group.by = "cell_type",label=T,repel = T)+NoLegend()
#
#
#
library(dplyr)
hpf10.W@active.ident <- as.factor(hpf10.W$cell_type)
hpf10.W.markers <- FindAllMarkers(hpf10.W, only.pos = TRUE, min.pct = 0.25)
hpf10.W.markers.top10 <- hpf10.W.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(hpf10.W, features = unique(hpf10.W.markers.top10$gene)) + NoLegend()
View(hpf10.W.markers)
#
hpf10.F@active.ident <- as.factor(hpf10.F$cell_type)
hpf10.F.markers <- FindAllMarkers(hpf10.F, only.pos = TRUE, min.pct = 0.25)
hpf10.F.markers.top10 <- hpf10.F.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(hpf10.F, features = unique(hpf10.F.markers.top10$gene)) + NoLegend()
View(hpf10.F.markers)
###
###
FeaturePlot(hpf10.F,features = "sox3")
FeaturePlot(hpf10.W,features = "sox3")
#
VlnPlot(hpf10.W,features = "zic3",group.by = "cell_type")+NoLegend()
VlnPlot(hpf10.W,features = "sox10",group.by = "cell_type")+NoLegend()
VlnPlot(hpf10.W,features = "pax2a",group.by = "cell_type")+NoLegend()
VlnPlot(hpf10.W,features = "aplnr2",group.by = "cell_type")+NoLegend()
VlnPlot(hpf10.F,features = "aplnr2",group.by = "cell_type")+NoLegend()
#
library(Seurat)
library(cowplot)
library(harmony)
hpf10 <- merge(hpf10.F,y=hpf10.W)
hpf10 <- 
  hpf10 <- hpf10 %>% Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = zhpf10.seurat@var.genes, npcs = 20, verbose = FALSE)
#
p1 <- DimPlot(object = hpf10, reduction = "pca", pt.size = .1, group.by = "Sample")
p2 <- VlnPlot(object = hpf10, features = "PC_1", group.by = "Sample", pt.size = .1)
plot_grid(p1,p2)
#
hpf10 <- hpf10 %>% 
  RunHarmony("Sample", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(hpf10, 'harmony')
harmony_embeddings[1:5, 1:5]
#
p1 <- DimPlot(object = hpf10, reduction = "harmony", pt.size = .1, group.by = "Sample")
p2 <- VlnPlot(object = hpf10, features = "harmony_1", group.by = "Sample", pt.size = .1)
plot_grid(p1,p2)
#
hpf10 <- hpf10 %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  RunTSNE(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20)
#
DimPlot(hpf10, reduction = "umap", group.by = "cell_type", pt.size = .1, split.by = 'Sample')
DimPlot(hpf10, reduction = "tsne", group.by = "cell_type", pt.size = .1, split.by = 'Sample')
#
DimPlot(hpf10, reduction = "umap", group.by = "cell_type",label = T,repel = T, pt.size = .1)+NoLegend()
DimPlot(hpf10, reduction = "tsne", group.by = "cell_type",label = T,repel = T, pt.size = .1)+NoLegend()
#
hpf10
#
addArchRThreads(threads = 1) 
hpf10
zhpf10
zhpf10 <- addGeneIntegrationMatrix(
  ArchRProj = zhpf10, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GIM_Tile_LSI.M1.1",
  reducedDims = "IterativeLSI_Tile_LM1",
  seRNA = hpf10,
  addToArrow = T,
  force = TRUE,
  dimsToUse = 2:100,
  #corCutOff = 0.75,
  UMAPParams = list(n_neighbors = 50, min_dist = 0.3, metric = "cosine", verbose = FALSE),
  nGenes = 3000,
  sampleCellsATAC = 17723,
  sampleCellsRNA = 11394,
  groupRNA = "cell_type",
  nameCell = "predictedCell_1",
  nameGroup = "predictedGroup_1",
  nameScore = "predictedScore_1")
#
#
addArchRThreads(threads = 1) 
hpf10.W
zhpf10
zhpf10 <- addGeneIntegrationMatrix(
  ArchRProj = zhpf10, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GIM_Tile_LSI.M1.2",
  reducedDims = "IterativeLSI_Tile_LM1",
  seRNA = hpf10.W,
  addToArrow = T,
  force = TRUE,
  dimsToUse = 2:100,
  #corCutOff = 0.75,
  UMAPParams = list(n_neighbors = 50, min_dist = 0.3, metric = "cosine", verbose = FALSE),
  nGenes = 3000,
  sampleCellsATAC = 17723,
  sampleCellsRNA = 4280,
  groupRNA = "cell_type",
  nameCell = "predictedCell_2",
  nameGroup = "predictedGroup_2",
  nameScore = "predictedScore_2")
#
#
#
#
plotEmbedding(ArchRProj = zhpf10, colorBy = "cellColData", name = "predictedGroup_1", embedding = "Tile_LM1_TSNE2")
plotEmbedding(ArchRProj = zhpf10, colorBy = "cellColData", name = "predictedGroup_2", embedding = "Tile_LM1_TSNE2")
#
# make one pseudo seurat file
x.num <- c(rep(1,800000),rep(0,800000),rep(2,800000))
head(x.num)
x.num.sample <- cbind(sample(x.num,length(zhpf10$cellNames)),
                      sample(x.num,length(zhpf10$cellNames)),
                      sample(x.num,length(zhpf10$cellNames)),
                      sample(x.num,length(zhpf10$cellNames)),
                      sample(x.num,length(zhpf10$cellNames)))
head(x.num.sample)
colnames(x.num.sample) <- c("F1","F2","F3","F4","F5")
rownames(x.num.sample) <- zhpf10$cellNames
#
zhpf10.seurat <- CreateSeuratObject(counts = t(x.num.sample), assay = "ATAC", project = "zhpf10", min.cells = 0, min.features = 0)
zhpf10.seurat
zhpf10.seurat <- NormalizeData(zhpf10.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
zhpf10.seurat <- FindVariableFeatures(zhpf10.seurat, selection.method = "vst", nfeatures = 4)
zhpf10.seurat <- ScaleData(zhpf10.seurat, features = rownames(zhpf10.seurat))
variablegene <- VariableFeatures(object = zhpf10.seurat)
zhpf10.seurat <- RunPCA(zhpf10.seurat, features = variablegene,npcs =2)
DimPlot(zhpf10.seurat)
#
x <- zhpf10@embeddings$Tile_LM1_TSNE2
y <- x$df
s.UMAP <- cbind(y$`IterativeLSI_Tile_LM1#TSNE_Dimension_1`,
                y$`IterativeLSI_Tile_LM1#TSNE_Dimension_2`)
#
head(s.UMAP )
rownames(s.UMAP ) <- rownames(y)
head(s.UMAP )
colnames(s.UMAP ) <- c("PC_1","PC_2")
head(s.UMAP )
class(s.UMAP )
#s.UMAP <- as.data.frame(s.UMAP)
zhpf10.seurat@reductions$pca@cell.embeddings <- s.UMAP 
DimPlot(zhpf10.seurat,label = T,repel = T,pt.size = 0.2)+NoLegend()
#
zhpf10.seurat$predictedGroup_1 <- zhpf10$predictedGroup_1
zhpf10.seurat$predictedGroup_2 <- zhpf10$predictedGroup_2
zhpf10.seurat$DoubletScore <- zhpf10$DoubletScore
zhpf10.seurat$DoubletEnrichment <- zhpf10$DoubletEnrichment
zhpf10.seurat$Tile_LM1.C3 <- zhpf10$Tile_LM1.C3
#
DimPlot(zhpf10.seurat,label = T,repel = T,pt.size = 0.2,group.by = "predictedGroup_1")+NoLegend()
DimPlot(zhpf10.seurat,label = T,repel = T,pt.size = 0.2,group.by = "predictedGroup_2")+NoLegend()
DimPlot(zhpf10.seurat,label = T,repel = T,pt.size = 0.2,group.by = "Tile_LM1.C3")+NoLegend()
#
FeaturePlot(zhpf10.seurat,features = "DoubletScore")
FeaturePlot(zhpf10.seurat,features = "DoubletEnrichment",max.cutoff = 8)
#
zhpf10 <- addImputeWeights(zhpf10,reducedDims = "IterativeLSI_Tile_LM1", corCutOff = 0.5)
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "glula", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "sox10", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "sox32", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
#
VlnPlot(hpf10.W,features = "dlb",group.by ="cell_type")+NoLegend()
VlnPlot(hpf10.W,features = "tbx1",group.by ="cell_type")+NoLegend()
VlnPlot(hpf10.W,features = "prdm1a",group.by ="cell_type")+NoLegend()
VlnPlot(hpf10.W,features = "dla",group.by ="cell_type")+NoLegend()
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "dla", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "dlb", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "neurod4", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "her13", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "sox2", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "sox19a", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "pnx", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "pnx", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "wnt11r", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "olig4", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "her3", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "hoxa9b", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "eve1", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "fn1b", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "msgn1", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "hsp90aa1.1", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "unc45b", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "klhl41b", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "nexn", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "mef2d", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "raraa", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "hnf1ba", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "sox10", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "tfap2a", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "zic2b", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "foxd3", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "tbx1", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "prdm1a", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "chic1", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "ripply1", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "cxcl12b", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "cxcl12a", 
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "fli1a",
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "aplnr2",
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
#
features <- list(NC1 = intersect(GRCz11.geneAnnotation$genes$symbol,c("foxd3", "sox10", "gdf6a", "lmx1bb")),
                 NC2 = intersect(GRCz11.geneAnnotation$genes$symbol,c("foxd3", "pax3a")))
#
zhpf10 <- addModuleScore(ArchRProj = zhpf10,useMatrix = "GeneScoreMatrix", features =features )
#
plotEmbedding(ArchRProj = zhpf10, colorBy = "cellColData",size = 0.0001, 
              name = "Module.NC1", embedding = "Tile_LM1_TSNE2",
              imputeWeights = getImputeWeights(zhpf10))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
#
addArchRThreads(threads = 8) 
markersGS <- getMarkerFeatures(
  ArchRProj = zhpf10, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Tile_LM1.C3",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
#
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$C6$name[1:20]
markerList$C7$name[1:20]
markerList$C27$name[1:20]
markerList$C28$name[1:20]
#
zhpf10.seurat@active.ident <- as.factor(zhpf10.seurat$Tile_LM1.C3)
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C1"="hpf10:Periderm")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C2"="hpf10:Prechordal plate")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C3"="hpf10:Neural plate posterior")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C4"="hpf10:Neural plate posterior")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C5"="hpf10:Neural crest")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C6"="hpf10:Doublets")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C7"="hpf10:Differentiating neurons?")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C8"="hpf10:Neural plate posterior")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C9"="hpf10:Neural plate posterior")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C10"="hpf10:Hindbrain")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C11"="hpf10:Midbrain")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C12"="hpf10:Diencephalon")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C13"="hpf10:Neural plate anterior")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C14"="hpf10:Neural anterior")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C15"="hpf10:Telencephalon")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C16"="hpf10:Neural crest")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C17"="hpf10:Tailbud mesoderm")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C18"="hpf10:Mesoderm adaxial cells")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C19"="hpf10:Tailbud spinal cord")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C20"="hpf10:Tailbud mesoderm")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C21"="hpf10:Tailbud mesoderm")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C22"="hpf10:Endoderm")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C23"="hpf10:Mesoderm lateral plate")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C24"="hpf10:Mesoderm lateral plate")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C25"="hpf10:Notochord")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C26"="hpf10:Tailbud mesoderm")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C27"="hpf10:Mesoderm lateral plate (tbx1+)")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C28"="hpf10:Mesoderm lateral plate (tbx1+)")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C29"="hpf10:YSL")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C30"="hpf10:Epidermal")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C31"="hpf10:Epidermal")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C32"="hpf10:Epidermal (gbx2+)")
zhpf10.seurat <- RenameIdents(zhpf10.seurat,"C33"="hpf10:Anterior neural ridge")
#
DimPlot(zhpf10.seurat,label = T,repel = T)+NoLegend()
#
zhpf10.seurat$ID <- zhpf10.seurat@active.ident
zhpf10$ID <- as.character(zhpf10.seurat$ID)
p <- plotBrowserTrack(ArchRProj = zhpf10, groupBy = "ID",
                      geneSymbol = "epcam", upstream = 50000,downstream = 50000,
                      ylim = c(0.001,0.999999))
grid::grid.newpage()
grid::grid.draw(p$epcam)
#
p <- plotBrowserTrack(ArchRProj = zhpf10, groupBy = "ID",
                      geneSymbol = "shhb", upstream = 10000,downstream = 3500,
                      ylim = c(0.001,0.999999))
grid::grid.newpage()
grid::grid.draw(p$shhb)
#
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "sox2",
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "epcam",
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "fn1b",
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "foxa3",
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "shhb",
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
#
getAvailableMatrices(zhpf10)
plotEmbedding(ArchRProj = zhpf10, colorBy = "GIM_Tile_LSI.M1.1", name = "sox2",
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GIM_Tile_LSI.M1.1", name = "epcam",
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GIM_Tile_LSI.M1.1", name = "fn1b",
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GIM_Tile_LSI.M1.1", name = "foxa3",
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
#
save.image("zhpf10.RData")
#
VlnPlot(hpf10.W,features = "shhb",group.by ="cell_type")+NoLegend()
VlnPlot(hpf10.F,features = "shhb",group.by ="cell_type")+NoLegend()
#
zhpf10.seurat$DoubletScore.log2 <- log2(zhpf10.seurat$DoubletScore+1)
zhpf10.seurat$DoubletEnrichment.log2 <- log2(zhpf10.seurat$DoubletEnrichment+1)
VlnPlot(zhpf10.seurat,features = "DoubletScore.log2",group.by ="ID",pt.size = 0)+NoLegend()
VlnPlot(zhpf10.seurat,features = "DoubletEnrichment.log2",group.by ="ID",pt.size = 0)+NoLegend()
#
#
# correlation analysis between scRNA-seq and scATAC-seq data set 
hpf10
table(hpf10$cell_type)
library(dplyr)
hpf10@active.ident <- as.factor(hpf10$cell_type)
hpf10.markers <- FindAllMarkers(hpf10, only.pos = TRUE, min.pct = 0.25)
hpf10.markers.top <- hpf10.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(hpf10, features = unique(hpf10.markers.top$gene)) + NoLegend()
##
#
hpf10@active.ident <- as.factor(hpf10$cell_type)
table(hpf10@active.ident)
hpf10 <- RenameIdents(hpf10,"anterior neural ridge"="Anterior neural ridge")
hpf10 <- RenameIdents(hpf10,"diencephalon"="Diencephalon")
hpf10 <- RenameIdents(hpf10,"diencephalon (aplnr2+)"="Diencephalon (aplnr2+)")
hpf10 <- RenameIdents(hpf10,"differentiating neurons"="Differentiating neurons")
hpf10 <- RenameIdents(hpf10,"endoderm"="Endoderm")
hpf10 <- RenameIdents(hpf10,"epidermal"="Epidermal")
hpf10 <- RenameIdents(hpf10,"epidermal (gbx2+)"="Epidermal (gbx2+)")
hpf10 <- RenameIdents(hpf10,"hindbrain"="Hindbrain")
hpf10 <- RenameIdents(hpf10,"mesoderm adaxial cells"="Mesoderm adaxial cells")
hpf10 <- RenameIdents(hpf10,"mesoderm lateral plate (tbx1+)"="Mesoderm lateral plate (tbx1+)")
hpf10 <- RenameIdents(hpf10,"mesoderm lateral plate"="Mesoderm lateral plate")
hpf10 <- RenameIdents(hpf10,"midbrain"="Midbrain")
hpf10 <- RenameIdents(hpf10,"neural anterior"="Neural anterior")
hpf10 <- RenameIdents(hpf10,"neural crest"="Neural crest")
hpf10 <- RenameIdents(hpf10,"neural plate anterior"="Neural plate anterior")
hpf10 <- RenameIdents(hpf10,"neural plate posterior"="Neural plate posterior")
hpf10 <- RenameIdents(hpf10,"notochord"="Notochord")
hpf10 <- RenameIdents(hpf10,"periderm"="Periderm")
hpf10 <- RenameIdents(hpf10,"prechordal plate"="Prechordal plate")
hpf10 <- RenameIdents(hpf10,"tailbud mesoderm"="Tailbud mesoderm")
hpf10 <- RenameIdents(hpf10,"tailbud spinal cord"="Tailbud spinal cord")
hpf10 <- RenameIdents(hpf10,"telencephalon"="Telencephalon")
#
hpf10$cell_type_FORcor <- hpf10@active.ident
#
levels(hpf10)
length(levels(hpf10))
ID.Keep <- levels(hpf10)[1:(length(levels(hpf10))-1)]
#
ID.Keep2 <- c("Mesoderm lateral plate","Telencephalon","Tailbud spinal cord","Tailbud mesoderm",
              "Prechordal plate","Periderm","Notochord","Neural plate posterior","Neural plate anterior",
              "Neural crest","Neural anterior","Midbrain",
              "Mesoderm lateral plate (tbx1+)","Mesoderm adaxial cells",
              "Hindbrain","Epidermal (gbx2+)","Epidermal","Endoderm",
              "Differentiating neurons","Diencephalon (aplnr2+)",
              "Diencephalon","Anterior neural ridge")
#
hpf10.2 <- subset(hpf10,idents=ID.Keep2)
#
hpf10.2@active.ident <- as.factor(hpf10.2$cell_type_FORcor)
levels(hpf10.2) <- cluster_rank 
hpf10.2.markers <- FindAllMarkers(hpf10.2, only.pos = TRUE, min.pct = 0.25)
hpf10.2.markers.top <- hpf10.2.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
DoHeatmap(hpf10.2, features = unique(hpf10.2.markers.top$gene)) + NoLegend()
##
#
zhpf10.seurat@active.ident <- as.factor(zhpf10.seurat$ID)
zhpf10.seurat$ID2 <- unlist(strsplit(as.character(zhpf10.seurat$ID),":"))[seq(2,2*17723,2)] ##removing hpf10: tag in names
zhpf10.seurat@active.ident <- as.factor(zhpf10.seurat$ID2)
#
zhpf10$ID2 <- as.character(zhpf10.seurat$ID2)
zhpf10$ID.standby <- as.character(zhpf10.seurat$ID.standby)
#
levels(zhpf10.seurat)
#
zhpf10.seurat.FORcor <- subset(zhpf10.seurat,idents=c("Doublets","YSL"),invert=T)
levels(zhpf10.seurat.FORcor)
#
Get.genescore.mat <- getMatrixFromProject(
  ArchRProj = zhpf10,
  useMatrix = "GeneScoreMatrix")
rownames(Get.genescore.mat@assays@data$GeneScoreMatrix) <- Get.genescore.mat@elementMetadata$name
Get.genescore.mat.matrix <- Get.genescore.mat@assays@data$GeneScoreMatrix
head(Get.genescore.mat.matrix)
#
zhpf10.meta.FORcor <- as.data.frame(zhpf10@cellColData)
#
Get.genescore.mat.matrix
cogenes <- intersect(rownames(Get.genescore.mat.matrix),rownames(hpf10))
cogenes
Get.genescore.mat.matrix.suerat <- CreateSeuratObject(counts = Get.genescore.mat.matrix[cogenes,], assay = "RNA", project = "GIMx")
rm(Get.genescore.mat.matrix) 
gc()
Get.genescore.mat.matrix.suerat <- NormalizeData(Get.genescore.mat.matrix.suerat)
Get.genescore.mat.matrix.suerat$ID <- zhpf10.meta.FORcor$ID2
Get.genescore.mat.matrix.suerat <- subset(Get.genescore.mat.matrix.suerat,idents=c("Doublets","YSL"),invert=T)
#
Get.genescore.mat.matrix.suerat@active.ident <- as.factor(Get.genescore.mat.matrix.suerat$ID)
table(Get.genescore.mat.matrix.suerat@active.ident)
length(table(Get.genescore.mat.matrix.suerat@active.ident))
#
levels(Get.genescore.mat.matrix.suerat) <- cluster_rank
Get.genescore.mat.matrix.suerat.averages <- AverageExpression(Get.genescore.mat.matrix.suerat, slot = 'data', return.seurat = TRUE)
Get.genescore.mat.matrix.suerat.averages
DoHeatmap(Get.genescore.mat.matrix.suerat.averages, features = unique(top3$gene),draw.lines = F) + NoLegend()+
  scale_fill_gradientn(colors = paletteContinuous("blueYellow"))
#
########
#
hpf10.2$cell_type_FORcor <- as.character(hpf10.2$cell_type_FORcor)
table(hpf10.2$cell_type_FORcor)
#
hpf10.2.averages <- AverageExpression(hpf10.2, slot = 'data', return.seurat = TRUE)
hpf10.2.averages
levels(hpf10.2.averages) <- cluster_rank
top500 <- hpf10.2.markers %>% group_by(cluster) %>% top_n(n = 37, wt = avg_log2FC)
length(unique(top3$gene))
DoHeatmap(hpf10.2.averages, features = unique(top500$gene),draw.lines = F) + NoLegend()+
  scale_fill_gradientn(colors = paletteContinuous("horizonExtra"))
#
#
Get.genescore.mat.matrix.suerat.averages.sub <- GetAssayData(Get.genescore.mat.matrix.suerat.averages, slot = "scale.data")[unique(top3$gene),]
dim(Get.genescore.mat.matrix.suerat.averages.sub)
colnames(Get.genescore.mat.matrix.suerat.averages.sub) <- paste("ATAC_",colnames(Get.genescore.mat.matrix.suerat.averages.sub),sep = "")
head(Get.genescore.mat.matrix.suerat.averages.sub)
hpf10.2.averages.sub <- GetAssayData(hpf10.2.averages, slot = "scale.data")[unique(top3$gene),]
dim(hpf10.2.averages.sub)
colnames(hpf10.2.averages.sub) <- paste("RNA_",colnames(hpf10.2.averages.sub),sep = "")
#
cc1 <- cbind(Get.genescore.mat.matrix.suerat.averages.sub,
             hpf10.2.averages.sub)
cc1.t <- cor(cc1,method = c("pearson"))
dim(cc1.t)
head(cc1.t)
#
pheatmap(cc1.t[23:44,1:22],border_color = NA,
         cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100))
###
pheatmap(cc1.t[23:44,1:22],border_color = NA,
         cluster_rows = T,cluster_cols = T,
         color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100))
###
pheatmap(cc1.t[23:44,1:22],border_color = NA,
         cluster_rows = F,cluster_cols = T,
         color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100))
#
pheatmap(cc1.t[1:22,1:22],border_color = NA,
         cluster_rows = T,cluster_cols = T,
         color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100))
pheatmap(cc1.t[23:44,23:44],border_color = NA,
         cluster_rows = T,cluster_cols = T,
         color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100))
#
out <- pheatmap::pheatmap(cor(hpf10.2.averages.sub))
out
rownames(hpf10.2.averages.sub[out$tree_row[["order"]],])
colnames(hpf10.2.averages.sub[,out$tree_col[["order"]]])
#
#
#cluster_rank <- unlist(strsplit(colnames(hpf10.2.averages.sub[,out$tree_col[["order"]]]),"RNA_"))[seq(2,2*22,2)]
cluster_rank 
#
#
#
# 
#########
# color set for hpf10h
# 
levels(zhpf10.seurat)
#
zhpf10.seurat@active.ident <- as.factor(zhpf10.seurat$ID.standby)
levels(zhpf10.seurat) 
cell.type.color.set2 <- c("#5E4FA2",# hpf10:Diencephalon (aplnr2+)
                          "#FDC086",# hpf10:Epidermal
                          "#66C2A5", #hpf10:Anterior neural ridge
                          "#DF8E32", #hpf10:Diencephalon #
                          "#3288BD", #hpf10:Differentiating neurons?
                          "#F7F7F7", #hpf10:Doublets
                          "#00BFFF", #hpf10:Endoderm
                          "#87CEFA", #hpf10:Epidermal (gbx2+) ##313695
                          "#DFC27D", #hpf10:Hindbrain
                          "#9933CC", #hpf10:Mesoderm adaxial cells
                          "#999999", #hpf10:Mesoderm lateral plate
                          "#48D1CC", #hpf10:Mesoderm lateral plate (tbx1+)
                          "#FF33CC", #hpf10:Midbrain
                          "#63D7D2", #hpf10:Neural anterior
                          "#FBB4AE", #hpf10:Neural crest
                          "#B3EE3A", #hpf10:Neural plate anterior
                          "#73C467", #hpf10:Neural plate posterior
                          "#9E0142", #hpf10:Notochord
                          "red", #hpf10:Periderm
                          "#9970AB", #hpf10:Prechordal plate
                          "#6495ED", #hpf10:Tailbud mesoderm
                          "#D6604D", #hpf10:Tailbud spinal cord
                          "#B35806", #hpf10:Telencephalon
                          "#0000FF" # hpf10:YSL
)
#
cell.type.color.set2.meta <- as.data.frame(cbind(levels(zhpf10.seurat),cell.type.color.set2))
cell.type.color.set2.meta
colnames(cell.type.color.set2.meta) <- c("cluster","color")
rownames(cell.type.color.set2.meta) <- cell.type.color.set2.meta$cluster
#
DimPlot(zhpf10.seurat,label = T,repel = T,cols = cell.type.color.set2,raster = T,pt.size = 0.00001)
DimPlot(zhpf10.seurat,label = T,repel = T,cols = cell.type.color.set2,pt.size = 0.00001)+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
#
p <- DimPlot(zhpf10.seurat,label = F,repel = T,cols = cell.type.color.set2,pt.size = 0.00001)+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/OBO/hpf10/Fig_save/hpf10h.ID.png",plot=p,device="png",dpi=300,units = "cm",width = 12,height = 12)
##
#
#
zhpf10.seurat$predictedGroup_1x <- as.character(zhpf10$predictedGroup_1x)
zhpf10.seurat@active.ident <- as.factor(zhpf10.seurat$predictedGroup_1x)
levels(zhpf10.seurat) 
##
levels(zhpf10.seurat) <- c("diencephalon (aplnr2+)",
                  "epidermal",
                  "anterior neural ridge",
                  "diencephalon",
                  "differentiating neurons",
                  "endoderm",
                  "epidermal (gbx2+)",
                  "hindbrain",
                  "mesoderm adaxial cells",
                  "mesoderm lateral plate",
                  "mesoderm lateral plate (tbx1+)",
                  "midbrain",
                  "neural anterior",
                  "neural crest",
                  "neural plate anterior",
                  "neural plate posterior",
                  "notochord",
                  "periderm",
                  "prechordal plate",
                  "tailbud mesoderm",
                  "tailbud spinal cord",
                  "telencephalon")
#
levels(zhpf10.seurat)
#
cell.type.color.set2.meta$color[c(-6,-24)]
DimPlot(zhpf10.seurat,label = T,repel = T,cols = cell.type.color.set2.meta$color[c(-6,-24)],pt.size = 0.00001)+NoLegend()
#
DimPlot(zhpf10.seurat,label = T,repel = T,cols = cell.type.color.set2.meta$color[c(-6,-24)],pt.size = 0.00001)
#
p <- DimPlot(zhpf10.seurat,label = F,repel = T,cols = cell.type.color.set2.meta$color[c(-6,-24)],pt.size = 0.00001)+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/OBO/hpf10/Fig_save/hpf10h.predicted.ID.png",plot=p,device="png",dpi=300,units = "cm",width = 12,height = 12)
##
####
# for correlation analysis between gene score inferred by chromatin accessibility and gene expression
# using top 500 DEG genes in hpf10
top500 <- hpf10.2.markers %>% group_by(cluster) %>% top_n(n = 37, wt = avg_log2FC)
length(unique(top3$gene))
#
cell.type.color.set2.meta2 <- cell.type.color.set2.meta
dim(cell.type.color.set2.meta2)
rownames(cell.type.color.set2.meta2) <- unlist(strsplit(rownames(cell.type.color.set2.meta2),":"))[seq(2,2*24,2)]
cell.type.color.set2.meta2[levels(hpf10.2.averages),]$color
#
DoHeatmap(Get.genescore.mat.matrix.suerat.averages, features = unique(top3$gene)[1:500], draw.lines =F,
          group.colors=cell.type.color.set2.meta2[levels(hpf10.2.averages),]$color) + NoLegend()+
  scale_fill_gradientn(colors = paletteContinuous("blueYellow"))
#######
DoHeatmap(hpf10.2.averages, features = unique(top3$gene)[1:500], draw.lines =F,
          group.colors=cell.type.color.set2.meta2[levels(hpf10.2.averages),]$color) + NoLegend()+
  scale_fill_gradientn(colors = paletteContinuous("horizonExtra"))
#
cogenes <- intersect(rownames(Get.genescore.mat.matrix.suerat.averages@assays$RNA@scale.data),
                     rownames(hpf10.2.averages@assays$RNA@scale.data))
zzz1 <- Get.genescore.mat.matrix.suerat.averages@assays$RNA@scale.data[unique(top3$gene)[1:500],]
zzz2 <- hpf10.2.averages@assays$RNA@scale.data[unique(top3$gene)[1:500],]
#
pheatmap(cor(t(cbind(zzz1,zzz2))),cluster_rows = F,cluster_cols = F,color  = paletteContinuous("horizonExtra"))
#
cor.tmp.x <- c()
for (i in 1:500) {
  cor.tmp <- cor.test(as.numeric(zzz1[i,]),as.numeric(zzz2[i,]))
  cor.tmp.x <- c(cor.tmp.x,cor.tmp$estimate)
}
boxplot(cor.tmp.x,outline=F)
median(cor.tmp.x)
mean(cor.tmp.x)
#
zzz1 <- Get.genescore.mat.matrix.suerat.averages@assays$RNA@scale.data[cogenes,]
zzz2 <- hpf10.2.averages@assays$RNA@scale.data[cogenes,]
#
cor.tmp.x <- c()
for (i in 1:length(cogenes)) {
  cor.tmp <- cor.test(as.numeric(zzz1[i,]),as.numeric(zzz2[i,]))
  cor.tmp.x <- c(cor.tmp.x,cor.tmp$estimate)
}
boxplot(cor.tmp.x,outline=F)
median(cor.tmp.x)
cor.tmp.x <- as.data.frame(cor.tmp.x)$cor.tmp.x
boxplot(na.omit(cor.tmp.x))
median(na.omit(cor.tmp.x))
################
################
###############
# showing some classical marker genes in the UMAP space
# based on Gene score matrix 
p <-plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix",size = 0.0001, name = "Sox2", 
                  embedding = "Tile_LM1_TSNE2",
                  imputeWeights = getImputeWeights(zhpf10),
                  pal =paletteContinuous("blueYellow"))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/OBO/hpf10/Fig_save/Sox2.genescore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix",size = 0.0001, name = "Epcam",
                  embedding = "Tile_LM1_TSNE2",
                  imputeWeights = getImputeWeights(zhpf10),
                  pal =paletteContinuous("blueYellow"))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/OBO/hpf10/Fig_save/Epcam.genescore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix",size = 0.0001, name = "foxa3",
                  embedding = "Tile_LM1_TSNE2",
                  imputeWeights = getImputeWeights(zhpf10),
                  pal =paletteContinuous("blueYellow"))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/OBO/hpf10/Fig_save/foxa3.genescore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix",size = 0.0001, name = "Fn1b",
                  embedding = "Tile_LM1_TSNE2",
                  imputeWeights = getImputeWeights(zhpf10),
                  pal =paletteContinuous("blueYellow"))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/OBO/hpf10/Fig_save/Fn1b.genescore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix",size = 0.0001, name = "shhb",
                  embedding = "Tile_LM1_TSNE2",
                  imputeWeights = getImputeWeights(zhpf10),
                  pal =paletteContinuous("blueYellow"))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/OBO/hpf10/Fig_save/shhb.genescore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
#################
# based on integrated gene expression matrix
p <-plotEmbedding(ArchRProj = zhpf10, colorBy = "GIM_Tile_LSI.M1.1", name = "Sox2",
                  size = 0.0001, embedding = "Tile_LM1_TSNE2", 
                  imputeWeights = getImputeWeights(zhpf10))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/OBO/hpf10/Fig_save/Sox2.GIM.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf10, colorBy = "GIM_Tile_LSI.M1x.1", name = "Epcam",
                  size = 0.0001, embedding = "Tile_LM1_TSNE2", 
                  imputeWeights = getImputeWeights(zhpf10))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/OBO/hpf10/Fig_save/Epcam.GIM.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf10, colorBy = "GIM_Tile_LSI.M1x.1", name = "foxa3",
                  size = 0.0001, embedding = "Tile_LM1_TSNE2", 
                  imputeWeights = getImputeWeights(zhpf10))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/OBO/hpf10/Fig_save/foxa3.GIM.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf10, colorBy = "GIM_Tile_LSI.M1x.1", name = "fn1b",
                  size = 0.0001, embedding = "Tile_LM1_TSNE2", 
                  imputeWeights = getImputeWeights(zhpf10))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/OBO/hpf10/Fig_save/fn1b.GIM.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf10, colorBy = "GIM_Tile_LSI.M1x.1", name = "shhb",
                  size = 0.0001, embedding = "Tile_LM1_TSNE2", 
                  imputeWeights = getImputeWeights(zhpf10))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/OBO/hpf10/Fig_save/shhb.GIM.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
#
#
#
######################
######################
######################
######################
######################
# track view in Figure 2G
#cell.type.color.set2.meta2[c(levels(hpf10.2.averages),"YSL"),]$color
p <- plotBrowserTrack(ArchRProj = zhpf10, groupBy = "ID2", geneSymbol = "sox2", 
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 1),ylim = c(0.001,0.999999),
                      upstream = 9000,downstream = 8000,
                      pal =cell.type.color.set2.meta2[c(levels(hpf10.2.averages),"YSL"),]$color,
                      useGroups =c(cluster_rank,"YSL"))
grid::grid.newpage()
grid::grid.draw(p$sox2)
plotPDF(plotList = p, name = "sox2.pdf", ArchRProj = zhpf10, addDOC = FALSE, width = 5, height = 8)
#
p <- plotBrowserTrack(ArchRProj = zhpf10, groupBy = "ID2", geneSymbol = "actb1", 
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 1),ylim = c(0.001,0.999999),
                      upstream = 5000,downstream = 2000,
                      pal =cell.type.color.set2.meta2[c(levels(hpf10.2.averages),"YSL"),]$color,
                      useGroups =c(cluster_rank,"YSL"))
grid::grid.newpage()
grid::grid.draw(p$actb1)
plotPDF(plotList = p, name = "actb1.pdf", ArchRProj = zhpf10, addDOC = FALSE, width = 5, height = 8)
#
p <- plotBrowserTrack(ArchRProj = zhpf10, groupBy = "ID2", geneSymbol = "shhb",
                      minCells = 25,useMatrix="GIM_Tile_LSI.M1.1",scCellsMax=400,
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 1),ylim = c(0.001,0.9999999),
                      upstream = 10000,downstream = 3500,
                      pal =cell.type.color.set2.meta2[c(levels(hpf10.2.averages),"YSL"),]$color,
                      useGroups =c(cluster_rank,"YSL")
)
grid::grid.newpage()
grid::grid.draw(p$shhb)
plotPDF(plotList = p, name = "shhb.pdf", ArchRProj = zhpf10, addDOC = FALSE, width = 5, height = 8)
#
p <- plotBrowserTrack(ArchRProj = zhpf10, groupBy = "ID2", geneSymbol = "pax2a", 
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 1),ylim = c(0.001,0.999999),
                      upstream = 20000,downstream = 40000,
                      pal =cell.type.color.set2.meta2[c(levels(hpf10.2.averages),"YSL"),]$color,
                      useGroups =c(cluster_rank,"YSL"))
grid::grid.newpage()
grid::grid.draw(p$pax2a)
plotPDF(plotList = p, name = "pax2a.pdf", ArchRProj = zhpf10, addDOC = FALSE, width = 5, height = 8)
#
p <- plotBrowserTrack(ArchRProj = zhpf10, groupBy = "ID2", geneSymbol = "hesx1", 
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 1),ylim = c(0.001,0.999999),
                      upstream = 2500,downstream = 8000,
                      pal =cell.type.color.set2.meta2[c(levels(hpf10.2.averages),"YSL"),]$color,
                      useGroups =c(cluster_rank,"YSL"))
grid::grid.newpage()
grid::grid.draw(p$hesx1)
plotPDF(plotList = p, name = "hesx1.pdf", ArchRProj = zhpf10, addDOC = FALSE, width = 5, height = 8)
#
p <- plotBrowserTrack(ArchRProj = zhpf10, groupBy = "ID2", geneSymbol = "ripply1", 
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 1),ylim = c(0.001,0.999999),
                      upstream = 8000,downstream = 4000,
                      pal =cell.type.color.set2.meta2[c(levels(hpf10.2.averages),"YSL"),]$color,
                      useGroups =c(cluster_rank,"YSL"))
grid::grid.newpage()
grid::grid.draw(p$ripply1)
plotPDF(plotList = p, name = "ripply1.pdf", ArchRProj = zhpf10, addDOC = FALSE, width = 5, height = 8)
#
p <- plotBrowserTrack(ArchRProj = zhpf10, groupBy = "ID2", geneSymbol = "rx3", 
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 1),ylim = c(0.001,0.999999),
                      upstream = 2500,downstream = 6500,
                      pal =cell.type.color.set2.meta2[c(levels(hpf10.2.averages),"YSL"),]$color,
                      useGroups =c(cluster_rank,"YSL"))
grid::grid.newpage()
grid::grid.draw(p$rx3)
plotPDF(plotList = p, name = "rx3.pdf", ArchRProj = zhpf10, addDOC = FALSE, width = 5, height = 8)
#
p <- plotBrowserTrack(ArchRProj = zhpf10, groupBy = "ID2", geneSymbol = "ctslb", 
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 1),ylim = c(0.001,0.999999),
                      upstream = 8000,downstream = 1500,
                      pal =cell.type.color.set2.meta2[c(levels(hpf10.2.averages),"YSL"),]$color,
                      useGroups =c(cluster_rank,"YSL"))
grid::grid.newpage()
grid::grid.draw(p$ctslb)
plotPDF(plotList = p, name = "ctslb.pdf", ArchRProj = zhpf10, addDOC = FALSE, width = 5, height = 8)
#
p <- plotBrowserTrack(ArchRProj = zhpf10, groupBy = "ID2", geneSymbol = "tbx6", 
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 1),ylim = c(0.001,0.999999),
                      upstream = 2000,downstream = 14000,
                      pal =cell.type.color.set2.meta2[c(levels(hpf10.2.averages),"YSL"),]$color,
                      useGroups =c(cluster_rank,"YSL"))
grid::grid.newpage()
grid::grid.draw(p$tbx6)
plotPDF(plotList = p, name = "tbx6.pdf", ArchRProj = zhpf10, addDOC = FALSE, width = 5, height = 8)
#
p <- plotBrowserTrack(ArchRProj = zhpf10, groupBy = "ID2", geneSymbol = "Tp63", 
                      plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 1),ylim = c(0.001,0.999999),
                      upstream = 12000,downstream = 135000,
                      pal =cell.type.color.set2.meta2[c(levels(hpf10.2.averages),"YSL"),]$color,
                      useGroups =c(cluster_rank,"YSL"))
grid::grid.newpage()
grid::grid.draw(p$tp63)
plotPDF(plotList = p, name = "tp63.pdf", ArchRProj = zhpf10, addDOC = FALSE, width = 5, height = 8)
#
p <- plotBrowserTrack(ArchRProj = zhpf10, groupBy = "ID2", geneSymbol = "sox32", plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 1),ylim = c(0.001,0.999999),
                      upstream = 12000,downstream = 15000,pal =cell.type.color.set2.meta2[c(levels(hpf10.2.averages),"YSL"),]$color,useGroups =c(cluster_rank,"YSL"))
grid::grid.newpage()
grid::grid.draw(p$sox32)
plotPDF(plotList = p, name = "sox32.pdf", ArchRProj = zhpf10, addDOC = FALSE, width = 5, height = 8)
##################################################################################################
##########
# cell ratio correlation analysis between scATAC-seq and scRNA-seq
table(zhpf10.seurat.FORcor$ID2)
table(zhpf10.seurat.FORcor$ID)
table(hpf10.2$cell_type_FORcor)
#
#
cellratio.metadata <- as.data.frame(cbind(table(zhpf10.seurat.FORcor$ID2),table(hpf10.2$cell_type_FORcor)))
head(cellratio.metadata)
dim(cellratio.metadata)
#
cellratio.metadata.tmp <- cellratio.metadata[rev(cluster_rank),]
dim(cellratio.metadata.tmp)
cellratio.metadata.tmp
#
cellratio.metadata.tmp[,1] <- cellratio.metadata.tmp[,1]/sum(cellratio.metadata.tmp[,1])
cellratio.metadata.tmp[,2] <- cellratio.metadata.tmp[,2]/sum(cellratio.metadata.tmp[,2])
#
plot(cellratio.metadata.tmp[,1],cellratio.metadata.tmp[,2])
#
x.cor <- cor(cellratio.metadata.tmp[,1],cellratio.metadata.tmp[,2],method = c("pearson"))
x.cor <- cor.test(cellratio.metadata.tmp[,1],cellratio.metadata.tmp[,2],method = c("pearson"))
x.cor ## 0.9738682 ## p-value = 2.406e-14
#
head(cellratio.metadata.tmp)
cellratio.metadata.tmp2 <- cbind(as.numeric(cellratio.metadata.tmp[,1]),as.numeric(cellratio.metadata.tmp[,2]))
cellratio.metadata.tmp2 <- data.frame(cellratio.metadata.tmp2)
colnames(cellratio.metadata.tmp2) <- c("ATAC","RNA")
head(cellratio.metadata.tmp2)
#
ggplot(cellratio.metadata.tmp2, aes(ATAC,RNA))+
  geom_point(size=2)+
  geom_point(colour=cell.type.color.set2.meta2[levels(hpf10.2.averages),]$color)
#
ggplot(cellratio.metadata.tmp2, aes(ATAC,RNA))+
  geom_point(size=4)+
  geom_smooth(method = "lm")+
  geom_point(colour=cell.type.color.set2.meta2[levels(hpf10.2.averages),]$color)+theme_test()
#
hist(zhpf10$predictedScore_1x,breaks = 40,col = "red")
#
############### 
# Peak calling
addArchRThreads(threads = 1) 
zhpf10 <- addGroupCoverages(ArchRProj = zhpf10, groupBy = "ID")
#
zhpf10 <- addReproduciblePeakSet(
  ArchRProj = zhpf10, 
  groupBy = "ID",
  genomeSize = 1.5e9,
  excludeChr = c("MT",paste("chr",chr_rm$V1,sep = "")), 
  genomeAnnotation = genomeAnnotation,
  geneAnnotation = GRCz11.geneAnnotation,
  pathToMacs2 = "/home/sunkeyong/.local/bin/macs2", force=T
)
#
zhpf10.peakset <- getPeakSet(zhpf10)
#
zhpf10 <- addPeakMatrix(ArchRProj = zhpf10, force =T )
# Motif annotation
# loading zebrafish motif collected from DANIO-CODE
load("/home/sunkeyong/ZEPA/Custom_danRer11_Motif/danRer11.Motif.RData")
#
zhpf10 <- addMotifAnnotations(ArchRProj = zhpf10, motifPWMs=danRer11.Motif, name = "Motif")
#
zhpf10 <- addBgdPeaks(zhpf10)
#
zhpf10 <- addDeviationsMatrix(
  ArchRProj = zhpf10, 
  peakAnnotation = "Motif",
  force = TRUE
)
#
getAvailableMatrices(zhpf10)
#
plotVarDev <- getVarDeviations(zhpf10, name = "MotifMatrix", plot = TRUE)
plotVarDev
#
#
motifs <- c("tfap","TP63","sox","fox","myf")
motifs <- c("myf","my")
motifs <- c("zic","hoxb9a")
motifs <- c("zeb1b")
markerMotifs <- getFeatures(zhpf10, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs
#
plotEmbedding(ArchRProj = zhpf10, colorBy = "MotifMatrix", name = "deviations:tp63",
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "MotifMatrix", name = "deviations:tfap2c",
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "MotifMatrix", name = "deviations:sox10",
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "MotifMatrix", name = "deviations:sox2",
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "MotifMatrix", name = "deviations:myod1",
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "MotifMatrix", name = "deviations:myog",
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "MotifMatrix", name = "deviations:myf6",
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "MotifMatrix", name = "deviations:zic1",
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "MotifMatrix", name = "deviations:hoxb9a",
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "MotifMatrix", name = "deviations:zeb1b",
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = "zeb1b",
              embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))
#
p <- plotBrowserTrack(ArchRProj = zhpf10, groupBy = "ID2",
                      geneSymbol = "zeb1b", upstream = 50000,downstream = 50000,
                      ylim = c(0.001,0.999999))
grid::grid.newpage()
grid::grid.draw(p$zeb1b)
#
#
#
seGroupMotif <- getGroupSE(ArchRProj = zhpf10, useMatrix = "MotifMatrix", groupBy = "ID2")
seGroupMotif
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs
corGSM_MM <- correlateMatrices(
  ArchRProj = zhpf10,
  useMatrix1 = "GeneScoreMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "IterativeLSI_Tile_LM1x",
  dimsToUse = 2:100
)
corGSM_MM
#
corGIM_MM <- correlateMatrices(
  ArchRProj = zhpf10,
  useMatrix1 = "GIM_Tile_LSI.M1x.1",
  useMatrix2 = "MotifMatrix",
  reducedDims = "IterativeLSI_Tile_LM1x",
  dimsToUse = 2:100
)
corGIM_MM
#
corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])
#
View(data.frame(corGSM_MM))
#
p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )
#
p
#
corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])
#
View(data.frame(corGIM_MM))
p <- ggplot(data.frame(corGIM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Expression") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGIM_MM$maxDelta)*1.05)
  )
p
#
#
#
#
#
zhpf10
zhpf10 <- addPeak2GeneLinks(
  ArchRProj = zhpf10,
  useMatrix = "GIM_Tile_LSI.M1x.1",
  reducedDims = "IterativeLSI_Tile_LM1x",
  dimsToUse = 2:100,
  cellsToUse = zhpf10$cellNames,
  k = 150,
  knnIteration = 500
)
#
p2g <- getPeak2GeneLinks(
  ArchRProj = zhpf10,
  corCutOff = 0.4,
  resolution = 1,
  returnLoops = FALSE,
  varCutOffATAC = 0.2,
  varCutOffRNA = 0.2
)
p2g
#
p <- plotBrowserTrack(ArchRProj = zhpf10, groupBy = "ID2",loops = getPeak2GeneLinks(zhpf10),
                      geneSymbol = "zeb1b", upstream = 50000,downstream = 50000,
                      ylim = c(0.001,0.999999))
grid::grid.newpage()
grid::grid.draw(p$zeb1b)
#
p <- plotPeak2GeneHeatmap(ArchRProj = zhpf10, groupBy = "ID2")+NoLegend()
p
#
#
##########################
# plot markers for YSL
markersGS <- getMarkerFeatures(
  ArchRProj = zhpf10, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "ID2",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
####
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$YSL$name[1:20]
# glula
#
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix",size = 0.0001, name = "glula",
              pal =paletteContinuous("blueYellow"), embedding = "Tile_LM1_TSNE2", imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix",size = 0.0001, name = "aldob",
              pal =paletteContinuous("blueYellow"), embedding = "Tile_LM1_TSNE2", imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix",size = 0.0001, name = "fbp1b",
              pal =paletteContinuous("blueYellow"), embedding = "Tile_LM1_TSNE2", imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix",size = 0.0001, name = "apoeb",
              pal =paletteContinuous("blueYellow"), embedding = "Tile_LM1_TSNE2", imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix",size = 0.0001, name = "apoa1b",
              pal =paletteContinuous("blueYellow"), embedding = "Tile_LM1_TSNE2", imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix",size = 0.0001, name = "afp4",
              pal =paletteContinuous("blueYellow"), embedding = "Tile_LM1_TSNE2", imputeWeights = getImputeWeights(zhpf10))
plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix",size = 0.0001, name = "fetub",
              pal =paletteContinuous("blueYellow"), embedding = "Tile_LM1_TSNE2", imputeWeights = getImputeWeights(zhpf10))

features.in.zhpf10 <- c("glula","aldob","fbp1b","apoeb","apobb.1","gata5","apobb.1","g6pca.2")
##
for (i in 1:length(features.in.zhpf10)) {
  p <- plotEmbedding(ArchRProj = zhpf10, colorBy = "GeneScoreMatrix", name = features.in.zhpf10[i],rastr = F,
                     embedding = "Tile_LM1_TSNE2",imputeWeights = getImputeWeights(zhpf10))+NoLegend()+
    theme(plot.title = element_blank())+
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
    theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
    theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
    theme(panel.background = element_rect( colour = "black", size = 0))+
    theme(panel.border = element_blank())
  ggsave(paste("/home/sunkeyong/ZEPA/OBO/hpf10/Fig_save_YSL/zhpf10.",features.in.zhpf10[i],".genescore.png",sep = ""),plot=p,device="png",dpi=300,units = "cm",width = 15,height = 15)
}
#
#
p <- plotBrowserTrack(ArchRProj = zhpf10, groupBy = "ID2", geneSymbol = "ripply1", plotSummary = c("bulkTrack", "geneTrack"),sizes = c(10, 1),ylim = c(0.001,0.999999),
                      upstream = 8000,downstream = 4000,pal =cell.type.color.set2.meta2[c(levels(hpf10.2.averages),"YSL"),]$color,useGroups =c(cluster_rank,"YSL"))
grid::grid.newpage()
grid::grid.draw(p$ripply1)
plotPDF(plotList = p, name = "ripply1.pdf", ArchRProj = zhpf10, addDOC = FALSE, width = 5, height = 8)
#
VlnPlot(hpf10.F,features = "dla",group.by = "cell_type")+NoLegend()
VlnPlot(hpf10.W,features = "dla",group.by = "cell_type")+NoLegend()
#
# analyze the shha shhb in the notochord for the reviewer 1
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf9.RData")
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf8.W.RData")
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf8.F.RData")
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf7.RData")
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf6.RData")
#
VlnPlot(hpf10.F,features = "shha",group.by = "cell_type")+NoLegend()
VlnPlot(hpf10.W,features = "shha",group.by = "cell_type")+NoLegend()
#
VlnPlot(hpf9,features = "shhb",group.by = "cell_type")+NoLegend()
VlnPlot(hpf8.W,features = "shhb",group.by = "cell_type")+NoLegend()
VlnPlot(hpf8.F,features = "shhb",group.by = "cell_type")+NoLegend()
VlnPlot(hpf6,features = "shhb",group.by = "cell_type")+NoLegend()
VlnPlot(hpf7,features = "shhb",group.by = "cell_type")+NoLegend()
# For Figure 2D left
hpf10$other <- "other"
hpf10@active.ident <- as.factor(hpf10$other)
DimPlot(hpf10,reduction = "tsne",cols = "#6FC8CB")+NoLegend()
p <- DimPlot(hpf10,label = F,repel = T,reduction = "tsne",cols = "#6FC8CB",pt.size = 0.00001)+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/OBO/hpf10.scRNAseq.tSNE.png",plot=p,device="png",dpi=300,units = "cm",width = 13,height = 13)
##
#
# for merge all scRNA-seq data for Figure S3A
# laoding all scRNA-seq data at each stage
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf3.8.RData")
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf4.3.RData")
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf4.8.RData")
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf5.3.RData")
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf6.RData")
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf7.RData")
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf8.RData")
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf9.RData")
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf10.RData")
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf11.RData")
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf12.RData")
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf14.RData")
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf18.RData")
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf24.RData")
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf48.RData")
# merge all sample
hpfall <- merge(hpf3.8,y=c(hpf4.3,hpf4.8,hpf5.3,hpf6,hpf7,
                           hpf8,hpf9,hpf10,hpf11,hpf12,
                           hpf14,hpf18,hpf24,hpf48))
#
hpfall

table(hpfall$Sample)
#
hpfall <- NormalizeData(hpfall, normalization.method = "LogNormalize", scale.factor = 10000)
hpfall <- FindVariableFeatures(hpfall, selection.method = "vst", nfeatures = 4000)
all.genes <- rownames(hpfall)
hpfall <- ScaleData(hpfall, features = all.genes)
hpfall <- RunPCA(hpfall, features = VariableFeatures(object = hpfall))
#
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = hpfall, reduction = "pca", pt.size = .1, group.by = "Sample")+NoLegend()
p2 <- VlnPlot(object = hpfall, features = "PC_1", group.by = "Sample", pt.size = .1)+NoLegend()
plot_grid(p1,p2)
#
options(repr.plot.height = 2.5, repr.plot.width = 6)
hpfall <- hpfall %>% 
  RunHarmony("Sample", plot_convergence = TRUE)
#
harmony_embeddings <- Embeddings(hpfall, 'harmony')
harmony_embeddings[1:5, 1:5]
#
options(repr.plot.height = 5, repr.plot.width = 12)
p1 <- DimPlot(object = hpfall, reduction = "harmony", pt.size = .1, group.by = "Sample")+NoLegend()
p2 <- VlnPlot(object = hpfall, features = "harmony_1", group.by = "Sample",pt.size = .1)+NoLegend()
plot_grid(p1,p2)
#
#
hpfall <- FindNeighbors(hpfall,reduction = "harmony", dims = 1:50)
hpfall <- RunUMAP(hpfall,reduction = "harmony", dims = 1:50,min.dist = 0.1)
hpfall <- RunTSNE(hpfall,reduction = "harmony", dims = 1:50)
#
DimPlot(hpfall)+NoLegend()
table(hpfall$group)
table(hpfall$stage)
#
hpfall$stage.group <- paste(hpfall$stage,hpfall$group,sep = "_")
table(hpfall$stage.group)
hpfall@active.ident <- as.factor(hpfall$stage.group)
levels(hpfall) <- c("3.8hpf_Farrell.etal","4.3hpf_Farrell.etal","4.8hpf_Farrell.etal","5.3hpf_Farrell.etal",
                    "6hpf_Farrell.etal","6hpf_Wagner.etal","7hpf_Farrell.etal",
                    "8hpf_Farrell.etal","8hpf_Wagner.etal","9hpf_Farrell.etal",
                    "10hpf_Farrell.etal","10hpf_Wagner.etal","11hpf_Farrell.etal",
                    "12hpf_Farrell.etal","14hpf_Wagner.etal",
                    "18hpf_Wagner.etal","24hpf_Farnsworth.etal","24hpf_Wagner.etal",
                    "48hpf_Farnsworth.etal")
#
New.ZtimeColor <- c("#FAD6A4","#F7B264","#E1520D","#AD3417",
                    "#F7C7D2","#EA66A1","#E50C84",
                    "#D2E2F2","#227DBD","#004695",
                    "#CFE6C0","#A8D392","#0C9E4A","#006D32",
                    "#C6C4E0","#857BB8","#694799",
                    "#EFBF00","#B78510","#935913")
DimPlot(hpfall,cols =New.ZtimeColor[1:19] )+NoLegend()
DimPlot(hpfall,group.by = "cell_type")+NoLegend()
DimPlot(hpfall,group.by = "cell_type",label = T,repel = T)+NoLegend()
#
p <- DimPlot(hpfall,label = F,repel = T,cols =New.ZtimeColor[1:19],pt.size = 0.00001)+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/OBO/hpf10/Fig_save/scRNAseq.all.UMAP.png",plot=p,device="png",dpi=300,units = "cm",width = 12,height = 12)
##
#
# For Figure S2A B
z24hpfmeta  <- readRDS("/home/sunkeyong/ZEPA/OBO/hpf24/V2/QualityControl/24hpf.SPATACseq.batch1/24hpf.SPATACseq.batch1-Pre-Filter-Metadata.rds")
z24hpfmeta <- as.data.frame(z24hpfmeta)
z24hpf.CG1 <- readRDS("/home/sunkeyong/ZEPA/Cell_genomics/bed/ArchR/QualityControl/24hpf_CG_rep1/24hpf_CG_rep1-Pre-Filter-Metadata.rds")
z24hpf.CG2 <- readRDS("/home/sunkeyong/ZEPA/Cell_genomics/bed/ArchR/QualityControl/24hpf_CG_rep2/24hpf_CG_rep2-Pre-Filter-Metadata.rds")
z24hpf.CG <- as.data.frame(rbind(z24hpf.CG1,z24hpf.CG2))
#
df
z24hpf.CG$log <- log10(z24hpf.CG[,2])
p <- ggPoint(
  x = z24hpf.CG[,12], 
  y = z24hpf.CG[,6], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(100), quantile(z24hpf.CG[,12], probs = 0.99999)),
  ylim = c(0, quantile(z24hpf.GC[,6], probs = 0.99999))
) + geom_hline(yintercept = 2, lty = "dashed") + geom_vline(xintercept = log10(1000), lty = "dashed")
p
#
z24hpfmeta$log <- log10(z24hpfmeta[,2])
p <- ggPoint(
  x = z24hpfmeta[,11], 
  y = z24hpfmeta[,6], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(100), quantile(z24hpfmeta[,11], probs = 0.99999)),
  ylim = c(0, quantile(z24hpfmeta[,6], probs = 0.99999))
) + geom_hline(yintercept = 8, lty = "dashed") + geom_vline(xintercept = log10(1500), lty = "dashed")
p
#
median(z24hpfmeta$nFrags)
median(z24hpfmeta$TSSEnrichment)
median(z24hpf.CG$nFrags)
median(z24hpf.CG$TSSEnrichment)
#
# for Figure S3B
hist(zhpf10$predictedScore_1,breaks = 40,col = "red")
#
# for Figure S3F
zhpf10.seurat$predictedGroup_1 <- zhpf10$predictedGroup_1
predictions <- table(zhpf10.seurat$ID, zhpf10.seurat$predictedGroup_1)
predictions2 <- predictions
predictions2[1:3,1:3]
dim(predictions2)
predictions2 <- predictions2[-24,]
predictions2 <- predictions2[-5,]
sum(diag(predictions2))/ncol(zhpf10.b1.seurat.0829)
#
predictions <- table(zhpf10.seurat$ID, zhpf10.seurat$predictedGroup_1)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
pheatmap(predictions,cluster_rows = F,cluster_cols = F,color = c("white","#4483B4"))
predictions <- as.data.frame(predictions)
p1 <- ggplot(predictions, aes(Var2, Var1, fill = Freq)) + geom_tile() + 
  scale_fill_gradient(name = "Fraction of cells",low = "white", high = "#4483B4")+ 
  xlab("Predicted cell type label (ATAC)")+ 
  ylab("Cell type annotation (RNA)") +
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p1
#
#
#
### 20230929 for Nature cell biology source data
colnames(predictions) <- c("Manual labeling","Predicted labeling","Fraction")
write.csv(predictions,file = "/home/sunkeyong/ZEPA/SourceData/FigS3F.csv")
#
FigureS3A <- as.data.frame(hpfall@reductions$umap@cell.embeddings)
FigureS3A$sample <- as.character(hpfall@active.ident)
head(FigureS3A)
write.csv(FigureS3A,file = "/home/sunkeyong/ZEPA/SourceData/FigS3A.csv")
#
write.csv(zhpf10$predictedScore_1,file = "/home/sunkeyong/ZEPA/SourceData/FigS3B.csv")
#
FigureS3G <- cbind(levels(hpf10.2.averages),cellratio.metadata.tmp2)
head(FigureS3G)
colnames(FigureS3G) <- c("Cell type","cell ratio in SPATAC-seq","cell ratio in scRNA-seq")
write.csv(FigureS3G,file = "/home/sunkeyong/ZEPA/SourceData/FigS3G.csv")
#
Fig2 <- zhpf10.seurat@meta.data
Fig2$UMAP1 <- zhpf10.seurat@reductions$pca@cell.embeddings[,1]
Fig2$UMAP2 <- zhpf10.seurat@reductions$pca@cell.embeddings[,2]
plot(Fig2$UMAP1,Fig2$UMAP2)
Fig2$IDF <- as.character(zhpf10.seurat@active.ident)
head(Fig2)
#
Fig2x <- as.data.frame(cbind(Fig2$Tile_LM1.C3,Fig2$Tile_LM1.C5,Fig2$IDF,Fig2$predictedGroup_1,Fig2$UMAP1,Fig2$UMAP2))
head(Fig2x)
colnames(Fig2x) <- c("resolution.C3","resolution.C5","ID","predicted labeling","tSNE1","tSNE2")
rownames(Fig2x) <- rownames(Fig2)
write.csv(Fig2x,file = "/home/sunkeyong/ZEPA/SourceData/Fig2_UMAP.csv")
#
write.csv(hpf10@reductions$tsne@cell.embeddings,file = "/home/sunkeyong/ZEPA/SourceData/Fig2D_scRNAseq.csv")
#
write.csv(Get.genescore.mat.matrix.suerat.averages@assays$RNA@scale.data[unique(top3$gene)[1:500],],
          file = "/home/sunkeyong/ZEPA/SourceData/Fig2I_scATACseq.csv")
write.csv(hpf10.2.averages@assays$RNA@scale.data[unique(top3$gene)[1:500],],
          file = "/home/sunkeyong/ZEPA/SourceData/Fig2I_scRNAseq.csv")
#
write.csv(cc1.t[23:44,1:22],file = "/home/sunkeyong/ZEPA/SourceData/Fig2H.csv")
#
#
#
#
#
# 
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################
# restart R 
# cell type annotation based on marker genes
# taking 72 hpf as an example 
setwd("/home/sunkeyong/ZEPA/OBO/hpf72")
#
###Zebrafish_
library(ArchR)
library(Seurat)
library("BSgenome.Drerio.UCSC.danRer11")
library(TxDb.Drerio.UCSC.danRer11.refGene)
library(org.Dr.eg.db)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)
library(pheatmap)
library(ComplexHeatmap)
library("RColorBrewer")
library("circlize")
#
# loading annotation files
load("/home/sunkeyong/ZEPA/GRCz11_anno/Lawson/storage/GRCz11.geneAnnotation.RData")
load("/home/sunkeyong/ZEPA/GRCz11_anno/Lawson/storage/genomeAnnotation.RData")
chr_rm <- read.table("/home/sunkeyong/ZEPA/GRCz11_anno/GRCz11_toplevel.chromosome.V2.txt",text = )
##
addArchRThreads(threads = 12) 
getArchRChrPrefix()
# load bed files
BedFiles <- c("/home/sunkeyong/ZEPA/All_arrow/72hpf.SPATACseq.batch1.bed.gz",
              "/home/sunkeyong/ZEPA/All_arrow/72hpf.SPATACseq.batch2.bed.gz",
              "/home/sunkeyong/ZEPA/All_arrow/72hpf.SPATACseq.batch4.bed.gz",
              "/home/sunkeyong/ZEPA/All_arrow/72hpf.10X.batch1.bed.gz")
names(BedFiles) = c("72hpf.SPATACseq.batch1",
                    "72hpf.SPATACseq.batch2",
                    "72hpf.SPATACseq.batch4",
                    "72hpf.10X.batch1")
# Create ArchR arrow file for downstream analysis
ArrowFiles <- createArrowFiles(
  inputFiles = BedFiles,
  sampleNames = names(BedFiles),
  minTSS = 4,
  minFrags = 500,
  offsetPlus = 0,
  offsetMinus = 0,
  geneAnnotation = GRCz11.geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  addTileMat = T,
  addGeneScoreMat = T,
  excludeChr =c("MT",paste("chr",chr_rm$V1,sep = "")))
# calculate Doublet score 
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)
#
zhpf72 <- ArchRProject(ArrowFiles = ArrowFiles,
                       geneAnnotation = GRCz11.geneAnnotation,
                       genomeAnnotation = genomeAnnotation,
                       copyArrows = T)
#
zhpf72
#
df <- getCellColData(zhpf72, select = c("log10(nFrags)", "TSSEnrichment"))
df
p <- ggPoint(
  x = df[,1], 
  y = df[,2], 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "Log10 Unique Fragments",
  ylabel = "TSS Enrichment",
  xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
  ylim = c(0, quantile(df[,2], probs = 0.99))
) + geom_hline(yintercept = 5, lty = "dashed") + geom_vline(xintercept = log10(1000), lty = "dashed")
p
#
zhpf72 <- addIterativeLSI(ArchRProj = zhpf72,
                          useMatrix = "TileMatrix", 
                          name = "Tile_IterativeLSI_M1", 
                          iterations = 3, 
                          clusterParams = list( #See Seurat::FindClusters
                            resolution = c(2,3), 
                            sampleCells = 10000, 
                            n.start = 10,
                            maxClusters = 100), 
                          varFeatures = 50000, 
                          dimsToUse = 1:100,
                          sampleCellsPre = 10000,
                          nPlot = 10000,
                          LSIMethod = 1,force = T) 
#
#
zhpf72 <- addHarmony(ArchRProj = zhpf72,
                     reducedDims = "Tile_IterativeLSI_M1",
                     name = "Tile_Harmony_M1",
                     groupBy = "Sample", #######
                     dimsToUse = 2:100,
                     force = T)
######
#
zhpf72 <- addUMAP(
  ArchRProj = zhpf72, 
  reducedDims = "Tile_Harmony_M1", 
  name = "Tile_Harmony_M1_UMAP1", 
  nNeighbors = 100, 
  dimsToUse = 1:100, 
  minDist = 0.1, 
  metric = "cosine",force = T
)
#
zhpf72 <- addUMAP(
  ArchRProj = zhpf72, 
  reducedDims = "Tile_Harmony_M1", 
  name = "Tile_Harmony_M1_UMAP2", 
  nNeighbors = 50, 
  dimsToUse = 1:100, 
  minDist = 0.1, 
  metric = "cosine",force = T
)
#
zhpf72 <- addTSNE(ArchRProj = zhpf72, 
                  reducedDims = "Tile_Harmony_M1", 
                  name = "Tile_Harmony_M1.TSNE1", 
                  perplexity = 200, dimsToUse = 1:100,
                  maxIterations = 1000,
                  force = T)
#
zhpf72 <- addTSNE(ArchRProj = zhpf72, 
                  reducedDims = "Tile_Harmony_M1", 
                  name = "Tile_Harmony_M1.TSNE2", 
                  perplexity = 100, dimsToUse = 1:100,
                  maxIterations = 1000,
                  force = T)
#
zhpf72 <- addClusters(input = zhpf72,name = "Tile_Harmony_M1.C1",
                      reducedDims = "Tile_Harmony_M1",
                      resolution = 1,method = "Seurat", 
                      dimsToUse = 1:100,maxClusters= 300,
                      force = T)
#
zhpf72 <- addClusters(input = zhpf72,name = "Tile_Harmony_M1.C3",
                      reducedDims = "Tile_Harmony_M1",
                      resolution = 3,method = "Seurat", 
                      dimsToUse = 1:100,maxClusters= 300,
                      force = T)
#
zhpf72 <- addClusters(input = zhpf72,name = "Tile_Harmony_M1.C5",
                      reducedDims = "Tile_Harmony_M1",
                      resolution = 5,method = "Seurat", 
                      dimsToUse = 1:100,maxClusters= 300,
                      force = T)
#
#
plotEmbedding(ArchRProj = zhpf72, colorBy = "cellColData", name = "Sample", embedding = "Tile_Harmony_M1.TSNE2")
plotEmbedding(ArchRProj = zhpf72, colorBy = "cellColData", name = "Sample", embedding = "Peak_Harmony_M1.TSNE2")
plotEmbedding(ArchRProj = zhpf72, colorBy = "cellColData", name = "Tile_Harmony_M1.C1", embedding = "Peak_Harmony_M1.TSNE2")
plotEmbedding(ArchRProj = zhpf72, colorBy = "cellColData", name = "Tile_Harmony_M1.C3", embedding = "Peak_Harmony_M1.TSNE2")
plotEmbedding(ArchRProj = zhpf72, colorBy = "cellColData", name = "Tile_Harmony_M1.C5", embedding = "Peak_Harmony_M1.TSNE2")
#
### Create one pseudo seurat file for quick processing
x.num <- c(rep(1,800000),rep(0,800000),rep(2,800000))
head(x.num)
x.num.sample <- cbind(sample(x.num,length(zhpf72$cellNames)),
                      sample(x.num,length(zhpf72$cellNames)),
                      sample(x.num,length(zhpf72$cellNames)),
                      sample(x.num,length(zhpf72$cellNames)),
                      sample(x.num,length(zhpf72$cellNames)))
head(x.num.sample)
colnames(x.num.sample) <- c("F1","F2","F3","F4","F5")
rownames(x.num.sample) <- zhpf72$cellNames
#
zhpf72.seurat <- CreateSeuratObject(counts = t(x.num.sample), assay = "ATAC", project = "zhpf72", min.cells = 0, min.features = 0)
zhpf72.seurat
zhpf72.seurat <- NormalizeData(zhpf72.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
zhpf72.seurat <- FindVariableFeatures(zhpf72.seurat, selection.method = "vst", nfeatures = 4)
zhpf72.seurat <- ScaleData(zhpf72.seurat, features = rownames(zhpf72.seurat))
variablegene <- VariableFeatures(object = zhpf72.seurat)
zhpf72.seurat <- RunPCA(zhpf72.seurat, features = variablegene,npcs =2,reduction.name = "tSNE",reduction.key = "tSNE_")
zhpf72.seurat
DimPlot(zhpf72.seurat)
#
zhpf72.seurat@reductions$pca@cell.embeddings[,1] <- zhpf72@embeddings$Tile_Harmony_M1.TSNE2$df[,1]
zhpf72.seurat@reductions$pca@cell.embeddings[,2] <- zhpf72@embeddings$Tile_Harmony_M1.TSNE2$df[,2]
#
DimPlot(zhpf72.seurat,label = T,repel = T,pt.size = 0.2)+NoLegend()
#
zhpf72.seurat$Tile_Harmony_M1.C1 <- zhpf72$Tile_Harmony_M1.C1
zhpf72.seurat$Tile_Harmony_M1.C3 <- zhpf72$Tile_Harmony_M1.C3
zhpf72.seurat$Tile_Harmony_M1.C5 <- zhpf72$Tile_Harmony_M1.C5
zhpf72.seurat$DoubletScore <- zhpf72$DoubletScore
zhpf72.seurat$DoubletEnrichment <- zhpf72$DoubletEnrichment
#
# imputating the gene score
zhpf72 <- addImputeWeights(zhpf72,reducedDims = "Tile_Harmony_M1")
# calculating the highly accessible genes
addArchRThreads(threads = 10) 
markersGS <- getMarkerFeatures(
  ArchRProj = zhpf72, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Tile_Harmony_M1.C3",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
#
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "glula", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "epcam", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "ly75", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "c1qa", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "cd63", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "hoxb3a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "rx1", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "gnat1", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "gngt2a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "six7", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "elavl3", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "vsx1", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "crx", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "cabp5b", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "tfap2a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "apoeb", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "ompa", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "rbpms2a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "pax6a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "rem1", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "lmo4a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "epcam", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "ihha", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "osr2", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "mkxa", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "lhx6", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "gsc", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "grem2b", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "dlx4a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "nkx3-1", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "tbx5a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "mymk", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "tbx1", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "myf5", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "apof", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "fpr1", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "gata5", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "gata6", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "tbx1", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "aplnr2", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "c7a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "tbx18", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "matn3b", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "and1", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "bgna", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "entpd5a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "phox2bb", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "eomesa", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "dlx1a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "npas4a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "kidins220a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "rx3", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "vsx2", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "hes2.2", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "foxg1b", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "shhb", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "shha", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "phox2bb", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "nkx6.2", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "nkx2.4a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "lrig1", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "efna2a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "fezf2", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "dlx4a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "slc13a4", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "itga11a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "entpd5a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "bgna", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "entpd5a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "slc9a3.1", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "chia.3", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "gck", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "aspn", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "hpdb", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "col9a1a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "ccdc3b", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "cyp1c1", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "slco2b1", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "otos", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "tnmd", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "epyc", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "sox19a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "vrk", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "sox19a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "sox19a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
#
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "Nr5a2", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "MotifMatrix", name = "z:nr5a2", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
#
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "Sox10", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix", name = "Tyr", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf72))
#
#
#
p <- plotBrowserTrack(ArchRProj = zhpf72, groupBy = "Tile_Harmony_M1.C3",
                      geneSymbol = "shha", upstream = 50000,downstream = 50000,
                      ylim = c(0.001,0.999999))
grid::grid.newpage()
grid::grid.draw(p$shha)
p <- plotBrowserTrack(ArchRProj = zhpf72, groupBy = "Tile_Harmony_M1.C3",
                      geneSymbol = "nr5a2", upstream = 150000,downstream = 20000,
                      ylim = c(0.001,0.999999))
grid::grid.newpage()
grid::grid.draw(p$nr5a2)
#
zhpf72.seurat@active.ident <- as.factor(zhpf72.seurat$Tile_Harmony_M1.C3)
#
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C1"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C2"="hpf72:Retina (RGC)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C3"="hpf72:Differentiating neurons (ganglion)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C4"="hpf72:Pronephric duct")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C5"="hpf72:Pharyngeal epidermis")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C6"="hpf72:Pharyngeal epidermis")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C7"="hpf72:Lateral line primordium")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C8"="hpf72:Otic placode")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C9"="hpf72:Periderm")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C10"="hpf72:Periderm")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C11"="hpf72:Gut?")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C12"="hpf72:Hatching gland")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C13"="hpf72:Hatching gland")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C14"="hpf72:Slow muscle cells")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C15"="hpf72:Slow muscle cells")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C16"="hpf72:Slow muscle cells")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C17"="hpf72:Slow muscle cells")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C18"="hpf72:Fast muscle cells")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C19"="hpf72:Fast muscle cells")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C20"="hpf72:Fast muscle cells")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C21"="hpf72:Fast muscle cells")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C22"="hpf72:Notochord")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C23"="hpf72:Notochord")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C24"="hpf72:Pigment cells")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C25"="hpf72:Pigment cells")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C26"="hpf72:Pigment cells")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C27"="hpf72:Pigment cells")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C28"="hpf72:Cranial ganglia")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C29"="hpf72:Cranial ganglia")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C30"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C31"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C32"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C33"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C34"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C35"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C36"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C37"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C38"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C39"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C40"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C41"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C42"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C43"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C44"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C45"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C46"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C47"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C48"="hpf72:Differentiating neurons (phox2bb+)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C49"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C50"="hpf72:Olfactory placode")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C51"="hpf72:Differentiating neurons (eomes+)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C52"="hpf72:Differentiating neurons (eomes+)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C53"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C54"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C55"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C56"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C57"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C58"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C59"="hpf72:Differentiating neurons")##Unknown (Interesting cluster)
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C60"="hpf72:Differentiating neurons")##Unknown (Interesting cluster)
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C61"="hpf72:Differentiating neurons (dlx1a+)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C62"="hpf72:Differentiating neurons (dlx1a+)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C63"="hpf72:Differentiating neurons (dlx1a+)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C64"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C65"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C66"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C67"="hpf72:Retina (Horizontal cells)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C68"="hpf72:Retina.Unknown")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C69"="hpf72:Retina (Cone bipolar cells)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C70"="hpf72:Retina (Cone bipolar cells)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C71"="hpf72:Retina (Cone bipolar cells)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C72"="hpf72:Retina (Cone bipolar cells)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C73"="hpf72:Retina (Cone bipolar cells)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C74"="hpf72:Retina (Cone bipolar cells)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C75"="hpf72:Retina (Amacrine)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C76"="hpf72:Retina (Amacrine)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C77"="hpf72:Retina (Amacrine)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C78"="hpf72:Retina (Amacrine)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C79"="hpf72:Retina (Muller glia)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C80"="hpf72:Retina neuroblasts")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C81"="hpf72:Retina neuroblasts")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C82"="hpf72:Midbrain")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C83"="hpf72:Midbrain")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C84"="hpf72:Diencephalon")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C85"="hpf72:Tailbud spinal cord")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C86"="hpf72:Tailbud spinal cord")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C87"="hpf72:Hindbrain dorsal")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C88"="hpf72:Hindbrain dorsal")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C89"="hpf72:Telencephalon")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C90"="hpf72:Midbrain ventral")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C91"="hpf72:Differentiating neurons")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C92"="hpf72:YSL")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C93"="hpf72:Endothelial")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C94"="hpf72:Cardiac muscle")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C95"="hpf72:Erythroid")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C96"="hpf72:Liver cells")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C97"="hpf72:Liver cells")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C98"="hpf72:Liver cells")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C99"="hpf72:Gut")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C100"="hpf72:Gut")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C101"="hpf72:Gut")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C102"="hpf72:Floorplate")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C103"="hpf72:Lens")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C104"="hpf72:Retina pigmented epithelium")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C105"="hpf72:Roofplate")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C106"="hpf72:Doublets")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C107"="hpf72:Doublets")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C108"="hpf72:Retina (Photoreceptor precursor cells)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C109"="hpf72:Retina (Photoreceptor precursor cells)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C110"="hpf72:Retina (Rods)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C111"="hpf72:Doublets")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C112"="hpf72:Retina (Cones)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C113"="hpf72:Retina (Cones)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C114"="hpf72:Retina (Cones)")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C115"="hpf72:Otic placode")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C116"="hpf72:Otic placode")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C117"="hpf72:Neural crest")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C118"="hpf72:Neural crest")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C119"="hpf72:Neural crest")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C120"="hpf72:Sclerotome?")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C121"="hpf72:Sclerotome?")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C122"="hpf72:Myoblast")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C123"="hpf72:Myoblast")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C124"="hpf72:Pectoral fin field")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C125"="hpf72:Heart field")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C126"="hpf72:Heart field")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C127"="hpf72:Heart field")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C128"="hpf72:Heart field")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C129"="hpf72:Heart field")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C130"="hpf72:Sclerotome")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C131"="hpf72:Myotome")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C132"="hpf72:Neural crest derived mesoderm")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C133"="hpf72:Neural crest derived mesoderm")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C134"="hpf72:Eye cornea")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C135"="hpf72:Pharyngeal arch")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C136"="hpf72:Pharyngeal arch")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C137"="hpf72:Neural crest derived mesoderm")
zhpf72.seurat <- RenameIdents(zhpf72.seurat,"C138"="hpf72:Neural crest derived mesoderm")
##
zhpf72.seurat$ID <- as.character(zhpf72.seurat@active.ident)
zhpf72$ID <- zhpf72.seurat$ID 
#
#
addArchRThreads(threads = 1) 
zhpf72 <- addGroupCoverages(ArchRProj = zhpf72, groupBy = "ID")
#
zhpf72 <- addReproduciblePeakSet(
  ArchRProj = zhpf72, 
  groupBy = "ID",
  genomeSize = 1.5e9,
  excludeChr = c("MT",paste("chr",chr_rm$V1,sep = "")), 
  genomeAnnotation = genomeAnnotation,
  geneAnnotation = GRCz11.geneAnnotation,
  pathToMacs2 = "/home/sunkeyong/.local/bin/macs2", force=T
)
#
zhpf72.peakset <- getPeakSet(zhpf72)
save(zhpf72.peakset,file="zhpf72.peakset.RData")
#
addArchRThreads(threads = 15) 
#
zhpf72 <- addPeakMatrix(ArchRProj = zhpf72, force =T )
##
zhpf72 <- addMotifAnnotations(ArchRProj = zhpf72, motifPWMs=danRer11.Motif, name = "Motif")
#
zhpf72 <- addBgdPeaks(zhpf72)
#
zhpf72 <- addDeviationsMatrix(
  ArchRProj = zhpf72, 
  peakAnnotation = "Motif",
  force = TRUE
)
#
getAvailableMatrices(zhpf72)
#
p <-plotEmbedding(ArchRProj = zhpf72, colorBy = "GeneScoreMatrix",size = 0.0001, name = "nr5a2",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf72))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/OBO/72hpf/Figure_save/nr5a2.GS.72hpf.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf72, colorBy = "MotifMatrix",size = 0.0001, name = "z:nr5a2",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf72))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/OBO/72hpf/Figure_save/nr5a2.motif.72hpf.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
#
#
DimPlot(zhpf72.seurat)+NoLegend()
zhpf72.seurat$Sample <- zhpf72$Sample
zhpf72.seurat$tech <- unlist(strsplit(as.character(zhpf72.seurat$Sample),".",fixed=T))[seq(2,3*ncol(zhpf72.seurat),3)]
zhpf72.seurat@active.ident <- as.factor(zhpf72.seurat$tech)
DimPlot(zhpf72.seurat,reduction = "tSNE",raster = F,shuffle = T)+NoLegend()
#
p <- DimPlot(zhpf72.seurat,label = F,repel = T,pt.size = 0.00001,reduction = "tSNE",raster = F)+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/OBO/72hpf/Figure_save/hpf72.tech.png",plot=p,device="png",dpi=300,units = "cm",width = 15,height = 15)
##
zhpf72.seurat@active.ident <- as.factor(zhpf72.seurat$ID)
DimPlot(zhpf72.seurat,reduction = "tSNE",raster = F,shuffle = T)+NoLegend()
#
p <- DimPlot(zhpf72.seurat,label = F,repel = T,pt.size = 0.00001,reduction = "tSNE",raster = F)+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/OBO/72hpf/Figure_save/hpf72.ID.png",plot=p,device="png",dpi=300,units = "cm",width = 15,height = 15)
##






