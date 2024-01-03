# Figure 5
# analyse the MOTIF
# taking the 24 hpf as an example
setwd("/home/sunkeyong/ZEPA/OBO/hpf24")
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
BedFiles <- c("/home/sunkeyong/ZEPA/All_arrow/24hpf.10X.batch1.bed.gz",
              "/home/sunkeyong/ZEPA/All_arrow/24hpf.10X.batch2.bed.gz",
              "/home/sunkeyong/ZEPA/All_arrow/24hpf.SPATACseq.batch1.bed.gz",
              "/home/sunkeyong/ZEPA/All_arrow/24hpf.SPATACseq.batch2.bed.gz",
              "/home/sunkeyong/ZEPA/All_arrow/24hpf.SPATACseq.batch3.bed.gz",
              "/home/sunkeyong/ZEPA/All_arrow/24hpf.SPATACseq.batch4.bed.gz",
              "/home/sunkeyong/ZEPA/All_arrow/24hpf.SPATACseq.batch5.bed.gz")
names(BedFiles) = c("24hpf.10X.batch1",
                    "24hpf.10X.batch2",
                    "24hpf.SPATACseq.batch1",
                    "24hpf.SPATACseq.batch2",
                    "24hpf.SPATACseq.batch3",
                    "24hpf.SPATACseq.batch4",
                    "24hpf.SPATACseq.batch5")
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
zhpf24 <- ArchRProject(ArrowFiles = ArrowFiles,
                       geneAnnotation = GRCz11.geneAnnotation,
                       genomeAnnotation = genomeAnnotation,
                       copyArrows = T)
#
zhpf24
#
df <- getCellColData(zhpf24, select = c("log10(nFrags)", "TSSEnrichment"))
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
zhpf24 <- addIterativeLSI(ArchRProj = zhpf24,
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
zhpf24 <- addHarmony(ArchRProj = zhpf24,
                     reducedDims = "Tile_IterativeLSI_M1",
                     name = "Tile_Harmony_M1",
                     groupBy = "Sample", #######
                     dimsToUse = 2:100,
                     force = T)
######
#
zhpf24 <- addUMAP(
  ArchRProj = zhpf24, 
  reducedDims = "Tile_Harmony_M1", 
  name = "Tile_Harmony_M1_UMAP1", 
  nNeighbors = 100, 
  dimsToUse = 1:100, 
  minDist = 0.1, 
  metric = "cosine",force = T
)
#
zhpf24 <- addUMAP(
  ArchRProj = zhpf24, 
  reducedDims = "Tile_Harmony_M1", 
  name = "Tile_Harmony_M1_UMAP2", 
  nNeighbors = 50, 
  dimsToUse = 1:100, 
  minDist = 0.1, 
  metric = "cosine",force = T
)
#
zhpf24 <- addTSNE(ArchRProj = zhpf24, 
                  reducedDims = "Tile_Harmony_M1", 
                  name = "Tile_Harmony_M1.TSNE1", 
                  perplexity = 200, dimsToUse = 1:100,
                  maxIterations = 1000,
                  force = T)
#
zhpf24 <- addTSNE(ArchRProj = zhpf24, 
                  reducedDims = "Tile_Harmony_M1", 
                  name = "Tile_Harmony_M1.TSNE2", 
                  perplexity = 100, dimsToUse = 1:100,
                  maxIterations = 1000,
                  force = T)
#
zhpf24 <- addClusters(input = zhpf24,name = "Tile_Harmony_M1.C1",
                      reducedDims = "Tile_Harmony_M1",
                      resolution = 1,method = "Seurat", 
                      dimsToUse = 1:100,maxClusters= 300,
                      force = T)
#
zhpf24 <- addClusters(input = zhpf24,name = "Tile_Harmony_M1.C3",
                      reducedDims = "Tile_Harmony_M1",
                      resolution = 3,method = "Seurat", 
                      dimsToUse = 1:100,maxClusters= 300,
                      force = T)
#
zhpf24 <- addClusters(input = zhpf24,name = "Tile_Harmony_M1.C5",
                      reducedDims = "Tile_Harmony_M1",
                      resolution = 5,method = "Seurat", 
                      dimsToUse = 1:100,maxClusters= 300,
                      force = T)
#
#
plotEmbedding(ArchRProj = zhpf24, colorBy = "cellColData", name = "Sample", embedding = "Tile_Harmony_M1.TSNE2")
plotEmbedding(ArchRProj = zhpf24, colorBy = "cellColData", name = "Sample", embedding = "Peak_Harmony_M1.TSNE2")
plotEmbedding(ArchRProj = zhpf24, colorBy = "cellColData", name = "Tile_Harmony_M1.C1", embedding = "Peak_Harmony_M1.TSNE2")
plotEmbedding(ArchRProj = zhpf24, colorBy = "cellColData", name = "Tile_Harmony_M1.C3", embedding = "Peak_Harmony_M1.TSNE2")
plotEmbedding(ArchRProj = zhpf24, colorBy = "cellColData", name = "Tile_Harmony_M1.C5", embedding = "Peak_Harmony_M1.TSNE2")
#
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/UMAP/hpf24.RData")
addArchRThreads(threads = 1) 
#
#
zhpf24 <- addGeneIntegrationMatrix(
  ArchRProj = zhpf24, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GIM",
  reducedDims = "Tile_Harmony_M1",
  seRNA = hpf24,
  addToArrow = T,
  force = TRUE,
  dimsToUse = 1:100,
  UMAPParams = list(n_neighbors = 50, min_dist = 0.1, metric = "cosine", verbose = FALSE),
  nGenes = 3000,
  sampleCellsATAC = 40000,
  sampleCellsRNA = 29526,
  groupRNA = "ID",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore")
#
### Create one pseudo seurat file for quick processing
x.num <- c(rep(1,800000),rep(0,800000),rep(2,800000))
head(x.num)
x.num.sample <- cbind(sample(x.num,length(zhpf24$cellNames)),
                      sample(x.num,length(zhpf24$cellNames)),
                      sample(x.num,length(zhpf24$cellNames)),
                      sample(x.num,length(zhpf24$cellNames)),
                      sample(x.num,length(zhpf24$cellNames)))
head(x.num.sample)
colnames(x.num.sample) <- c("F1","F2","F3","F4","F5")
rownames(x.num.sample) <- zhpf24$cellNames
#
zhpf24.seurat <- CreateSeuratObject(counts = t(x.num.sample), assay = "ATAC", project = "zhpf24", min.cells = 0, min.features = 0)
zhpf24.seurat
zhpf24.seurat <- NormalizeData(zhpf24.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
zhpf24.seurat <- FindVariableFeatures(zhpf24.seurat, selection.method = "vst", nfeatures = 4)
zhpf24.seurat <- ScaleData(zhpf24.seurat, features = rownames(zhpf24.seurat))
variablegene <- VariableFeatures(object = zhpf24.seurat)
zhpf24.seurat <- RunPCA(zhpf24.seurat, features = variablegene,npcs =2,reduction.name = "tSNE",reduction.key = "tSNE_")
zhpf24.seurat
DimPlot(zhpf24.seurat)
#
zhpf24.seurat@reductions$pca@cell.embeddings[,1] <- zhpf24@embeddings$Tile_Harmony_M1.TSNE2$df[,1]
zhpf24.seurat@reductions$pca@cell.embeddings[,2] <- zhpf24@embeddings$Tile_Harmony_M1.TSNE2$df[,2]
#
DimPlot(zhpf24.seurat,label = T,repel = T,pt.size = 0.2)+NoLegend()
#
zhpf24.seurat$Tile_Harmony_M1.C1 <- zhpf24$Tile_Harmony_M1.C1
zhpf24.seurat$Tile_Harmony_M1.C3 <- zhpf24$Tile_Harmony_M1.C3
zhpf24.seurat$Tile_Harmony_M1.C5 <- zhpf24$Tile_Harmony_M1.C5
zhpf24.seurat$DoubletScore <- zhpf24$DoubletScore
zhpf24.seurat$DoubletEnrichment <- zhpf24$DoubletEnrichment
zhpf24.seurat$predictedGroup <- zhpf24$predictedGroup
#
# imputating the gene score
zhpf24 <- addImputeWeights(zhpf24,reducedDims = "Tile_Harmony_M1")
# calculating the highly accessible genes
addArchRThreads(threads = 10) 
markersGS <- getMarkerFeatures(
  ArchRProj = zhpf24, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Tile_Harmony_M1.C3",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
#
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "phox2bb", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "hoxb3a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "fxr2", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "mylz3", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "myog", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "desma", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "myf5", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "six1b", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "stm", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "meox1", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "twist1a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "gata5", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "gata6", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "tbx15", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "tbx5a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "col2a1a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "pitx3", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "alx4a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "ppl", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "pmp22a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "foxc1b", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "isl1", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "elavl3", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "ndrg1a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "sox10", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "grem2b", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "grem2a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "mcamb", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "neurod1", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "msc", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "epd", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "myh6", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "hsd3b1", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "pkd1l2a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "pkd2l1", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "nppc", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "foxi3a", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix", name = "foxg1b", 
              embedding = "Tile_Harmony_M1.TSNE2",imputeWeights = getImputeWeights(zhpf24))
#
zhpf24.seurat@active.ident <- as.factor(zhpf24.seurat$Tile_Harmony_M1.C3)
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C1"="hpf24:Hatching gland")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C2"="hpf24:Hatching gland")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C3"="hpf24:Periderm")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C4"="hpf24:Periderm")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C5"="hpf24:Differentiating neurons")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C6"="hpf24:Differentiating neurons")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C7"="hpf24:Differentiating neurons")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C8"="hpf24:Differentiating neurons")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C9"="hpf24:Differentiating neurons")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C10"="hpf24:Differentiating neurons")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C11"="hpf24:Differentiating neurons")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C12"="hpf24:Differentiating neurons")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C13"="hpf24:Differentiating neurons")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C14"="hpf24:Differentiating neurons")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C15"="hpf24:Cranial ganglia")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C16"="hpf24:Cranial ganglia")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C17"="hpf24:Cardiac muscle")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C18"="hpf24:Slow muscle cells")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C19"="hpf24:Fast muscle cells")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C20"="hpf24:Fast muscle cells")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C21"="hpf24:Slow muscle cells")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C22"="hpf24:Fast muscle cells")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C23"="hpf24:Fast muscle cells")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C24"="hpf24:Notochord")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C25"="hpf24:Notochord")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C26"="hpf24:Notochord")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C27"="hpf24:Floorplate")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C28"="hpf24:Floorplate")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C29"="hpf24:Roofplate")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C30"="hpf24:Tailbud mesoderm")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C31"="hpf24:Tailbud mesoderm")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C32"="hpf24:Differentiating neurons (dlx1a+)")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C33"="hpf24:Differentiating neurons (eomesa+)")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C34"="hpf24:Diencephalon")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C35"="hpf24:Telencephalon")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C36"="hpf24:Telencephalon")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C37"="hpf24:Telencephalon")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C38"="hpf24:Midbrain")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C39"="hpf24:Midbrain")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C40"="hpf24:Midbrain ventral")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C41"="hpf24:Midbrain ventral")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C42"="hpf24:Telencephalon")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C43"="hpf24:Midbrain")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C44"="hpf24:Midbrain")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C45"="hpf24:Midbrain")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C46"="hpf24:Posterior ventral (phox2bb+)")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C47"="hpf24:Differentiating neurons")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C48"="hpf24:Tailbud spinal cord")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C49"="hpf24:Tailbud spinal cord")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C50"="hpf24:Tailbud spinal cord")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C51"="hpf24:Tailbud spinal cord")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C52"="hpf24:Differentiating neurons")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C53"="hpf24:Differentiating neurons")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C54"="hpf24:Differentiating neurons")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C55"="hpf24:Hindbrain dorsal")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C56"="hpf24:Hindbrain dorsal")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C57"="hpf24:Hindbrain dorsal")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C58"="hpf24:Hindbrain dorsal")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C59"="hpf24:Hindbrain dorsal")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C60"="hpf24:Hindbrain dorsal")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C61"="hpf24:Hindbrain dorsal")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C62"="hpf24:Hindbrain ventral")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C63"="hpf24:Hindbrain ventral")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C64"="hpf24:Hindbrain ventral")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C65"="hpf24:Hindbrain ventral")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C66"="hpf24:Hindbrain ventral")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C67"="hpf24:Doublets")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C68"="hpf24:Hindbrain ventral")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C69"="hpf24:Hindbrain ventral")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C70"="hpf24:Lens")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C71"="hpf24:Retina pigmented epithelium")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C72"="hpf24:Telencephalon")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C73"="hpf24:Retina neuroblasts")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C74"="hpf24:Optic cup")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C75"="hpf24:Optic cup")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C76"="hpf24:Optic cup")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C77"="hpf24:Retina?")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C78"="hpf24:Neural crest")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C79"="hpf24:Neural crest")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C80"="hpf24:Neural crest")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C81"="hpf24:Pharyngeal arch")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C82"="hpf24:Vessel progenitor")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C83"="hpf24:Cardiac neural crest")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C84"="hpf24:Cranial neural crest")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C85"="hpf24:Cranial neural crest")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C86"="hpf24:Pectoral fin field")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C87"="hpf24:Heart field")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C88"="hpf24:Heart field")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C89"="hpf24:Doublets")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C90"="hpf24:Doublets")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C91"="hpf24:Doublets")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C92"="hpf24:Myotome")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C93"="hpf24:Myotome")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C94"="hpf24:Sclerotome")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C95"="hpf24:Otic placode")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C96"="hpf24:Lateral line primordium")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C97"="hpf24:Epidermal (foxi3a+)")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C98"="hpf24:Olfactory placode")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C99"="hpf24:Doublets")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C100"="hpf24:Epidermal")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C101"="hpf24:Epidermal")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C102"="hpf24:Pharyngeal epidermis")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C103"="hpf24:Pharyngeal epidermis")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C104"="hpf24:YSL")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C105"="hpf24:Pronephric duct")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C106"="hpf24:Gut")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C107"="hpf24:Gut")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C108"="hpf24:Endothelial")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C109"="hpf24:Erythroid")
zhpf24.seurat <- RenameIdents(zhpf24.seurat,"C110"="hpf24:Immune cells")
#
zhpf24.seurat$ID <- as.character(zhpf24.seurat@active.ident)
zhpf24$ID <- zhpf24.seurat$ID 
############### 
# Peak calling
addArchRThreads(threads = 1) 
zhpf24 <- addGroupCoverages(ArchRProj = zhpf24, groupBy = "ID")
#
zhpf24 <- addReproduciblePeakSet(
  ArchRProj = zhpf24, 
  groupBy = "ID",
  genomeSize = 1.5e9,
  excludeChr = c("MT",paste("chr",chr_rm$V1,sep = "")), 
  genomeAnnotation = genomeAnnotation,
  geneAnnotation = GRCz11.geneAnnotation,
  pathToMacs2 = "/home/sunkeyong/.local/bin/macs2", force=T
)
#
zhpf24.peakset <- getPeakSet(zhpf24)
#
zhpf24 <- addPeakMatrix(ArchRProj = zhpf24, force =T )
# Motif annotation
# loading zebrafish motif collected from DANIO-CODE
load("/home/sunkeyong/ZEPA/Custom_danRer11_Motif/danRer11.Motif.RData")
#
zhpf24 <- addMotifAnnotations(ArchRProj = zhpf24, motifPWMs=danRer11.Motif, name = "Motif")
#
zhpf24 <- addBgdPeaks(zhpf24)
#
zhpf24 <- addDeviationsMatrix(
  ArchRProj = zhpf24, 
  peakAnnotation = "Motif",
  force = TRUE
)
#
getAvailableMatrices(zhpf24)
#
plotVarDev <- getVarDeviations(zhpf24, name = "MotifMatrix", plot = TRUE)
plotVarDev
#
#
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:ebf3a",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/ebf3a.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:etv2",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/etv2.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:tp63",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/tp63.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:foxa",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/foxa.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:flia",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/flia.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:mef2d",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/mef2d.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:zic1",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/zic1.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:grhl1",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/grhl1.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:gata5",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/gata5.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:hnf4a",quantCut = c(0.45, 0.99),
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/hnf4g.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:pax8",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/pax8.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:nr5a2",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/nr5a2.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:tfap2c",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/tfap2c.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:pou4f1",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/pou4f1.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:rx3",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/rx3.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:sox2",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/sox2.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:myog",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/myog.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:tbx5a",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/tbx5a.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:rfx3",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/rfx3.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:rx3",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/rx3.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:rx1",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/rx1.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:sox19a",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/sox19a.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:patz1",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/patz1.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:tfap4",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/tfap4.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:pou4f1",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/pou4f1.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:zeb1a",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/zeb1a.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <-plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:zeb1b",
                  embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_featureplot/zeb1b.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:zeb1b",
              embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix",size = 0.0001, name = "zeb1b",
              embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()
plotEmbedding(ArchRProj = zhpf24, colorBy = "GIM_Peak_LSI.M1.2",size = 0.0001, name = "zeb1b",
              embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()
#
plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:nr5a2",
              embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()
plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix",size = 0.0001, name = "nr5a2",
              embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()
plotEmbedding(ArchRProj = zhpf24, colorBy = "GIM_Peak_LSI.M1.2",size = 0.0001, name = "nr5a2",
              embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()
#
VlnPlot(hpf24.DC.Allon,features = "nr5a2",group.by = "cell_type")+NoLegend()
#
p <- plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:sox10",
                   embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_motif_featureplot/sox10.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix",size = 0.0001, name = "sox10",
                   embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_motif_featureplot/sox10.GS.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = zhpf24, colorBy = "GIM_Peak_LSI.M1.2",size = 0.0001, name = "sox10",continuousSet = "blueYellow",
                   embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_motif_featureplot/sox10.GIM.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
#
p <- plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:etv2",
                   embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_motif_featureplot/etv2.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix",size = 0.0001, name = "etv2",
                   embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_motif_featureplot/etv2.GS.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)

p <- plotEmbedding(ArchRProj = zhpf24, colorBy = "GIM_Peak_LSI.M1.2",size = 0.0001, name = "etv2",continuousSet = "blueYellow",
                   embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_motif_featureplot/etv2.GIM.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = zhpf24, colorBy = "MotifxMatrix",size = 0.0001, name = "z:nr5a2",
                   embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_motif_featureplot/nr5a2.motif.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = zhpf24, colorBy = "GeneScoreMatrix",size = 0.0001, name = "nr5a2",
                   embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_motif_featureplot/nr5a2.GS.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)

p <- plotEmbedding(ArchRProj = zhpf24, colorBy = "GIM_Peak_LSI.M1.2",size = 0.0001, name = "nr5a2",continuousSet = "blueYellow",
                   embedding = "Tile_Harmony_M1.TSNE2", imputeWeights = getImputeWeights(zhpf24))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/hpf24/motif_motif_featureplot/nr5a2.GIM.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
#
motifPositions <- getPositions(zhpf24,name="Motifx")
#
motifs <- c("sox7")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs
seFoot <- getFootprints(
  ArchRProj = zhpf24, 
  positions = motifPositions[markerMotifs], 
  groupBy = "ID")
#
plotFootprints(
  seFoot = seFoot,
  ArchRProj = zhpf24, 
  normMethod = "divide",
  plotName = "Footprints-Subtract-Bias.sox7",
  addDOC = F,
  pal =color.for.hpf24,
  smoothWindow = 5,
  height = 8,
  width = 6
)
#
# for Figure 6 flt1 and cdh5 track view
color.for.hpf24 <- c("#6A3D9A","#00BFFF",  "#B3E2CD", "#E6AB02", "#48D1CC", ##1-5
                     "#F4CAE4", "#EE82EE", "#00FFFF", "#99CCCC", "#ABD9E9", ##6-10
                     "#1A9850", "#9933CC", "#C2A5CF", "#87CEFA", "#E08214", ##11-15
                     "#9970AB", "#DE77AE", "#FDCDAC", "#3288BD", "#000000", ##16-20
                     "#66C2A5", "#FF00CC", "#F1E2CC", "#CC99FF", "#FDC086", ##21-25
                     "#9E0142", "#B3CDE3", "#D9D9D9", "#8DA0CB", "#999999", ##26-30
                     "#DFC27D", "#ABDDA4", "#FBB4AE", "#B35806", "#D6604D", ##31-35
                     "#8C510A", "#B15928", "#0099FF", "#B2ABD2", "#3288BD", ##36-40
                     "#276419", "#5E4FA2", "#D1E5F0", "#66BD63", "#FF33CC", ##41-45
                     "#B3EE3A")
#
p <- plotBrowserTrack(ArchRProj = zhpf24, groupBy = "ID", geneSymbol = "flt1", ylim =c(0.01,0.999),
                      upstream = 60000,downstream = 10000,pal =color.for.hpf24,
                      plotSummary = c("bulkTrack","featureTrack", "geneTrack"),
                      sizes = c(10, 1,1),)
grid::grid.newpage()
grid::grid.draw(p$flt1)
plotPDF(plotList = p, name = "track-flt1", ArchRProj = zhpf24,addDOC = FALSE, width = 5, height = 5)
#
p <- plotBrowserTrack(ArchRProj = zhpf24, groupBy = "ID", geneSymbol = "cdh5", ylim =c(0.01,1),
                      upstream = 17000,downstream = 35000,pal =color.for.hpf24,
                      plotSummary = c("bulkTrack","featureTrack", "geneTrack"),
                      sizes = c(10, 1,1),)
grid::grid.newpage()
grid::grid.draw(p$cdh5)
plotPDF(plotList = p, name = "track-cdh5", ArchRProj = zhpf24,addDOC = FALSE, width = 5, height = 5)
#
#
z24hpf.MotifMatrix <- getMatrixFromProject(
  ArchRProj = zhpf24,
  useMatrix = "MotifMatrix")
#
z24hpf.plotVarDev <- getVarDeviations(zhpf24, name = "MotifMatrix", plot = TRUE)
#
save(z24hpf.MotifMatrix,file="z24hpf.MotifMatrix.RData")
save(z24hpf.plotVarDev,file="z24hpf.plotVarDev.RData")
# save motif matrix 
# performed all above steps for all stages and save motif matrix
#
# resatrt R studio
# motif scRNA-scATAC-correlation
setwd("/home/sunkeyong/ZEPA/Motif")
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
#laoding all motif matrix
load("/home/sunkeyong/ZEPA/OBO/hpf5/z5hpf.MotifMatrix.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf6/z6hpf.MotifMatrix.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf7/z7hpf.MotifMatrix.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf8/z8hpf.MotifMatrix.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf9/z9hpf.MotifMatrix.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf10/z10hpf.MotifMatrix.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf11/z11hpf.MotifMatrix.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf12/z12hpf.MotifMatrix.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf14/z14hpf.MotifMatrix.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf18/z18hpf.MotifMatrix.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf24/z24hpf.MotifMatrix.RData")
#
rownames(z5hpf.MotifMatrix@assays@data$z) <- z5hpf.MotifMatrix@elementMetadata$name
rownames(z6hpf.MotifMatrix@assays@data$z) <- z6hpf.MotifMatrix@elementMetadata$name
rownames(z7hpf.MotifMatrix@assays@data$z) <- z7hpf.MotifMatrix@elementMetadata$name
rownames(z8hpf.MotifMatrix@assays@data$z) <- z8hpf.MotifMatrix@elementMetadata$name
rownames(z9hpf.MotifMatrix@assays@data$z) <- z9hpf.MotifMatrix@elementMetadata$name
rownames(z10hpf.MotifMatrix@assays@data$z) <- z10hpf.MotifMatrix@elementMetadata$name
rownames(z11hpf.MotifMatrix@assays@data$z) <- z11hpf.MotifMatrix@elementMetadata$name
rownames(z12hpf.MotifMatrix@assays@data$z) <- z12hpf.MotifMatrix@elementMetadata$name
rownames(z14hpf.MotifMatrix@assays@data$z) <- z14hpf.MotifMatrix@elementMetadata$name
rownames(z18hpf.MotifMatrix@assays@data$z) <- z18hpf.MotifMatrix@elementMetadata$name
rownames(z24hpf.MotifMatrix@assays@data$z) <- z24hpf.MotifMatrix@elementMetadata$name
#
z5hpf.Get.z.matrix <- z5hpf.MotifMatrix@assays@data$z
z6hpf.Get.z.matrix <- z6hpf.MotifMatrix@assays@data$z
z7hpf.Get.z.matrix <- z7hpf.MotifMatrix@assays@data$z
z8hpf.Get.z.matrix <- z8hpf.MotifMatrix@assays@data$z
z9hpf.Get.z.matrix <- z9hpf.MotifMatrix@assays@data$z
z10hpf.Get.z.matrix <- z10hpf.MotifMatrix@assays@data$z
z11hpf.Get.z.matrix <- z11hpf.MotifMatrix@assays@data$z
z12hpf.Get.z.matrix <- z12hpf.MotifMatrix@assays@data$z
z14hpf.Get.z.matrix <- z14hpf.MotifMatrix@assays@data$z
z18hpf.Get.z.matrix <- z18hpf.MotifMatrix@assays@data$z
z24hpf.Get.z.matrix <- z24hpf.MotifMatrix@assays@data$z
#
z24hpf.Get.z.matrix[1:3,1:3]
#
rownames(z24hpf.MotifMatrix@assays@data$z) <- z24hpf.MotifMatrix@elementMetadata$name
#
z5hpf.Get.z.matrix <- z5hpf.MotifMatrix@assays@data$z
z6hpf.Get.z.matrix <- z6hpf.MotifMatrix@assays@data$z
z7hpf.Get.z.matrix <- z7hpf.MotifMatrix@assays@data$z
z8hpf.Get.z.matrix <- z8hpf.MotifMatrix@assays@data$z
z9hpf.Get.z.matrix <- z9hpf.MotifMatrix@assays@data$z
z10hpf.Get.z.matrix <- z10hpf.MotifMatrix@assays@data$z
z11hpf.Get.z.matrix <- z11hpf.MotifMatrix@assays@data$z
z12hpf.Get.z.matrix <- z12hpf.MotifMatrix@assays@data$z
z14hpf.Get.z.matrix <- z14hpf.MotifMatrix@assays@data$z
z18hpf.Get.z.matrix <- z18hpf.MotifMatrix@assays@data$z
z24hpf.Get.z.matrix <- z24hpf.MotifMatrix@assays@data$z
#
z24hpf.Get.z.matrix[1:3,1:3]
#
ZEPA.motif <- CreateSeuratObject(counts = cbind(z5hpf.Get.z.matrix,
                                                z6hpf.Get.z.matrix,
                                                z7hpf.Get.z.matrix,
                                                z8hpf.Get.z.matrix,
                                                z9hpf.Get.z.matrix,
                                                z10hpf.Get.z.matrix,
                                                z11hpf.Get.z.matrix,
                                                z12hpf.Get.z.matrix,
                                                z14hpf.Get.z.matrix,
                                                z18hpf.Get.z.matrix,
                                                z24hpf.Get.z.matrix), 
                                 assay = "deviationsZscore", project = "zebrafish.deviations",
                                 min.cells = 0, min.features = 0)
#
ZEPA.motif
#
rm(z5hpf.MotifMatrix)
rm(z6hpf.MotifMatrix)
rm(z7hpf.MotifMatrix)
rm(z8hpf.MotifMatrix)
rm(z9hpf.MotifMatrix)
rm(z10hpf.MotifMatrix)
rm(z11hpf.MotifMatrix)
rm(z12hpf.MotifMatrix)
rm(z14hpf.MotifMatrix)
rm(z18hpf.MotifMatrix)
rm(z24hpf.MotifMatrix)
#
rm(z5hpf.Get.z.matrix)
rm(z6hpf.Get.z.matrix)
rm(z7hpf.Get.z.matrix)
rm(z8hpf.Get.z.matrix)
rm(z9hpf.Get.z.matrix)
rm(z10hpf.Get.z.matrix)
rm(z11hpf.Get.z.matrix)
rm(z12hpf.Get.z.matrix)
rm(z14hpf.Get.z.matrix)
rm(z18hpf.Get.z.matrix)
rm(z24hpf.Get.z.matrix)
#
Cells(ZEPA.motif)
#
load("/home/sunkeyong/ZEPA/All_cells/ZEPA.allsub.seurat.RData")
ZEPA.allsub.seura.meta <- as.data.frame(ZEPA.allsub.seurat@meta.data)
head(ZEPA.allsub.seura.meta)
ZEPA.allsub.seura.meta2 <- ZEPA.allsub.seura.meta[Cells(ZEPA.motif),]
dim(ZEPA.allsub.seura.meta2)
#
ZEPA.motif$ID <- ZEPA.allsub.seura.meta2$ID
ZEPA.motif$batch <- ZEPA.allsub.seura.meta2$batch
ZEPA.motif$time <- ZEPA.allsub.seura.meta2$time
#
View(as.data.frame(table(ZEPA.motif$ID)))
self.ID <- as.data.frame(table(ZEPA.motif$ID))
write.csv(self.ID,file = "self.ID.csv",row.names = F)
#
Co.ID <- read.csv("CoID.csv")
ZEPA.motif@active.ident <- as.factor(ZEPA.motif$ID)
ZEPA.motif.sub <- subset(ZEPA.motif,idents = unique(Co.ID$SelfID))
ZEPA.motif.sub
#
length(table(ZEPA.motif.sub@active.ident))
length(table(ZEPA.motif@active.ident))
#
Co.ID.Motif <- Co.ID[!duplicated(Co.ID[,"SelfID"]),]
head(Co.ID.Motif)
dim(Co.ID.Motif)
table(Co.ID.Motif$SelfID==Co.ID.Motif$CoID)
length(unique(Co.ID$CoID))
#
ZEPA.motif.sub.meta <- as.data.frame(ZEPA.motif.sub@meta.data)
ZEPA.motif.sub.meta$SelfID <- ZEPA.motif.sub.meta$ID
ZEPA.motif.sub.meta$cells <- rownames(ZEPA.motif.sub.meta)
ZEPA.motif.sub.meta.x <- merge(ZEPA.motif.sub.meta,Co.ID.Motif,by="SelfID")
head(ZEPA.motif.sub.meta.x)
rownames(ZEPA.motif.sub.meta.x) <- ZEPA.motif.sub.meta.x$cells
ZEPA.motif.sub.meta.x2 <- ZEPA.motif.sub.meta.x[rownames(ZEPA.motif.sub@meta.data),]
ZEPA.motif.sub$CoID <- ZEPA.motif.sub.meta.x2$CoID
ZEPA.motif.sub@active.ident <- as.factor(ZEPA.motif.sub$CoID)
length(table(ZEPA.motif.sub@active.ident))
#
#
ZEPA.motif.sub@active.ident <- as.factor(ZEPA.motif.sub$CoID)
ZEPA.motif.sub <- ScaleData(ZEPA.motif.sub, features = rownames(ZEPA.motif.sub))
ZEPA.motif.sub.avg <- AverageExpression(ZEPA.motif.sub,return.seurat = T,slot = "counts")
ZEPA.motif.sub.avg
#
DoHeatmap(ZEPA.motif.sub.avg,group.by = "Co.ID", features = "myog",label = T)+NoLegend()
#
#
load("/home/sunkeyong/ZEPA/Zebrafish_TOME/reprocess/allmerge/hpfallmerge.RData")
#
hpfall
table(hpfall$stage)
DimPlot(hpfall,reduction = "harmony")+NoLegend()
DimPlot(hpfall,reduction = "umap",group.by = "stage")+NoLegend()
DimPlot(hpfall,reduction = "umap1",group.by = "stage")+NoLegend()
DimPlot(hpfall,reduction = "umap2",group.by = "stage")+NoLegend()
DimPlot(hpfall,reduction = "umap2",group.by = "cell_type2")+NoLegend()
DimPlot(hpfall,reduction = "tsne",group.by = "stage")+NoLegend()
DimPlot(hpfall,reduction = "tsne1",group.by = "stage")+NoLegend()
#
table(hpfall$cell_type)
table(hpfall$group)
table(hpfall$cell_type2)
table(hpfall$stage)
#
hpfall$cell_typex <- paste(hpfall$stage,hpfall$cell_type,sep = "_")
table(hpfall$cell_typex,hpfall$stage)
table(hpfall$group,hpfall$stage)
View(as.data.frame(table(hpfall$cell_typex)))
View(table(hpfall$cell_typex,hpfall$group))
#
SD.ID <- as.data.frame(table(hpfall$cell_typex))
write.csv(SD.ID,file = "SD.ID.csv",row.names = F)
#
length(table(hpfall@active.ident))
hpfall@active.ident <- as.factor(hpfall$cell_typex)
hpfall <- subset(hpfall,idents = unique(Co.ID$Var1))
length(table(hpfall@active.ident))
#
Co.ID.RNA <- Co.ID[!duplicated(Co.ID[,"Var1"]),]
hpfall.meta <- as.data.frame(hpfall@meta.data)
hpfall.meta$Var1 <- hpfall.meta$cell_typex
hpfall.meta$cells <- rownames(hpfall.meta)
hpfall.meta.x <- merge(hpfall.meta,Co.ID.RNA,by="Var1")
dim(hpfall.meta.x)
rownames(hpfall.meta.x) <- hpfall.meta.x$cells
hpfall.meta.x2 <- hpfall.meta.x[rownames(hpfall@meta.data),]
#
hpfall$CoID <- hpfall.meta.x2$CoID
hpfall@active.ident <- as.factor(hpfall$CoID)
length(table(hpfall@active.ident))
#
hpfall.avg <- AverageExpression(hpfall,return.seurat = T,slot = "data")
hpfall.avg
ZEPA.motif.sub.avg
#
hpfall.avg$Co.ID <- as.character(hpfall.avg@active.ident)
ZEPA.motif.sub.avg$Co.ID <- as.character(ZEPA.motif.sub.avg@active.ident)
#
hpfall.avg$time <- unlist(strsplit(as.character(hpfall.avg@active.ident),":"))[seq(1,2*200,2)]
ZEPA.motif.sub.avg$time <- unlist(strsplit(as.character(ZEPA.motif.sub.avg@active.ident),":"))[seq(1,2*200,2)]
#
hpfall.avg@active.ident <- as.factor(hpfall.avg$time)
ZEPA.motif.sub.avg@active.ident <- as.factor(ZEPA.motif.sub.avg$time)
table(hpfall.avg@active.ident)
#
ZEPA.motif.sub.avg <- ScaleData(ZEPA.motif.sub.avg, features = rownames(ZEPA.motif.sub.avg))
#
co.TFnames <- intersect(rownames(hpfall.avg),rownames(ZEPA.motif.sub.avg))
#
x1 <- subset(hpfall.avg,idents = "hpf24")@assays$RNA@data[co.TFnames,]
y1 <- subset(ZEPA.motif.sub.avg,idents = "hpf24")@assays$deviationsZscore@counts[co.TFnames,]
#
DoHeatmap(subset(ZEPA.motif.sub.avg,idents="hpf24"),group.by = "Co.ID", features = "sox19a",slot = "counts")
DoHeatmap(subset(ZEPA.motif.sub.avg,idents="hpf24"),group.by = "Co.ID", features = "sox19a",label = T)
DoHeatmap(subset(ZEPA.motif.sub.avg,idents="hpf24"),group.by = "Co.ID", features = "myog",label = T,slot = "counts")
DoHeatmap(subset(ZEPA.motif.sub.avg,idents="hpf24"),group.by = "Co.ID", features = "myog",label = T)
#
dim(x1)
dim(y1)
x1[1:3,1:3]
#
colnames(x1)==colnames(y1)
#
timepoint <- c("hpf5","hpf6","hpf7","hpf8","hpf9",
               "hpf10","hpf11","hpf12","hpf14","hpf18","hpf24")
#
cor.spearman.pvalue <- c()
cor.spearman.rho <- c()
cor.pearson.pvalue <- c()
cor.pearson.cor <- c()
for (j in 1:11) {
  cor.tmp <- c()
  cor.tmp.spearman.pvalue <- c()
  cor.tmp.spearman.rho <- c()
  cor.tmp.pearson.pvalue <- c()
  cor.tmp.pearson.cor <- c()
  x1 <- subset(hpfall.avg,idents = timepoint[j])@assays$RNA@data[co.TFnames,]
  y1 <- subset(ZEPA.motif.sub.avg,idents = timepoint[j])@assays$deviationsZscore@counts[co.TFnames,]
  #
  for (i in 1:836) {
    cor.tmp <- cor.test(x1[i,],y1[i,],method = "spearman",exact=FALSE)
    cor.tmp.spearman.rho <- c(cor.tmp.spearman.rho,cor.tmp$estimate)
    cor.tmp.spearman.pvalue <- c(cor.tmp.spearman.pvalue,cor.tmp$p.value)
    cor.tmp <- cor.test(x1[i,],y1[i,],exact=FALSE)
    cor.tmp.pearson.pvalue <- c(cor.tmp.pearson.pvalue,cor.tmp$p.value)
    cor.tmp.pearson.cor <- c(cor.tmp.pearson.cor,cor.tmp$estimate)
  }
  cor.spearman.pvalue <- c(cor.spearman.pvalue,cor.tmp.spearman.pvalue)
  cor.spearman.rho <- c(cor.spearman.rho,cor.tmp.spearman.rho)
  cor.pearson.pvalue <- c(cor.pearson.pvalue,cor.tmp.pearson.pvalue)
  cor.pearson.cor <- c(cor.pearson.cor,cor.tmp.pearson.cor)
}
#
plot(x1[i,],y1[i,])
cor.test(x1[i,],y1[i,],method = "spearman",exact=FALSE)
cor.test(x1[i,],y1[i,],exact=FALSE)
#
cor.spearman.rho[1]
cor.spearman.rho[836*1+1]
cor.spearman.rho[836*2+1]
cor.spearman.rho[836*4+1]
#
cor.pearson.cor[1]
cor.pearson.cor[836*1+1]
cor.pearson.cor[836*2+1]
cor.pearson.cor[836*4+1]
#
plot(cor.spearman.rho,cor.pearson.cor)
#
cor.final <- cbind(cor.spearman.rho,cor.spearman.pvalue,cor.pearson.cor,cor.pearson.pvalue)
cor.final <- as.data.frame(cor.final)
cor.final$time <- rep(timepoint,each=836)
cor.final$motif <- rep(co.TFnames,times=11)
#
head(cor.final)
#
load("/home/sunkeyong/ZEPA/OBO/hpf5/z5hpf.plotVarDev.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf6/z6hpf.plotVarDev.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf7/z7hpf.plotVarDev.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf8/z8hpf.plotVarDev.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf9/z9hpf.plotVarDev.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf10/z10hpf.plotVarDev.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf11/z11hpf.plotVarDev.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf12/z12hpf.plotVarDev.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf14/z14hpf.plotVarDev.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf18/z18hpf.plotVarDev.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf24/z24hpf.plotVarDev.RData")
#
z5hpf.plotVarDev.V <- as.data.frame(z5hpf.plotVarDev$data)
z6hpf.plotVarDev.V <- as.data.frame(z6hpf.plotVarDev$data)
z7hpf.plotVarDev.V <- as.data.frame(z7hpf.plotVarDev$data)
z8hpf.plotVarDev.V <- as.data.frame(z8hpf.plotVarDev$data)
z9hpf.plotVarDev.V <- as.data.frame(z9hpf.plotVarDev$data)
z10hpf.plotVarDev.V <- as.data.frame(z10hpf.plotVarDev$data)
z11hpf.plotVarDev.V <- as.data.frame(z11hpf.plotVarDev$data)
z12hpf.plotVarDev.V <- as.data.frame(z12hpf.plotVarDev$data)
z14hpf.plotVarDev.V <- as.data.frame(z14hpf.plotVarDev$data)
z18hpf.plotVarDev.V <- as.data.frame(z18hpf.plotVarDev$data)
z24hpf.plotVarDev.V <- as.data.frame(z24hpf.plotVarDev$data)
#
rownames(z5hpf.plotVarDev.V) <- z5hpf.plotVarDev.V$name
rownames(z6hpf.plotVarDev.V) <- z6hpf.plotVarDev.V$name
rownames(z7hpf.plotVarDev.V) <- z7hpf.plotVarDev.V$name
rownames(z8hpf.plotVarDev.V) <- z8hpf.plotVarDev.V$name
rownames(z9hpf.plotVarDev.V) <- z9hpf.plotVarDev.V$name
rownames(z10hpf.plotVarDev.V) <- z10hpf.plotVarDev.V$name
rownames(z11hpf.plotVarDev.V) <- z11hpf.plotVarDev.V$name
rownames(z12hpf.plotVarDev.V) <- z12hpf.plotVarDev.V$name
rownames(z14hpf.plotVarDev.V) <- z14hpf.plotVarDev.V$name
rownames(z18hpf.plotVarDev.V) <- z18hpf.plotVarDev.V$name
rownames(z24hpf.plotVarDev.V) <- z24hpf.plotVarDev.V$name
#
z5hpf.plotVarDev.V <- z5hpf.plotVarDev.V[co.TFnames,]
z6hpf.plotVarDev.V <- z6hpf.plotVarDev.V[co.TFnames,]
z7hpf.plotVarDev.V <- z7hpf.plotVarDev.V[co.TFnames,]
z8hpf.plotVarDev.V <- z8hpf.plotVarDev.V[co.TFnames,]
z9hpf.plotVarDev.V <- z9hpf.plotVarDev.V[co.TFnames,]
z10hpf.plotVarDev.V <- z10hpf.plotVarDev.V[co.TFnames,]
z11hpf.plotVarDev.V <- z11hpf.plotVarDev.V[co.TFnames,]
z12hpf.plotVarDev.V <- z12hpf.plotVarDev.V[co.TFnames,]
z14hpf.plotVarDev.V <- z14hpf.plotVarDev.V[co.TFnames,]
z18hpf.plotVarDev.V <- z18hpf.plotVarDev.V[co.TFnames,]
z24hpf.plotVarDev.V <- z24hpf.plotVarDev.V[co.TFnames,]
#
rm(z5hpf.plotVarDev.V);rm(z5hpf.plotVarDev)
rm(z6hpf.plotVarDev.V);rm(z6hpf.plotVarDev)
rm(z7hpf.plotVarDev.V);rm(z7hpf.plotVarDev)
rm(z8hpf.plotVarDev.V);rm(z8hpf.plotVarDev)
rm(z9hpf.plotVarDev.V);rm(z9hpf.plotVarDev)
rm(z10hpf.plotVarDev.V);rm(z10hpf.plotVarDev)
rm(z11hpf.plotVarDev.V);rm(z11hpf.plotVarDev)
rm(z12hpf.plotVarDev.V);rm(z12hpf.plotVarDev)
rm(z14hpf.plotVarDev.V);rm(z14hpf.plotVarDev)
rm(z18hpf.plotVarDev.V);rm(z18hpf.plotVarDev)
rm(z24hpf.plotVarDev.V);rm(z24hpf.plotVarDev)
gc()
#
#
combinedVars <- c(z5hpf.plotVarDev.V$combinedVars,
                  z6hpf.plotVarDev.V$combinedVars,
                  z7hpf.plotVarDev.V$combinedVars,
                  z8hpf.plotVarDev.V$combinedVars,
                  z9hpf.plotVarDev.V$combinedVars,
                  z10hpf.plotVarDev.V$combinedVars,
                  z11hpf.plotVarDev.V$combinedVars,
                  z12hpf.plotVarDev.V$combinedVars,
                  z14hpf.plotVarDev.V$combinedVars,
                  z18hpf.plotVarDev.V$combinedVars,
                  z24hpf.plotVarDev.V$combinedVars)
#
combinedMeans <- c(z5hpf.plotVarDev.V$combinedMeans,
                   z6hpf.plotVarDev.V$combinedMeans,
                   z7hpf.plotVarDev.V$combinedMeans,
                   z8hpf.plotVarDev.V$combinedMeans,
                   z9hpf.plotVarDev.V$combinedMeans,
                   z10hpf.plotVarDev.V$combinedMeans,
                   z11hpf.plotVarDev.V$combinedMeans,
                   z12hpf.plotVarDev.V$combinedMeans,
                   z14hpf.plotVarDev.V$combinedMeans,
                   z18hpf.plotVarDev.V$combinedMeans,
                   z24hpf.plotVarDev.V$combinedMeans)
#
cor.final$combinedVars <- combinedVars
cor.final$combinedMeans <- combinedMeans
#
head(cor.final)
#
xyin <- cor.final[which(cor.final$time=="hpf24"),]
head(xyin)
plot(xyin$cor.spearman.rho,xyin$combinedVars)
plot(xyin$cor.pearson.cor,xyin$combinedVars)
#
plot(xyin$cor.spearman.rho,log2(xyin$combinedVars+1))
plot(xyin$cor.pearson.cor,log2(xyin$combinedVars+1))
#
boxplot(cor.final$cor.spearman.rho)
boxplot(cor.final$cor.pearson.cor)
#
x1 <- subset(hpfall.avg,idents = "hpf24")@assays$RNA@data[co.TFnames,]
y1 <- subset(ZEPA.motif.sub.avg,idents = "hpf24")@assays$deviationsZscore@counts[co.TFnames,]
#scale.data | counts | data
head(x1)
x.taget <- "myog"
cor.test(x1[x.taget,],y1[x.taget,],method = "spearman",exact=FALSE)
plot(x1[x.taget,],y1[x.taget,])
#
cor.final.24hpf <- cor.final[which(cor.final$time=="hpf24"),]
dim(cor.final.24hpf)
plot(cor.final.24hpf$cor.spearman.rho,cor.final.24hpf$combinedVars)
plot(cor.final.24hpf$cor.pearson.cor,cor.final.24hpf$combinedVars)
#
cor.final.24hpf$TFRegulator <- c("NO")
head(cor.final.24hpf)
#
cor.final.24hpf$TFRegulator[which(cor.final.24hpf$cor.pearson.cor > 0.5 & cor.final.24hpf$cor.pearson.pvalue < 0.01 & cor.final.24hpf$combinedVars > quantile(cor.final.24hpf$combinedVars, 0.75))] <- "YES1"
cor.final.24hpf$TFRegulator[which(cor.final.24hpf$cor.pearson.cor < -0.5 & cor.final.24hpf$cor.pearson.pvalue < 0.01 & cor.final.24hpf$combinedVars > quantile(cor.final.24hpf$combinedVars, 0.75))] <- "YES2"
#
cor.final.24hpf$TFRegulator[which(cor.final.24hpf$cor.pearson.cor > 0.5 & cor.final.24hpf$cor.pearson.pvalue < 0.01 & cor.final.24hpf$combinedVars > 1)] <- "YES1"
cor.final.24hpf$TFRegulator[which(cor.final.24hpf$cor.pearson.cor < -0.5 & cor.final.24hpf$cor.pearson.pvalue < 0.01 & cor.final.24hpf$combinedVars > 1)] <- "YES2"
#
cor.final.24hpf$TFRegulator[which(cor.final.24hpf$cor.pearson.cor > 0.5 & cor.final.24hpf$cor.pearson.pvalue < 0.01)] <- "YES"
cor.final.24hpf$TFRegulator[which(cor.final.24hpf$cor.pearson.cor < -0.5 & cor.final.24hpf$cor.pearson.pvalue < 0.01)] <- "NO"
#
cor.final.24hpf$TFRegulator[which(cor.final.24hpf$cor.pearson.cor > 0.5)] <- "YES"
cor.final.24hpf$TFRegulator[which(cor.final.24hpf$cor.pearson.cor < -0.5)] <- "NO"
#
table(cor.final.24hpf$TFRegulator)
#cor.final.24hpf$TFRegulator[which(cor.final.24hpf$cor.pearson.cor < -0.5)] <- "YES"
p <- ggplot(data.frame(cor.final.24hpf), aes(cor.pearson.cor, combinedVars, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="#7F7F7F", "YES1"="#CD2626", "YES2"="#009944")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(cor.final.24hpf$combinedVars)*1.05)
  )
p
#
write.table(cor.final,file = "ZEPA.allmotif.cor.final.csv",quote=F, sep = "\t",row.names = F,col.names = T)
#


