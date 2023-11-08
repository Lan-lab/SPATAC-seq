setwd("/home/sunkeyong/MOPA_project")  
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
library("RColorBrewer")
library("circlize")
#
addArchRGenome("mm10")
addArchRThreads(threads = 56) 
inputFiles <- c("/home/sunkeyong/MOPA_project/fragment_total/E1A.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E1B.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E1C.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E1D.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E1E.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E1F.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E1G.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E1H.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E2A.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E2B.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E2C.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E2D.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E2E.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E2F.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E2G.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E2H.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E7A.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E7B.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E7C.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E7D.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E7E.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E7F.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E7G.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E7H.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E8A.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E8B.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E8C.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E8D.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E8E.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E8F.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E8G.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E8H.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E11A.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E11B.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E11C.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E11D.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E11E.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E11F.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E11G.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E11H.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E17A.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E17B.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E17C.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E17D.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E17E.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E17F.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E17G.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E17H.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E18A.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E18B.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E18C.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E18D.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E18E.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E18F.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E18G.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E18H.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E19A.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E19B.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E19C.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E19D.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E19E.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E19F.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E19G.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E19H.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E153.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E154.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E155.bed.gz",
                "/home/sunkeyong/MOPA_project/fragment_total/E156.bed.gz")
#
names(inputFiles) <- c("E1A","E1B","E1C","E1D","E1E","E1F","E1G","E1H",
                       "E2A","E2B","E2C","E2D","E2E","E2F","E2G","E2H",
                       "E7A","E7B","E7C","E7D","E7E","E7F","E7G","E7H",
                       "E8A","E8B","E8C","E8D","E8E","E8F","E8G","E8H",
                       "E11A","E11B","E11C","E11D","E11E","E11F","E11G","E11H",
                       "E17A","E17B","E17C","E17D","E17E","E17F","E17G","E17H",
                       "E18A","E18B","E18C","E18D","E18E","E18F","E18G","E18H",
                       "E19A","E19B","E19C","E19D","E19E","E19F","E19G","E19H",
                       "E153","E154","E155","E156")
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, 
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)
#
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, 
  knnMethod = "UMAP", 
  LSIMethod = 1,force = T
)
#
MOPA <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  copyArrows = TRUE)
##
idxPass <- which(MOPA$DoubletEnrichment <= 4)
cellsPass <- MOPA$cellNames[idxPass]
MOPA <- MOPA[cellsPass, ]
#
MOPA <- addIterativeLSI(
  ArchRProj = MOPA,
  useMatrix = "TileMatrix", 
  name = "TileMatrix_IterativeLSI", 
  iterations = 3, 
  clusterParams = list( 
    resolution = c(2,3), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 50000, 
  dimsToUse = 1:100
)
MOPA <- addClusters(
  input = MOPA,
  reducedDims = "TileMatrix_IterativeLSI",
  method = "Seurat",
  name = "TileMatrix_IterativeLSI.Cluster.R4",
  resolution = 4
)
MOPA <- addGroupCoverages(ArchRProj = MOPA, groupBy = "TileMatrix_IterativeLSI.Cluster.R4")
MOPA <- addReproduciblePeakSet(
  ArchRProj = MOPA, 
  groupBy = "TileMatrix_IterativeLSI.Cluster.R4",
  pathToMacs2 = "/home/sunkeyong/.local/bin/macs2", force=T
)
#
MOPA <- addPeakMatrix(ArchRProj = MOPA, force =T)
MOPA <- addIterativeLSI(
  ArchRProj = MOPA,
  useMatrix = "PeakMatrix", 
  name = "PeakMatrix_IterativeLSI", 
  iterations = 3, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(3,4), 
    sampleCells = 10000, 
    n.start = 10,
    maxClusters = 100
  ), 
  varFeatures = 100000, 
  dimsToUse = 1:100,
  force = T
)
MOPA <- addClusters(
  input = MOPA,
  reducedDims = "PeakMatrix_IterativeLSI",
  method = "Seurat",
  name = "PeakMatrix_IterativeLSI.Cluster.R4",
  resolution = 3
)
MOPA <- addTSNE(
  ArchRProj = MOPA, 
  reducedDims = "PeakMatrix_IterativeLSI", 
  name = "PeakMatrix_IterativeLSI_tSNE", 
  perplexity = 100, 
  dimsToUse = 1:100,
  maxIterations = 1000,
  force = T
)
MOPA <- addImputeWeights(MOPA,reducedDims = "PeakMatrix_IterativeLSI")
#
plotEmbedding(
  ArchRProj = MOPA, 
  colorBy = "GeneScoreMatrix", 
  name = "Myod1", 
  embedding = "PeakMatrix_IterativeLSI_tSNE",
  imputeWeights = getImputeWeights(MOPA)
)
#
# Create pseudo-Seurat object for rapid processing metadata
#
x.num <- c(rep(1,800000),rep(0,800000),rep(2,800000))
head(x.num)
x.num.sample <- cbind(sample(x.num,length(MOPA$cellNames)),
                      sample(x.num,length(MOPA$cellNames)),
                      sample(x.num,length(MOPA$cellNames)),
                      sample(x.num,length(MOPA$cellNames)),
                      sample(x.num,length(MOPA$cellNames)))
head(x.num.sample)
colnames(x.num.sample) <- c("F1","F2","F3","F4","F5")
rownames(x.num.sample) <- z.Notochordx$cellNames
#
MOPA.seurat_type <- CreateSeuratObject(counts = t(x.num.sample), assay = "ATAC", project = "MOPA", min.cells = 0, min.features = 0)
MOPA.seurat_type
MOPA.seurat_type <- NormalizeData(MOPA.seurat_type, normalization.method = "LogNormalize", scale.factor = 10000)
MOPA.seurat_type <- FindVariableFeatures(MOPA.seurat_type, selection.method = "vst", nfeatures = 4)
MOPA.seurat_type <- ScaleData(MOPA.seurat_type, features = rownames(MOPA.seurat_type))
variablegene <- VariableFeatures(object = MOPA.seurat_type)
MOPA.seurat_type <- RunPCA(MOPA.seurat_type, features = variablegene,npcs =2,reduction.name = "tSNE",reduction.key = "tSNE_")
MOPA.seurat_type
DimPlot(MOPA.seurat_type)+NoLegend()
#
x <- MOPA@embeddings$PeakMatrix_IterativeLSI_tSNE
y <- x$df
s.TSNE <- cbind(y$`PeakMatrix_IterativeLSI#TSNE_Dimension_1`,
                y$`PeakMatrix_IterativeLSI#TSNE_Dimension_2`)
#
head(s.TSNE)
rownames(s.TSNE) <- rownames(y)
head(s.TSNE)
colnames(s.TSNE) <- c("tSNE_1","tSNE_2")
head(s.TSNE)
class(s.TSNE)
MOPA.seurat_type@reductions$tSNE@cell.embeddings <- s.TSNE
#
#
MOPA.seurat_type$PeakMatrix_IterativeLSI.Cluster.R4 <- MOPA$PeakMatrix_IterativeLSI.Cluster.R4
MOPA.seurat_type@active.ident <- as.factor(MOPA.seurat_type$PeakMatrix_IterativeLSI.Cluster.R4)
DimPlot(MOPA.seurat_type,label = T,repel = T)+NoLegend()
#
#
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C1"="Definitive erythroid lineage")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C2"="Definitive erythroid lineage")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C3"="Definitive erythroid lineage")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C4"="Primitive erythroid lineage")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C5"="Primitive erythroid lineage")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C6"="Definitive erythroid lineage")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C7"="Primitive erythroid lineage")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C8"="Definitive erythroid lineage")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C9"="Hepatocytes")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C10"="Hepatocytes")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C11"="Hepatocytes")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C12"="Hepatocytes")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C13"="Hepatocytes")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C14"="Hepatocytes")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C15"="White blood cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C16"="Endothelial cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C17"="Endothelial cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C18"="Neural progenitor cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C19"="Neural progenitor cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C20"="Neural progenitor cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C21"="Neural progenitor cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C22"="Granule neurons")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C23"="Granule neurons")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C24"="Neural progenitor cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C25"="Neural progenitor cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C26"="Inhibitory interneurons")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C27"="Inhibitory interneurons")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C28"="Inhibitory interneurons")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C29"="Inhibitory neuron progenitors")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C30"="Inhibitory neuron progenitors")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C31"="Postmitotic premature neurons")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C32"="Inhibitory neuron progenitors")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C33"="Excitatory neurons")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C34"="Inhibitory neurons")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C35"="Excitatory neurons")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C36"="Cholinergic neurons")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C37"="Cholinergic neurons")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C38"="Premature oligodendrocyte")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C39"="Premature oligodendrocyte")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C40"="Inhibitory interneurons")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C41"="Oligodendrocyte Progenitors")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C42"="Oligodendrocyte Progenitors")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C43"="Radial glia")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C44"="Radial glia")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C45"="Neural Tube")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C46"="Radial glia")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C47"="Radial glia")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C48"="Notochord cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C49"="Notochord cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C50"="Notochord cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C51"="Notochord cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C52"="Ependymal cell")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C53"="Radial glia")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C54"="Melanocytes")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C55"="Radial glia")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C56"="Neural progenitor cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C57"="Radial glia")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C58"="Isthmic organizer cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C59"="Oligodendrocyte Progenitors")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C60"="Isthmic organizer cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C61"="Radial glia")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C62"="Radial glia")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C63"="Schwann cell precursor")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C64"="Endothelial cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C65"="Chondroctye progenitors")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C66"="Cardiac muscle lineages")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C67"="Intermediate Mesoderm")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C68"="Intermediate Mesoderm")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C69"="Intermediate Mesoderm")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C70"="Intermediate Mesoderm")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C71"="Intermediate Mesoderm")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C72"="Intermediate Mesoderm")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C73"="Epithelial cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C74"="Intermediate Mesoderm")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C75"="Intermediate Mesoderm")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C76"="Jaw and tooth progenitors")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C77"="Jaw and tooth progenitors")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C78"="Jaw and tooth progenitors")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C79"="Jaw and tooth progenitors")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C80"="Limb mesenchyme")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C81"="Connective tissue progenitors")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C82"="Connective tissue progenitors")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C83"="Limb mesenchyme")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C84"="Limb mesenchyme")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C85"="Connective tissue progenitors")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C86"="Connective tissue progenitors")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C87"="Connective tissue progenitors")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C88"="Chondrocytes & osteoblasts")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C89"="Intermediate Mesoderm")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C90"="Chondrocytes & osteoblasts")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C91"="Chondrocytes & osteoblasts")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C92"="Connective tissue progenitors")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C93"="Chondrocytes & osteoblasts")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C94"="Chondrocytes & osteoblasts")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C95"="Chondroctye progenitors")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C96"="Chondrocytes & osteoblasts")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C97"="Chondrocytes & osteoblasts")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C98"="Early mesenchyme")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C99"="Early mesenchyme")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C100"="Chondrocytes & osteoblasts")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C101"="Chondrocytes & osteoblasts")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C102"="Epithelial cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C103"="Epithelial cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C104"="Epithelial cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C105"="Epithelial cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C106"="Epithelial cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C107"="Epithelial cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C108"="Myocytes")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C109"="Connective tissue progenitors")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C110"="Epithelial cells")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C111"="Connective tissue progenitors")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C112"="Chondrocytes & osteoblasts")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C113"="Connective tissue progenitors")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C114"="Myocytes")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C115"="Myocytes")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C116"="Myocytes")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C117"="Myocytes")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C118"="Sensory neurons")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C119"="Sensory neurons")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C120"="Sensory neurons")
MOPA.seurat_type <- RenameIdents(MOPA.seurat_type,"C121"="Sensory neurons")
MOPA.seurat_type$ID <- as.character(MOPA.seurat_type@active.ident)
MOPA$ID <- as.character(MOPA.seurat_type@active.ident)
#
MOPA <- addMotifAnnotations(ArchRProj = MOPA, motifSet = "cisbp", name = "Motif",force=T)
MOPA <- addBgdPeaks(MOPA,force=T)
MOPA <- addDeviationsMatrix(
  ArchRProj = MOPA, 
  peakAnnotation = "Motif",
  force = TRUE
)
#
time.raw <- read.csv("MOPA_Round1_index.csv")
head(time.raw)
dim(time.raw)
R1index <- substring(unlist(strsplit(rownames(MOPA.seurat_type@meta.data),"#"))[seq(2,2*dim(MOPA.seurat_type@meta.data)[1],2)],1,16)
head(R1index)
length(R1index)
BC.meta <- as.data.frame(rownames(MOPA.seurat_type@meta.data)) 
head(BC.meta)
dim(BC.meta)
colnames(BC.meta) <- "BC.total"
BC.meta$R1index <- R1index
head(BC.meta)
rownames(BC.meta) <- BC.meta$BC.total
head(BC.meta)
#
BC.meta.merge <- merge(BC.meta,time.raw,by="R1index")
head(BC.meta.merge)
dim(BC.meta.merge)
rownames(BC.meta.merge) <- BC.meta.merge$BC.total
BC.meta.merge.F <- BC.meta.merge[rownames(MOPA.seurat_type@meta.data),]
head(BC.meta.merge.F )
head(rownames(MOPA.Myocytes.seurat@meta.data))
#
MOPA$embryo <- BC.meta.merge.F$embryo
MOPA$Timepoint <- BC.meta.merge.F$Timepoint
#
MOPA.seurat_type$embryo <- BC.meta.merge.F$embryo
MOPA.seurat_type$Timepoint <- BC.meta.merge.F$Timepoint
#
MOPA
MOPA$BW.ID <- paste(MOPA$ID,"_in_",MOPA$Timepoint,sep = "")
table(MOPA$BW.ID)
length(table(MOPA$BW.ID))
#
getGroupBW(
  ArchRProj = MOPA,
  groupBy = "BW.ID",
  normMethod = "ReadsInTSS",
  tileSize = 100,
  maxCells = 1000,
  ceiling = 4,
  verbose = TRUE,
  threads = getArchRThreads(),
  logFile = createLogFile("getGroupBW")
)
#
features <- list(NeuralprogenitorcellsScore = c("Prmt8", "Gadd45g", "Cdkn1c"),
                 MyocyteScore = c("Neb", "Myh3", "Tpm2", "Acta2"),
                 EndothelialScore = c("Ptprb", "Pecam1", "Vwf"),
                 PrematureoligodendrocyteScore = c("Fut9", "Id4", "Pcdh19", "Cdon", "Emx1"),
                 EpiScore = c("EPCAM","TRP63","GRHL2"),
                 InhibitoryneuronsScore = c("Pax2", "Slc6a5"),
                 IntermediateMesodermScore = c("Wt1", "Mylk"),
                 OligodendrocyteProgenitorsScore = c("Slc17a6", "Sox1","Olig2", "Nkx2-1"))
#
MOPA <- addModuleScore(ArchRProj = MOPA,useMatrix = "GeneScoreMatrix", features =features )
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "cellColData",size = 0.0001, name = "Module.NeuralprogenitorcellsScore", embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave("./Mfeatureplot/NeuralprogenitorcellsScore.png",plot=p,device="png",
       dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotBrowserTrack(
  ArchRProj = MOPA, 
  groupBy = "ID", 
  geneSymbol = "Gapdh", 
  upstream = 50000,
  downstream = 50000
)
grid::grid.newpage()
grid::grid.draw(p$Gapdh)
#
color.use.em <- c("#FDCDAC","#FFFF99","#B3CDE3","#6495ED","#00FFFF", ## 1--5 
                  "#FBB4AE","#F4CAE4","#ABDDA4","#B15928","#D6604D", ## 6--10
                  "#999999","#00BFFF","#CC99FF","#48D1CC","#FFED6F", ## 11--15 
                  "#9E0142","#0099FF","#B3E2CD","#276419","#66C2A5", ## 16--20
                  "#B3EE3A","#6A3D9A","#87CEFA","#ABD9E9","#B2ABD2", ## 21--25 
                  "#0000FF","#FF33CC","#66BD63","#E08214","#DE77AE", ## 26--30
                  "#9970AB","#EE82EE","#99CCCC") 
#
table(MOPA.seurat_type@active.ident)
#
MOPA.seurat_type@active.ident <- as.factor(MOPA.seurat_type$ID)
p <- DimPlot(MOPA.seurat_type,label =F,repel = T,cols =color.use.em,label.size = 6,pt.size = 0.000001,raster=F)+NoLegend()+
  theme(plot.title = element_text(color="black", size=14))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("MOPA.celltype.define.png",plot=p,device="png",dpi=500,units = "cm",width = 16,height = 16)
#
MOPA.seurat_type@active.ident <- as.factor(MOPA.seurat_type$Timepoint)
Timepoint.color <- c("#98304E","#DBCCE0","#48D1CC","#6861A9")
#
p <- DimPlot(MOPA.seurat_type,label =F,repel = T,cols =Timepoint.color,label.size = 6,pt.size = 0.1)+NoLegend()+
  theme(plot.title = element_text(color="black", size=14))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("MOPA.Timepoint.png",plot=p,device="png",dpi=400,units = "cm",width = 13,height = 13)
###
table(MOPA.seurat_type$embryo)
barplot(table(MOPA.seurat_type$embryo),space =0)
#
#
# download the MOCA datasets
MOCA.metadata <- read.csv("MOCA_cell_annotate.csv")
MOCA.metadata[1:2,1:2]
#
MOCA.metadata2 <- subset(MOCA.metadata,doublet_cluster=="FALSE")
dim(MOCA.metadata2)
MOCA.metadata2 <- subset(MOCA.metadata2,development_stage!="9.5")
dim(MOCA.metadata2)
MOCA.m <- as.matrix.data.frame(t(table(MOCA.metadata2$development_stage,MOCA.metadata2$Main_cell_type)))
MOCA.m2 <- as.data.frame(t(table(MOCA.metadata2$development_stage,MOCA.metadata2$Main_cell_type)))
rownames(MOCA.m) <- as.character(MOCA.m2$Var1[1:38])
colnames(MOCA.m) <- c("E10.5","E11.5","E12.5","E13.5")
MOCA.m
MOCA.m <- MOCA.m[rownames(MOAP.ratio),]
#
MOCA.m3 <- MOCA.m
MOCA.m3[,1] <- MOCA.m3[,1]/sum(MOCA.m3[,1])*100
MOCA.m3[,2] <- MOCA.m3[,2]/sum(MOCA.m3[,2])*100
MOCA.m3[,3] <- MOCA.m3[,3]/sum(MOCA.m3[,3])*100
MOCA.m3[,4] <- MOCA.m3[,4]/sum(MOCA.m3[,4])*100
#
MOPA.seurat_type
MOAP.ratio <- as.data.frame(t(table(MOPA.seurat_type$Timepoint,MOPA.seurat_type$ID)))
rownames(MOAP.ratio) <- as.character(MOAP.ratio2$Var1[1:33])
colnames(MOAP.ratio) <- c("E10.5","E11.5","E12.5","E13.5")
MOAP.ratio
#
MOAP.ratio <- MOAP.ratio
MOAP.ratio[,1] <- MOAP.ratio[,1]/sum(MOAP.ratio[,1])*100
MOAP.ratio[,2] <- MOAP.ratio[,2]/sum(MOAP.ratio[,2])*100
MOAP.ratio[,3] <- MOAP.ratio[,3]/sum(MOAP.ratio[,3])*100
MOAP.ratio[,4] <- MOAP.ratio[,4]/sum(MOAP.ratio[,4])*100
#
plot(MOAP.ratio[,1],MOCA.m3[,1])
plot(MOAP.ratio[,2],MOCA.m3[,2])
plot(MOAP.ratio[,3],MOCA.m3[,3])
plot(MOAP.ratio[,4],MOCA.m3[,4])
#
cor.test(MOAP.ratio,MOCA.m3);cor.test(MOAP.ratio,MOCA.m3,method ="spearman")
#
cor.test(MOAP.ratio[,1],MOCA.m3[,1]);cor.test(MOAP.ratio[,1],MOCA.m3[,1],method ="spearman")
cor.test(MOAP.ratio[,2],MOCA.m3[,2]);cor.test(MOAP.ratio[,2],MOCA.m3[,2],method ="spearman")
cor.test(MOAP.ratio[,3],MOCA.m3[,3]);cor.test(MOAP.ratio[,3],MOCA.m3[,3],method ="spearman")
cor.test(MOAP.ratio[,4],MOCA.m3[,4]);cor.test(MOAP.ratio[,4],MOCA.m3[,4],method ="spearman")
#
Both_tech_ratio <- as.data.frame(cbind(c(MOAP.ratio[,1],MOAP.ratio[,2],MOAP.ratio[,3],MOAP.ratio[,4]),c(MOCA.m3[,1],MOCA.m3[,2],MOCA.m3[,3],MOCA.m3[,4])))
colnames(Both_tech_ratio) <- c("ATAC","RNA")
Both_tech_ratio$celltype <- rownames(Both_tech_ratio)
Both_tech_ratio$timepoint <- rep(c("E10.5","E11.5","E12.5","E13.5"),each=33)
Both_tech_ratio$celltype.c <- rep(1:33,4)
#
Both_tech_ratio2 <- Both_tech_ratio
Both_tech_ratio2$timepoint <- as.factor(Both_tech_ratio2$timepoint)
ggplot(Both_tech_ratio2, aes(ATAC, RNA))+
  geom_point(size=2)+
  geom_smooth(method = "lm")+
  geom_point(colour=rep(color.use.em,4))+theme_test()
#
Timepoint.color <- c("#98304E","#DBCCE0","#48D1CC","#6861A9")
cellratio_in_MOPA <- table(MOPA.seurat_type$ID,MOPA.seurat_type$Timepoint)
head(cellratio_in_MOPA)
cellratio_in_MOPA[,1] <- cellratio_in_MOPA[,1]/sum(cellratio_in_MOPA[,1])
cellratio_in_MOPA[,2] <- cellratio_in_MOPA[,2]/sum(cellratio_in_MOPA[,2])
cellratio_in_MOPA[,3] <- cellratio_in_MOPA[,3]/sum(cellratio_in_MOPA[,3])
cellratio_in_MOPA[,4] <- cellratio_in_MOPA[,4]/sum(cellratio_in_MOPA[,4])
#
p <- ggplot(cellratio_in_MOPA,aes(x=Var1,weight=Freq,fill=Var2))+
  geom_bar( position = "stack")+
  scale_fill_manual( values = Timepoint.color)
p
#
save.image("MOPA.RData")
#
# correlation analysis with MOCA single cell RNA-seq in the figure 1I and Figure S1H
color.use.em <- c("#FDCDAC","#FFFF99","#B3CDE3","#6495ED","#00FFFF", ## 1--5 
                  "#FBB4AE","#F4CAE4","#ABDDA4","#B15928","#D6604D", ## 6--10
                  "#999999","#00BFFF","#CC99FF","#48D1CC","#FFED6F", ## 11--15 
                  "#9E0142","#0099FF","#B3E2CD","#276419","#66C2A5", ## 16--20
                  "#B3EE3A","#6A3D9A","#87CEFA","#ABD9E9","#B2ABD2", ## 21--25 
                  "#0000FF","#FF33CC","#66BD63","#E08214","#DE77AE", ## 26--30
                  "#9970AB","#EE82EE","#99CCCC") 
#
#
Get.genescore.mat <- getMatrixFromProject(
  ArchRProj = MOPA,
  useMatrix = "GeneScoreMatrix")
rownames(Get.genescore.mat@assays@data$GeneScoreMatrix) <- Get.genescore.mat@elementMetadata$name
Get.genescore.mat.matrix <- Get.genescore.mat@assays@data$GeneScoreMatrix
head(Get.genescore.mat.matrix)
dim(Get.genescore.mat.matrix)
Get.genescore.mat.matrix.seurat <- CreateSeuratObject(counts = Get.genescore.mat.matrix, assay = "ATAC.GS", project = "mem.GS.matrix", min.cells = 0, min.features = 0)
#
MOPA.metadata <- as.data.frame(MOPA@cellColData)
head(MOPA.metadata)
MOPA.metadata <- MOPA.metadata[rownames(Get.genescore.mat.matrix.seurat@meta.data),]
Get.genescore.mat.matrix.seurat$emid <- MOPA.metadata$emid
#
Get.genescore.mat.matrix.seurat@active.ident <- as.factor(Get.genescore.mat.matrix.seurat$emid)
levels(Get.genescore.mat.matrix.seurat)
levels(Get.genescore.mat.matrix.seurat) <- cluster.averages.sub.cor.rank.levels
#
Get.genescore.mat.matrix.seurat <- NormalizeData(Get.genescore.mat.matrix.seurat)
Get.genescore.mat.matrix.seurat.averages <- AverageExpression(Get.genescore.mat.matrix.seurat, slot = 'data', return.seurat = TRUE)
#
TOME_M.V2.downsample.33.pseudo.F
table(TOME_M.V2.downsample.33.pseudo.F@active.ident)
# find markers for every cluster compared to all remaining cells, report only the positive
# ones
levels(TOME_M.V2.downsample.33.pseudo.F) <- cluster.averages.sub.cor.rank.levels
TOME_M.V2.downsample.33.pseudo.F.markers <- FindAllMarkers(TOME_M.V2.downsample.33.pseudo.F, only.pos = TRUE, min.pct = 0.25)
top3 <- TOME_M.V2.downsample.33.pseudo.F.markers %>% group_by(cluster) %>% top_n(n = 25, wt = avg_log2FC)
DoHeatmap(TOME_M.V2.downsample.33.pseudo.F, features = unique(top3$gene)) + NoLegend()
###
TOME_M.V2.downsample.33.pseudo.F.averages <- AverageExpression(TOME_M.V2.downsample.33.pseudo.F, slot = 'data', return.seurat = TRUE)
TOME_M.V2.downsample.33.pseudo.F.averages
#
top3 <- TOME_M.V2.downsample.33.pseudo.F.markers %>% group_by(cluster) %>% top_n(n = 21, wt = avg_log2FC)
length(unique(top3$gene))
DoHeatmap(TOME_M.V2.downsample.33.pseudo.F.averages, features = unique(top3$gene),draw.lines = F) + NoLegend()+
  scale_fill_gradientn(colors = paletteContinuous("solarExtra"))
#
cluster.averages.sub.cor.rank.levels <- c("Granule neurons","Inhibitory interneurons",
                                          "Postmitotic premature neurons","Inhibitory neuron progenitors",
                                          "Inhibitory neurons","Excitatory neurons","Cholinergic neurons","Sensory neurons","Neural progenitor cells",
                                          "Notochord cells","Neural Tube","Radial glia","Oligodendrocyte Progenitors","Premature oligodendrocyte","Isthmic organizer cells",
                                          "Ependymal cell","Schwann cell precursor","Melanocytes","Connective tissue progenitors",
                                          "Chondrocytes & osteoblasts","Jaw and tooth progenitors","Chondroctye progenitors",
                                          "Intermediate Mesoderm","Early mesenchyme","Limb mesenchyme","Cardiac muscle lineages","Myocytes","Epithelial cells",
                                          "Hepatocytes","Endothelial cells","White blood cells",
                                          "Definitive erythroid lineage","Primitive erythroid lineage")
#
#
Get.genescore.mat.matrix.seurat.averages.sub <- GetAssayData(Get.genescore.mat.matrix.seurat.averages, slot = "scale.data")[unique(top3$gene),]
dim(Get.genescore.mat.matrix.seurat.averages.sub)
colnames(Get.genescore.mat.matrix.seurat.averages.sub) <- paste("ATAC_",colnames(Get.genescore.mat.matrix.seurat.averages.sub),sep = "")
head(Get.genescore.mat.matrix.seurat.averages.sub)
TOME_M.V2.downsample.33.pseudo.F.averages.sub <- GetAssayData(TOME_M.V2.downsample.33.pseudo.F.averages, slot = "scale.data")[unique(top3$gene),]
dim(TOME_M.V2.downsample.33.pseudo.F.averages.sub)
colnames(TOME_M.V2.downsample.33.pseudo.F.averages.sub) <- paste("RNA_",colnames(TOME_M.V2.downsample.33.pseudo.F.averages.sub),sep = "")
#
cc1 <- cbind(Get.genescore.mat.matrix.seurat.averages.sub,
             TOME_M.V2.downsample.33.pseudo.F.averages.sub)
cc1.t <- cor(cc1,method = c("spearman"))
dim(cc1.t)
head(cc1.t)
#
pheatmap(cc1.t)
#
pheatmap::pheatmap(cc1.t[34:66,1:33],border_color = NA,
                   cluster_rows = F,
                   cluster_cols = F,
                   color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100))
###
pheatmap::pheatmap(cc1.t[34:66,1:33],border_color = NA,
                   cluster_rows = T,
                   cluster_cols = T,
                   color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100))
###
##
####
# Based on the results of 33 main cell types
## Then clustering the sub-clusters
# taking granule neuron as an example
MOPA.Granule <- projHeme1[Cells(subset(MOPA.seurat_type,idents="Granule")),]
MOPA.Granule
#
saveArchRProject(ArchRProj = MOPA.Granule, outputDirectory = "MOPA.Granule.sub_cluster", load = F,dropCells = T)
#
#
MOPA.Granule <- addIterativeLSI(
  ArchRProj = MOPA.Granule,
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
MOPA.Granule <- addClusters(
  input = MOPA.Granule,
  reducedDims = "IterativeLSI_Tile",
  method = "Seurat",
  name = "TileMatrix_R5",
  resolution = 5, 
  dimsToUse = 1:100,
  maxClusters= 300,
  force = T,
  sampleCells = NULL
)
#
MOPA.Granule <- addUMAP(
  ArchRProj = MOPA.Granule, 
  reducedDims = "IterativeLSI_Tile", 
  name = "TileMatrix_UMAP_1", 
  nNeighbors = 100, 
  dimsToUse = 1:100, 
  minDist = 0.4, 
  metric = "cosine",force = T
)
#
MOPA.Granule <- addTSNE(
  ArchRProj = MOPA.Granule, 
  reducedDims = "IterativeLSI_Tile", 
  name = "TileMatrix_TSNE_1", 
  perplexity = 50, 
  dimsToUse = 1:100,
  maxIterations = 1500,
  force = T
)
#
#
MOPA.Granule <- addGroupCoverages(ArchRProj = MOPA.Granule, groupBy = "TileMatrix_R5")
addArchRThreads(threads = 15) 
MOPA.Granule <- addReproduciblePeakSet(
  ArchRProj = MOPA.Granule, 
  groupBy = "TileMatrix_R5",
  pathToMacs2 = "/home/sunkeyong/.local/bin/macs2", force=T
)
#
addArchRThreads(threads = 25)
MOPA.Granule <- addPeakMatrix(MOPA.Granule,force = T)
getAvailableMatrices(MOPA.Granule)
#
MOPA.Granule <- addIterativeLSI(
  ArchRProj = MOPA.Granule,
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
MOPA.Granule <- addClusters(
  input = MOPA.Granule,
  reducedDims = "PeakMatrix_LSI",
  method = "Seurat",
  name = "PeakMatrix_LSI_R08",
  resolution = 0.8, 
  dimsToUse = 1:100,
  maxClusters= 300,
  force = T,
  sampleCells = NULL
)
#
MOPA.Granule <- addUMAP(
  ArchRProj = MOPA.Granule, 
  reducedDims = "PeakMatrix_LSI", 
  name = "PeakMatrix_LSI_UMAP_1", 
  nNeighbors = 100, 
  dimsToUse = 1:100, 
  minDist = 0.4, 
  metric = "cosine",force = T
)
#
###### create pseudo seurat file
x.num <- c(rep(1,800000),rep(0,800000),rep(2,800000))
head(x.num)
x.num.sample <- cbind(sample(x.num,length(MOPA.Granule$cellNames)),
                      sample(x.num,length(MOPA.Granule$cellNames)),
                      sample(x.num,length(MOPA.Granule$cellNames)),
                      sample(x.num,length(MOPA.Granule$cellNames)),
                      sample(x.num,length(MOPA.Granule$cellNames)))
head(x.num.sample)
colnames(x.num.sample) <- c("F1","F2","F3","F4","F5")
rownames(x.num.sample) <- MOPA.Granule$cellNames
#
MOPA.Granule.seurat <- CreateSeuratObject(counts = t(x.num.sample), assay = "ATAC", project = "MOPA.Granule", min.cells = 0, min.features = 0)
MOPA.Granule.seurat
MOPA.Granule.seurat <- NormalizeData(MOPA.Granule.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
MOPA.Granule.seurat <- FindVariableFeatures(MOPA.Granule.seurat, selection.method = "vst", nfeatures = 4)
MOPA.Granule.seurat <- ScaleData(MOPA.Granule.seurat, features = rownames(MOPA.Granule.seurat))
variablegene <- VariableFeatures(object = MOPA.Granule.seurat)
MOPA.seurat_type <- RunPCA(MOPA.seurat_type, features = variablegene,npcs =2,reduction.name = "UMAP",reduction.key = "UMAP_")
DimPlot(MOPA.Granule.seurat)
#
x <- MOPA.Granule@embeddings$PeakMatrix_LSI_UMAP_1
y <- x$df
s.UMAP <- cbind(y$`PeakMatrix_LSI#UMAP_Dimension_1`,y$`PeakMatrix_LSI#UMAP_Dimension_2`)
#
head(s.UMAP )
rownames(s.UMAP ) <- rownames(y)
head(s.UMAP )
colnames(s.UMAP ) <- c("PC_1","PC_2")
head(s.UMAP )
class(s.UMAP )
MOPA.Granule.seurat@reductions$pca@cell.embeddings <- s.UMAP 
DimPlot(MOPA.Granule.seurat,label = T,repel = T,pt.size = 0.2)+NoLegend()
# Then define the cell types based on the results of label transfer as shown in the Figure 4.
# The cell type of Granule are shown in the Figure 6
save.image("MOPA.Granule.RData")
#