# Figure 3 
# integrating all cells from 20 stages
setwd("/home/sunkeyong/ZEPA/All_cells")
#
###Zebrafish
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
ZEPA.all.arrow.list <- c("/home/sunkeyong/ZEPA/All_arrow/10hpf.SPATACseq.batch1.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/11hpf.SPATACseq.batch3.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/11hpf.SPATACseq.batch5.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/12hpf.SPATACseq.batch1.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/12hpf.SPATACseq.batch3.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/14hpf.SPATACseq.batch1.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/14hpf.SPATACseq.batch3.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/14hpf.SPATACseq.batch5.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/18hpf.SPATACseq.batch1.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/18hpf.SPATACseq.batch2.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/18hpf.SPATACseq.batch3.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/20hpf.SPATACseq.batch3.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/20hpf.SPATACseq.batch5.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/22hpf.SPATACseq.batch3.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/22hpf.SPATACseq.batch4.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/22hpf.SPATACseq.batch5.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/24hpf.10X.batch1.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/24hpf.10X.batch2.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/24hpf.SPATACseq.batch1.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/24hpf.SPATACseq.batch2.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/24hpf.SPATACseq.batch3.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/24hpf.SPATACseq.batch4.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/24hpf.SPATACseq.batch5.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/30hpf.SPATACseq.batch1.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/30hpf.SPATACseq.batch2.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/30hpf.SPATACseq.batch4.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/30hpf.SPATACseq.batch5.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/34hpf.SPATACseq.batch1.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/34hpf.SPATACseq.batch2.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/38hpf.SPATACseq.batch1.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/38hpf.SPATACseq.batch2.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/42hpf.SPATACseq.batch1.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/42hpf.SPATACseq.batch2.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/48hpf.10X.batch1.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/48hpf.SPATACseq.batch1.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/48hpf.SPATACseq.batch2.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/48hpf.SPATACseq.batch4.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/4hpf.SPATACseq.batch3.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/4hpf.SPATACseq.batch5.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/5hpf.SPATACseq.batch3.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/5hpf.SPATACseq.batch5.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/6hpf.SPATACseq.batch3.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/6hpf.SPATACseq.batch5.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/72hpf.10X.batch1.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/72hpf.SPATACseq.batch1.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/72hpf.SPATACseq.batch2.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/72hpf.SPATACseq.batch4.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/7hpf.SPATACseq.batch3.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/7hpf.SPATACseq.batch5.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/8hpf.SPATACseq.batch1.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/8hpf.SPATACseq.batch3.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/8hpf.SPATACseq.batch5.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/9hpf.SPATACseq.batch3.arrow",
                         "/home/sunkeyong/ZEPA/All_arrow/9hpf.SPATACseq.batch5.arrow")
#
ZEPA <- ArchRProject(ArrowFiles = ZEPA.all.arrow.list$V1,
                     geneAnnotation = GRCz11.geneAnnotation,
                     genomeAnnotation = genomeAnnotation,
                     copyArrows = T)
ZEPA
ZEPA$batch <- unlist(strsplit(unlist(strsplit(ZEPA$cellNames,"#"))[seq(1,2*length(ZEPA$cellNames),2)],"hpf."))[seq(2,2*length(ZEPA$cellNames),2)]
ZEPA$time <- unlist(strsplit(unlist(strsplit(ZEPA$cellNames,"#"))[seq(1,2*length(ZEPA$cellNames),2)],".",fixed=T))[seq(1,3*length(ZEPA$cellNames),3)]
table(ZEPA$batch)
table(ZEPA$time)
ZEPA$batch.plus.time <- paste(ZEPA$batch,ZEPA$time,sep = "_")
#
# loading all cell annotaion results of each stage
load("/home/sunkeyong/ZEPA/OBO/hpf4/zhpf4.seurat.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf5/zhpf5.seurat.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf6/zhpf6.seurat.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf7/zhpf7.seurat.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf8/zhpf8.seurat.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf9/zhpf9.seurat.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf10/zhpf10.seurat.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf11/zhpf11.seurat.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf12/zhpf12.seurat.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf14/zhpf14.seurat.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf18/zhpf18.seurat.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf20/zhpf20.seurat.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf22/zhpf22.seurat.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf24/zhpf24.seurat.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf30/zhpf30.seurat.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf34/zhpf34.seurat.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf38/zhpf38.seurat.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf42/zhpf42.seurat.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf48/zhpf48.seurat.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf72/zhpf72.seurat.RData")
#
ZEPA.allsub.seurat <- merge(x=zhpf4.seurat,y=c(zhpf5.seurat,zhpf6.seurat,zhpf7.seurat,zhpf8.seurat,zhpf9.seurat,
                                               zhpf10.seurat,zhpf11.seurat,zhpf12.seurat,zhpf14.seurat,
                                               zhpf18.seurat,zhpf20.seurat,zhpf22.seurat,zhpf24.seurat,
                                               zhpf30.seurat,zhpf34.seurat,zhpf38.seurat,
                                               zhpf42.seurat,zhpf48.seurat,zhpf72.seurat))
#
ZEPA.allsub.seurat$batch <- unlist(strsplit(unlist(strsplit(Cells(ZEPA.allsub.seurat),"#"))[seq(1,2*length(Cells(ZEPA.allsub.seurat)),2)],"hpf."))[seq(2,2*length(Cells(ZEPA.allsub.seurat)),2)]
ZEPA.allsub.seurat$time <- unlist(strsplit(unlist(strsplit(Cells(ZEPA.allsub.seurat),"#"))[seq(1,2*length(Cells(ZEPA.allsub.seurat)),2)],".",fixed=T))[seq(1,3*length(Cells(ZEPA.allsub.seurat)),3)]
#
table(ZEPA.allsub.seurat$batch)
table(ZEPA.allsub.seurat$time)
table(ZEPA$batch)
table(ZEPA$time)
#
table(ZEPA.allsub.seurat$ID)
length(table(ZEPA.allsub.seurat$ID))
#
ZEPA.allsub.seurameta <- as.data.frame(ZEPA.allsub.seurat@meta.data)
head(ZEPA.allsub.seurameta)
ZEPA.allsub.seurameta2 <- ZEPA.allsub.seurameta[ZEPA$cellNames,]
ZEPA$ID <- as.character(ZEPA.allsub.seurameta2$ID)
length(table(ZEPA$ID))
table(ZEPA$ID,ZEPA$batch.plus.time)
table(ZEPA$ID,ZEPA$time)
#
##### Peak calling using macs2
addArchRThreads(threads = 15) 
ZEPA <- addGroupCoverages(ArchRProj = ZEPA, groupBy = "ID",useLabels=T)
#
ZEPA <- addReproduciblePeakSet(
  ArchRProj = ZEPA, 
  groupBy = "ID",
  genomeSize = 1.5e9,
  excludeChr = c("MT",paste("chr",chr_rm$V1,sep = "")), 
  genomeAnnotation = genomeAnnotation,
  geneAnnotation = GRCz11.geneAnnotation,
  pathToMacs2 = "/home/sunkeyong/.local/bin/macs2", force=T
)
# create peak matrix
ZEPA <- addPeakMatrix(ArchRProj = ZEPA, force =T )
#
ZEPA <- addIterativeLSI(ArchRProj = ZEPA,useMatrix = "PeakMatrix", name = "Peak_LSI_LM1", iterations = 3, 
                        clusterParams = list(resolution = c(3,4), sampleCells = 10000, n.start = 10,maxClusters = 200), 
                        sampleCellsPre = 10000,
                        UMAPParams = list(n_neighbors = 50, min_dist = 0.1, metric = "cosine", verbose =
                                            FALSE, fast_sgd = TRUE),
                        nPlot = 10000,corCutOff = 0.3,
                        varFeatures = 100000, dimsToUse = 1:100,LSIMethod = 1,force=T)
#
ZEPA <- addHarmony(ArchRProj = ZEPA,reducedDims = "Peak_LSI_LM1",groupBy = "batch",
                   name = "Peak_LSI_LM1_Harmony.batch",dimsToUse = 2:100,force = T)
#
ZEPA <- addUMAP(ArchRProj = ZEPA, reducedDims = "Peak_LSI_LM1_Harmony.batch",
                name = "Peak_LSI_LM1_Harmony.batch.UMAP", nNeighbors = 100, 
                dimsToUse = 1:100, minDist = 0.3, metric = "cosine",force = T)
ZEPA <- addUMAP(ArchRProj = ZEPA, reducedDims = "Peak_LSI_LM1",
                name = "Peak_LSI_LM1.UMAP", nNeighbors = 100, 
                dimsToUse = 1:100, minDist = 0.3, metric = "cosine",force = T)
#
ZEPA <- addClusters(input = ZEPA,name = "Peak_LSI_LM1_Harmony.batch.C3",reducedDims = "Peak_LSI_LM1_Harmony.batch",
                    resolution = 3,method = "Seurat", dimsToUse = 1:100,maxClusters= 300,force = T)
#
ZEPA <- addImputeWeights(ZEPA,reducedDims = "Peak_LSI_LM1_Harmony.batch")
#
#
#
markersGS <- getMarkerFeatures(
  ArchRProj = ZEPA, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Peak_LSI_LM1_Harmony.batch.C3",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$C1
#
# create one pseudo seurat file for quick processing 
x.num <- c(rep(1,800000),rep(0,800000),rep(2,800000))
head(x.num)
x.num.sample <- cbind(sample(x.num,length(ZEPA$cellNames)),
                      sample(x.num,length(ZEPA$cellNames)),
                      sample(x.num,length(ZEPA$cellNames)),
                      sample(x.num,length(ZEPA$cellNames)),
                      sample(x.num,length(ZEPA$cellNames)))
head(x.num.sample)
colnames(x.num.sample) <- c("F1","F2","F3","F4","F5")
rownames(x.num.sample) <- ZEPA$cellNames
#
ZEPA.Seurat <- CreateSeuratObject(counts = t(x.num.sample), assay = "ATAC", project = "ZEPA", min.cells = 0, min.features = 0)
ZEPA.Seurat
ZEPA.Seurat <- NormalizeData(ZEPA.Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
ZEPA.Seurat <- FindVariableFeatures(ZEPA.Seurat, selection.method = "vst", nfeatures = 4)
ZEPA.Seurat <- ScaleData(ZEPA.Seurat, features = rownames(ZEPA.Seurat))
variablegene <- VariableFeatures(object = ZEPA.Seurat)
ZEPA.Seurat <- RunPCA(ZEPA.Seurat, features = variablegene,npcs =2,reduction.name = "UMAP")
ZEPA.Seurat <- RunPCA(ZEPA.Seurat, features = variablegene,npcs =2,reduction.name = "UMAPraw")
ZEPA.Seurat
#
ZEPA.Seurat@reductions$UMAP@cell.embeddings[,1] <- ZEPA@embeddings$Peak_LSI_LM1_Harmony.batch.UMAP$df[,1]
ZEPA.Seurat@reductions$UMAP@cell.embeddings[,2] <- ZEPA@embeddings$Peak_LSI_LM1_Harmony.batch.UMAP$df[,2]
ZEPA.Seurat@reductions$UMAPraw@cell.embeddings[,1] <- ZEPA@embeddings$Peak_LSI_LM1.UMAP$df[,1]
ZEPA.Seurat@reductions$UMAPraw@cell.embeddings[,2] <- ZEPA@embeddings$Peak_LSI_LM1.UMAP$df[,2]
#
ZEPA.Seurat$Peak_LSI_LM1_Harmony.batch.C3 <- as.character(ZEPA$Peak_LSI_LM1_Harmony.batch.C3)
ZEPA.Seurat$batch <- as.character(ZEPA$batch)
ZEPA.Seurat$time <- as.character(ZEPA$time)
#
ZEPA.Seurat@active.ident <- as.factor(ZEPA.Seurat$Peak_LSI_LM1.3x_Harmony.batch.C3)
#
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C1"="Hatching gland")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C2"="Hatching gland")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C3"="Hatching gland")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C4"="Hatching gland")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C5"="Periderm")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C6"="Epidermal (foxi3a+)")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C7"="Periderm")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C8"="EVL")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C9"="Periderm")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C10"="Fast muscle cells")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C11"="Fast muscle cells")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C12"="Fast muscle cells")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C13"="Fast muscle cells")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C14"="Fast muscle cells")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C15"="Fast muscle cells")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C16"="Slow muscle cells")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C17"="Slow muscle cells")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C18"="Fast muscle cells")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C19"="Myoblast")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C20"="Tailbud mesoderm")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C21"="Tailbud mesoderm")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C22"="Differentiating neurons (Ganglion)")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C23"="Retina (RGC)")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C24"="Differentiating neurons (Ganglion)")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C25"="Cranial ganglia")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C26"="Differentiating neurons")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C27"="Differentiating neurons")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C28"="Differentiating neurons")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C29"="Differentiating neurons")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C30"="Differentiating neurons (phox2bb+)")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C31"="Differentiating neurons")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C32"="Differentiating neurons")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C33"="Differentiating neurons")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C34"="Differentiating neurons")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C35"="Differentiating neurons")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C36"="Doublets")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C37"="Differentiating neurons")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C38"="Differentiating neurons")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C39"="Differentiating neurons (dlx1a+)")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C40"="Differentiating neurons (eomes+)")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C41"="Differentiating neurons")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C42"="Differentiating neurons")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C43"="Differentiating neurons")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C44"="Tailbud mesoderm")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C45"="Tailbud spinal cord")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C46"="Diencephalon")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C47"="Diencephalon (aplnr2+)")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C48"="Epiblast")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C49"="Epidermal (gbx2+)")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C50"="Telencephalon")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C51"="Telencephalon")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C51"="Neural anterior")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C53"="Telencephalon")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C54"="Midbrain")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C55"="Midbrain")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C56"="Midbrain")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C57"="Midbrain")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C58"="Midbrain ventral")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C59"="Differentiating neurons")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C60"="Fast muscle cells")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C61"="Doublets")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C62"="Retina neurogenesis")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C63"="Hindbrain")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C64"="Hindbrain dorsal")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C65"="Hindbrain ventral")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C66"="Hindbrain dorsal")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C67"="Hindbrain dorsal")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C68"="Neural plate anterior")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C69"="Neural plate posterior")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C70"="Differentiating neurons")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C71"="Differentiating neurons")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C72"="Differentiating neurons")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C73"="Hindbrain ventral")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C74"="Hindbrain ventral")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C75"="Hindbrain ventral")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C76"="Tailbud spinal cord")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C77"="Tailbud spinal cord")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C78"="Hindbrain ventral")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C79"="Tailbud spinal cord")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C80"="Pigment cells")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C81"="Pigment cells")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C82"="Pigment cells")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C83"="Pigment cells")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C84"="Neural crest")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C85"="Neural crest")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C86"="Retina (Muller glia)")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C87"="Optic cup")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C88"="Optic cup")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C89"="Retina neurogenesis")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C90"="Retina neuroblasts")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C91"="Lens")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C92"="Retina (Horizontal cells)")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C93"="Retina (Amacrine)")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C94"="Retina (Amacrine)")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C95"="Retina pigmented epithelium")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C96"="Roofplate")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C97"="Margin involuted")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C98"="Retina (Rods)")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C99"="Retina (Cones)")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C100"="Retina (Photoreceptor precursor cells)")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C101"="Retina (Cone bipolar cells)")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C102"="DEL")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C103"="Margin")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C104"="Differentiating neurons")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C105"="Mesoderm adaxial cells")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C106"="Mesoderm lateral plate")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C107"="Pronephric duct")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C108"="Gut")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C109"="Liver cells")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C110"="Gut")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C111"="YSL")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C112"="YSL")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C113"="Immune cells")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C114"="DEL")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C115"="DEL")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C116"="Cardiac muscle")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C117"="Heart field")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C118"="Immune cells")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C119"="Dorsal margin involuted")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C120"="Mesoderm progenitors")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C121"="Mesoderm lateral plate")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C122"="Mesoderm progenitors")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C123"="Mesoderm lateral plate (tbx1+)")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C124"="Otic placode")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C125"="Olfactory placode")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C126"="Olfactory placode")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C127"="Gut")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C128"="Lateral line primordium")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C129"="Anterior neural ridge")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C130"="Prechordal plate")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C131"="Olfactory placode")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C132"="Otic placode")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C133"="Otic placode")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C134"="Pharyngeal epidermal")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C135"="Epidermal")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C136"="Epidermal")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C137"="Epidermal")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C138"="Blood island")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C139"="Endothelial cells")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C140"="Non dorsal margin involuted")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C141"="Floorplate")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C142"="Endoderm")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C143"="Notochord")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C144"="Mesoderm lateral plate (fli1a+)")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C145"="Neural crest derived mesoderm")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C146"="Heart field")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C147"="Cranial neural crest")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C148"="Pharyngeal arch")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C149"="Cardiac neural crest")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C150"="Heart field")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C151"="Pectoral fin field")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C152"="Vessel progenitor")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C153"="Myoblast")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C154"="Mesoderm lateral plate (tbx1+)")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C155"="Sclerotome")
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat,"C156"="Myotome")
#
levels(ZEPA.Seurat) <- sort(names(table(ZEPA.Seurat@active.ident)))
ZEPA.Seurat$ID.82.celltype <- as.character(ZEPA.Seurat@active.ident)
ZEPA$ID.82.celltype <- ZEPA.Seurat$ID.82.celltype
#
ZEPA.82celltype.color <- c("#D19198","#E0D2DA","#80A794","#A7924E","#77C4E5",
                           "#DFC17A","#7886AD","#CBDDED","#F2ECA5","#DA7A75",
                           "#D8A7C5","#8B6CA2","#C29185","#88B2AE","#BD6870",
                           "#D5E0B2","3DABE68","#74AABC","#EFB687","#DC8C50",
                           "#AB7099","#91766A","#EDD4CB","#E6BF8D","#91C58B",
                           "#ECEAAC","#C6ACA1","#724C73","#B9ACD0","#8C9B64",
                           "#609392","#D5AF6B","#76AC7D","#A1568F","#91598E",
                           "#B75859","#BA7490","#787174","#D5A05F","#ADC27E",
                           "#6D9B5F","#B0A8A8","#9094C5","#E899A1","#BA9D83",
                           "#C4CFA0","#D9C27C","#7EADAB","#917E9C","#BF9DBC",
                           "#D8BF96","#9D8D89","#DE623E","#CA8678","#EAB7C7",
                           "#AD5B4B","#5E82BF","#CACAE1","#D36179","#926FA4",
                           "#5597CF","#99BFD9","#B0C4CD","#4F9598","#53894C",
                           "#C66364","#8DAE93","#526DB1","#526DB1","#AED081",
                           "#699758","#F6DAB4","#C16839","#CB5B4A","#9DCBC1",
                           "#A6ABB1","#A58B82","#DD9345","#67696B","#FDCDAC",
                           "#00BFFF","#9970AB")
DimPlot(ZEPA.Seurat,cols = ZEPA.82celltype.color,
        label = T,repel = T,psize = 0.0000002,
        shuffle=T,raster=T,)+NoLegend()
#
p <- DimPlot(ZEPA.Seurat,cols = ZEPA.82celltype.color,
             label = F,repel = T,psize = 0.0000002,reduction = "UMAP",
             shuffle=T,raster=T,raster.dpi = c(1012, 1012))+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.texx=element_blank(),axis.texy=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
#
ggsave("/home/sunkeyong/ZEPA/All_cells/ZEPA.82celltypes.withBatchCorrection.png",plot=p,device="png",dpi=900,units = "cm",width = 30,height = 30)
#
p <- DimPlot(ZEPA.Seurat,cols = ZEPA.82celltype.color,
             label = F,repel = T,psize = 0.0000002,reduction = "UMAPraw",
             shuffle=T,raster=T,raster.dpi = c(1012, 1012))+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.texx=element_blank(),axis.texy=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
#
ggsave("/home/sunkeyong/ZEPA/All_cells/ZEPA.82celltypes.UMAPraw.png",plot=p,device="png",dpi=900,units = "cm",width = 30,height = 30)
#
#
features.in.ZEPA <- c("glula", "sox19a", "epcam", "tbx1", "sox2", 
                      "rx1", "elavl3", "sox32", "tbx16", "fn1b", 
                      "ntd5", "sox10", "shha", "shhb", "efna2a", 
                      "lrig1", "fezf2", "rx3", "stm", "myf5", 
                      "myod1", "flt1", "nkx3-1", "dmrt2a", "tbx5a", 
                      "hoxb3a","hoxb13a", "gata5", "tmem98", "gata6", 
                      "angptl6", "gck", "nkx2.4a", "nkx2.4b", "eomesa",
                      "foxa2", "sst6", "foxa", "foxg1a", "foxg1b",
                      "kidins220a", "fezf1", "shox2", "otx2a", "olig3")
for (i in 1:length(features.in.ZEPA)) {
  p <- plotEmbedding(ArchRProj = ZEPA, colorBy = "GeneScoreMatrix", name = features.in.ZEPA[i],rastr = T,
                     embedding = "Peak_LSI_LM1_Harmony.batch.UMAP",imputeWeights = getImputeWeights(ZEPA))+NoLegend()+
    theme(plotitle = element_blank())+
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
    theme(axis.texx=element_blank(),axis.texy=element_blank())+
    theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
    theme(panel.background = element_rect( colour = "black", size = 0))+
    theme(panel.border = element_blank())
  ggsave(paste("/home/sunkeyong/ZEPA/All_cells/featureplot/ZEPA.",features.in.ZEPA[i],".genescore.png",sep = ""),plot=p,device="png",dpi=300,units = "cm",width = 15,height = 15)
}
#
p <- plotBrowserTrack(ArchRProj = ZEPA, groupBy = "ID.82.celltype",pal=ZEPA.82celltype.color,
                      geneSymbol = "sp5l", upstream = 10000,downstream = 10000,
                      ylim = c(0.001,0.999999))
grid::grid.newpage()
grid::grid.draw(p$sp5l)
#
p <- plotBrowserTrack(ArchRProj = ZEPA, groupBy = "ID.82.celltype",pal=ZEPA.82celltype.color,
                      geneSymbol = "gad2", upstream = 20000,downstream = 20000,
                      ylim = c(0.001,0.999999))
grid::grid.newpage()
grid::grid.draw(p$gad2)
#
p <- plotBrowserTrack(ArchRProj = ZEPA, groupBy = "ID.82.celltype",pal=ZEPA.82celltype.color,
                      geneSymbol = "slc2a3a", upstream = 30000,downstream = 20000,
                      ylim = c(0.001,0.999999))
grid::grid.newpage()
grid::grid.draw(p$slc2a3a)
#
ZEPA.Seurat
#
time.levels <- c("4hpf","5hpf","6hpf","7hpf",
                 "8hpf","9hpf","10hpf",
                 "11hpf","12hpf","14hpf",
                 "18hpf","20hpf","22hpf","24hpf",
                 "30hpf","34hpf","38hpf",
                 "42hpf","48hpf","72hpf")
ZEPA.timeColor <- c("#FAD6A4","#F7B264","#E1520D","#AD3417",
                    "#F7C7D2","#EA66A1","#E50C84",
                    "#D2E2F2","#227DBD","#004695",
                    "#CFE6C0","#A8D392","#0C9E4A","#006D32",
                    "#C6C4E0","#857BB8","#694799",
                    "#EFBF00","#B78510","#935913")
ZEPA.Seurat@active.ident <- as.factor(ZEPA.Seurat$time)
levels(ZEPA.Seurat) <- time.levels
#
p <- DimPlot(ZEPA.Seurat,cols = ZEPA.timeColor,label = F,repel = T,pt.size = 0.0000002,reduction = "UMAP",
             shuffle=T,raster=T,raster.dpi = c(1012, 1012))+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/All_cells/ZEPA.stages.withBatchCorrection.png",plot=p,device="png",dpi=900,units = "cm",width = 30,height = 30)
#
p <- DimPlot(ZEPA.Seurat,cols = ZEPA.timeColor,label = F,repel = T,pt.size = 0.0000002,reduction = "UMAPraw",
             shuffle=T,raster=T,raster.dpi = c(1012, 1012))+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/All_cells/ZEPA.stages.UMAPraw.png",plot=p,device="png",dpi=900,units = "cm",width = 30,height = 30)
#
ZEPA.Seurat$FRIP <- ZEPA$FRIP
ZEPA.Seurat$nFrags <- ZEPA$nFrags
ZEPA.Seurat$Sample <- ZEPA$Sample
ZEPA.Seurat$DoubletEnrichment <- ZEPA$DoubletEnrichment
#
ZEPA.Seurat$nFrags.log2 <- log2(ZEPA.Seurat$nFrags+1)
ZEPA.Seurat$DoubletEnrichment.log2 <- log2(ZEPA.Seurat$DoubletEnrichment+1)
#
p <- FeaturePlot(ZEPA.Seurat,features = "nFrags.log2",pt.size = 0.0000002,reduction = "UMAP",
                 raster=T,raster.dpi = c(1012, 1012))+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/All_cells/ZEPA.nFrag.log2.png",plot=p,device="png",dpi=900,units = "cm",width = 30,height = 30)
#
p <- FeaturePlot(ZEPA.Seurat,features = "DoubletEnrichment.log2",pt.size = 0.0000002,reduction = "UMAP",
                 raster=T,raster.dpi = c(1012, 1012))+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())+
  theme(legend.title=element_blank())
p
ggsave("/home/sunkeyong/ZEPA/All_cells/ZEPA./DoubletEnrichment.log2.png",plot=p,device="png",dpi=900,units = "cm",width = 30,height = 30)
#
p <- FeaturePlot(ZEPA.Seurat,features = "FRIP",pt.size = 0.0000002,reduction = "UMAP",
                 raster=T,raster.dpi = c(1012, 1012))+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())+
  theme(legend.title=element_blank())
p
ggsave("/home/sunkeyong/ZEPA/All_cells/FRIP.png",plot=p,device="png",dpi=900,units = "cm",width = 30,height = 30)
#
#
ZEPA.Seurat$TSSEnrichment <- ZEPA$TSSEnrichment
p <- FeaturePlot(ZEPA.Seurat,features = "TSSEnrichment",pt.size = 0.0000002,reduction = "UMAP",
                 raster=T,raster.dpi = c(1012, 1012))+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())+
  theme(legend.title=element_blank())
p
ggsave("/home/sunkeyong/ZEPA/All_cells/TSSEnrichment.png",plot=p,device="png",dpi=900,units = "cm",width = 30,height = 30)
#
#
#
ID.for.36M <- read.csv("ID.info.csv")
head(ID.for.36M)
#
renameID <- ID.for.36M$Var4
names(renameID) <- ID.for.36M$Var1
ZEPA.Seurat@active.ident <- as.factor(ZEPA.Seurat$ID.82.celltype)
ZEPA.Seurat <- RenameIdents(ZEPA.Seurat, renameID)
length(table(ZEPA.Seurat@active.ident))
ZEPA.Seurat$ID.36M <- as.factor(ZEPA.Seurat@active.ident)
#
levels.ZEPA.36M <- c("Cranial ganglia","Diencephalon","DNs","Early progenitor",
                     "Endoderm","Endothelial","Epiblast","Epidermal", 
                     "Epidermal (foxi3a+)","Erythroid","Floorplate","Gut",
                     "Hatching gland","Heart","Hindbrain","Immune cells",
                     "Lens","Mesoderm","Mesoderm progenitor","Midbrain",
                     "Muscle","Neural anterior","Neural crest","Neural plate anterior",
                     "Neural plate posterior","Notochord","Olfactory placode","Otic placode",
                     "Periderm","Pharyngeal epidermal","Pronephric duct","Retina",
                     "Roofplate","Tailbud spinal cord","Telencephalon","YSL")
color.for.Main <- c("#FDCDAC","#00BFFF","#9970AB","#F4CAE4",
                    "#0000FF","#3288BD","#66BD63","#FBB4AE",
                    "#FF33CC","#C2A5CF","#00FFFF","#D6604D",
                    "#1A9850","#99CCCC","#EE82EE","#FF7F00",
                    "#9E0142","#B3CDE3","#87CEFA","#6495ED",
                    "#8DA0CB","#B3EE3A","#0099FF","#DE77AE",
                    "#CC99FF","#313695","#276419","#5E4FA2",
                    "#E6AB02","#F6E8C3","#B15928","#66C2A5",
                    "#FF00CC","#9933CC","#48D1CC","#DFC27D")
levels(ZEPA.Seurat) <- levels.ZEPA.36M
p <- DimPlot(ZEPA.Seurat,cols = color.for.Main,label = F,repel = T,pt.size = 0.0000002,reduction = "UMAP",
             shuffle=T,raster=T,raster.dpi = c(1012, 1012))+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/All_cells/ZEPA.Allcell.Main.ID.png",plot=p,device="png",dpi=900,units = "cm",width = 30,height = 30)
#
p <- DimPlot(ZEPA.Seurat,cols = color.for.Main,label = F,repel = T,pt.size = 0.0000002,reduction = "UMAPraw",
             shuffle=T,raster=T,raster.dpi = c(1012, 1012))+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/All_cells/ZEPA.Allcell.Main.ID.UMAPraw.png",plot=p,device="png",dpi=900,units = "cm",width = 30,height = 30)
#
#####
# analyse the cell ratio of hatching gland using Daniocell data
Daniocell.meta <- read.csv("/home/sunkeyong/ZEPA/Daniocell/GSE223922_Sur2023_metadata.tsv",sep = "\t")
head(Daniocell.meta)
table(Daniocell.meta$tissue.name)
table(Daniocell.meta$tissue.figure)
table(Daniocell.meta$stage.group)
table(Daniocell.meta$stage.integer)
table(Daniocell.meta$clust)
table(Daniocell.meta$tissue)
table(Daniocell.meta$original.top.clust)
table(Daniocell.meta$MS_SNR)
table(Daniocell.meta$MS_Seurat)
#
Daniocell.meta.HG <- Daniocell.meta[which(Danioce1.50.52.51.02.00.0ll.meta$tissue.figure=="hatching gland"),]
barplot(as.numeric(table(Daniocell.meta.HG$stage.group))/as.numeric(table(Daniocell.meta$stage.group))[2:14],
        col = New.ZtimeColor[1:13])
barplot(as.numeric(table(Daniocell.meta.HG$stage.integer))/as.numeric(table(Daniocell.meta$stage.integer)))
#
#
#
################################################################
# for analysis in retina cells
# restart
## integrating all cells from 20 stages
setwd("/home/sunkeyong/ZEPA/retina")
#
###Zebrafish
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
BedFiles <- c("/home/sunkeyong/ZEPA/retina/retina.bed.gz")
names(BedFiles) = c("retina")
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
retina <- ArchRProject(ArrowFiles = ArrowFiles,
                       geneAnnotation = GRCz11.geneAnnotation,
                       genomeAnnotation = genomeAnnotation,
                       copyArrows = T)
retina
#
retina <- addIterativeLSI(ArchRProj = retina,useMatrix = "TileMatrix", name = "Tile_LSI_LM1.N1", iterations = 3, 
                          clusterParams = list(resolution = c(1,2), sampleCells = 10000, n.start = 10,maxClusters = 200), 
                          sampleCellsPre = 10000,
                          UMAPParams = list(n_neighbors = 50, min_dist = 0.1, metric = "cosine", verbose =
                                              FALSE, fast_sgd = TRUE),
                          nPlot = 10000,
                          varFeatures = 100000, dimsToUse = 1:100,LSIMethod = 1,force=T)
retina <- addHarmony(ArchRProj = retina,reducedDims = "Tile_LSI_LM1.N1",groupBy = "batch",
                     name = "Tile_LSI_LM1.N1.batch",dimsToUse = 2:100,force = T)
retina <- addUMAP(ArchRProj = retina, reducedDims = "Tile_LSI_LM1.N1.batch",
                  name = "Tile_LSI_LM1.N1.batch.UMAP1", nNeighbors = 50, 
                  dimsToUse = 1:100, minDist = 0.5, metric = "cosine",force = T)
#
#
plotEmbedding(ArchRProj = retina, colorBy = "cellColData", name = "time", embedding = "Tile_LSI_LM1.N1.batch.UMAP1")+NoLegend()
# 
### Create one pseudo seurat file for quick processing
x.num <- c(rep(1,800000),rep(0,800000),rep(2,800000))
head(x.num)
x.num.sample <- cbind(sample(x.num,length(retina$cellNames)),
                      sample(x.num,length(retina$cellNames)),
                      sample(x.num,length(retina$cellNames)),
                      sample(x.num,length(retina$cellNames)),
                      sample(x.num,length(retina$cellNames)))
head(x.num.sample)
colnames(x.num.sample) <- c("F1","F2","F3","F4","F5")
rownames(x.num.sample) <- retina$cellNames
#
retina.seurat <- CreateSeuratObject(counts = t(x.num.sample), assay = "ATAC", project = "z72hpf", min.cells = 0, min.features = 0)
retina.seurat
retina.seurat <- NormalizeData(retina.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
retina.seurat <- FindVariableFeatures(retina.seurat, selection.method = "vst", nfeatures = 4)
retina.seurat <- ScaleData(retina.seurat, features = rownames(retina.seurat))
variablegene <- VariableFeatures(object = retina.seurat)
retina.seurat <- RunPCA(retina.seurat, features = variablegene,npcs =2,reduction.name = "UMAP",reduction.key = "UMAP_")
retina.seurat
#
retina.seurat@reductions$UMAP@cell.embeddings[,1] <- retina@embeddings$Tile_LSI_LM1.N1.batch.UMAP1$df[,1]
retina.seurat@reductions$UMAP@cell.embeddings[,2] <- retina@embeddings$Tile_LSI_LM1.N1.batch.UMAP1$df[,2]
##
load("/home/sunkeyong/ZEPA/All_cells/ZEPA.allsub.seurat.RData")
ZEPA.allsub.seurat
#
ZEPA.allsub.seurat$cells2 <- unlist(strsplit(rownames(ZEPA.allsub.seurat@meta.data),"#"))[seq(2,2*ncol(ZEPA.allsub.seurat),2)]
head(ZEPA.allsub.seurat$cells2)
ZEPA.allsub.seurat.meta <- as.data.frame(ZEPA.allsub.seurat@meta.data)
rownames(ZEPA.allsub.seurat.meta) <- ZEPA.allsub.seurat.meta$cells2
ZEPA.allsub.seurat.meta2 <- ZEPA.allsub.seurat.meta[retina$cellNames2,]
dim(ZEPA.allsub.seurat.meta2)
#
retina$batch <- ZEPA.allsub.seurat.meta2$batch
retina$time <- ZEPA.allsub.seurat.meta2$time
table(retina$batch)
table(retina$time)
#
retina.seurat$batch <- ZEPA.allsub.seurat.meta2$batch
retina.seurat$time <- ZEPA.allsub.seurat.meta2$time
retina.seurat$ID <- ZEPA.allsub.seurat.meta2$ID
retina.seurat$ID.notime <- unlist(strsplit(retina.seurat$ID,":"))[seq(2,2*ncol(retina.seurat),2)]
#
ZEPA.timeColor <- c("#FAD6A4","#F7B264","#E1520D","#AD3417",
                    "#F7C7D2","#EA66A1","#E50C84",
                    "#D2E2F2","#227DBD","#004695",
                    "#CFE6C0","#A8D392","#0C9E4A","#006D32",
                    "#C6C4E0","#857BB8","#694799",
                    "#EFBF00","#B78510","#935913")#
#
retina.seurat@active.ident <- as.factor(retina.seurat$time)
#
p <- DimPlot(retina.seurat,cols = New.ZtimeColor[10:20],pt.size = 0.0001,label = F,repel = T,
             shuffle=T,raster=F)+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/retina/retina.time.png",plot=p,device="png",dpi=400,units = "cm",width = 15,height = 15)
#
retina <- addImputeWeights(retina,reducedDims = "Tile_LSI_LM1.N1")
#
features.in.retina <- c("rx1","vsx1","lmo4a","pax6a","rem1","rx3","cavin2a","ca14","hyal6",
                        "vsx2","hes2.2","foxg1b","pax6b","foxg1d","rlbp1b","neurod1",
                        "crx","gngt2a","gnat1","ompa","tfap2a","rhol","nr2e3","six7","apln","nr2f1b",
                        "apoeb","apoea","rbpms2a","cabp5b","fabp11a","pmela","dct")
for (i in 1:length(features.in.retina)) {
  p <- plotEmbedding(ArchRProj = retina, colorBy = "GeneScoreMatrix", name = features.in.retina[i],rastr = F,
                     embedding = "Tile_LSI_LM1.N1.batch.UMAP1",imputeWeights = getImputeWeights(retina))+NoLegend()+
    theme(plot.title = element_blank())+
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
    theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
    theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
    theme(panel.background = element_rect( colour = "black", size = 0))+
    theme(panel.border = element_blank())
  ggsave(paste("/home/sunkeyong/ZEPA/retina/featureplot/retina.",features.in.retina[i],".genescore.png",sep = ""),plot=p,device="png",dpi=300,units = "cm",width = 15,height = 15)
}
#
### for Nature cell biology source data
Fig3C <- as.data.frame(retina.seurat@reductions$UMAP@cell.embeddings)
colnames(Fig3C) <- c("UMAP1","UMAP2")
Fig3C$time <- as.character(retina.seurat$time)
head(Fig3C)
write.csv(Fig3C,file = "/home/sunkeyong/ZEPA/SourceData/Fig3C.csv")
#
############################################################################################################
############################################################################################################
## for the construction of the develeopmental tree of zebrafish development
# restart R studio
# Taking 4 hpf and 5 hpf as an example
setwd("/home/sunkeyong/ZEPA/Tree/zhpf4to5")
####
###
###Zebrafish
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
# using the arrow files of 4 hpf and 5 hpf
ArrowFiles <- c("/home/sunkeyong/ZEPA/All_arrow/4hpf.SPATACseq.batch3.arrow",
                "/home/sunkeyong/ZEPA/All_arrow/4hpf.SPATACseq.batch5.arrow",
                "/home/sunkeyong/ZEPA/All_arrow/5hpf.SPATACseq.batch3.arrow",
                "/home/sunkeyong/ZEPA/All_arrow/5hpf.SPATACseq.batch5.arrow")
#
zhpf4to5 <- ArchRProject(ArrowFiles = ArrowFiles,
                         geneAnnotation = GRCz11.geneAnnotation,
                         genomeAnnotation = genomeAnnotation,
                         copyArrows = T)
#
zhpf4to5
#
zhpf4to5 <- addIterativeLSI(
  ArchRProj = zhpf4to5,
  useMatrix = "TileMatrix", 
  name = "Tile_IterativeLSI_M1", 
  iterations = 3, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(2,3), 
    sampleCells = 10000, 
    n.start = 10,
    maxClusters = 100
  ), 
  varFeatures = 50000, 
  dimsToUse = 1:100,
  sampleCellsPre = 10000,
  nPlot = 10000,
  LSIMethod = 1,force = T
) 
#
#
zhpf4to5 <- addHarmony(
  ArchRProj = zhpf4to5,
  reducedDims = "Tile_IterativeLSI_M1",
  name = "Tile_Harmony_M1",
  groupBy = "Sample", #######
  dimsToUse = 2:100,force = T
)
#
zhpf4to5 <- addUMAP(
  ArchRProj = zhpf4to5, 
  reducedDims = "Tile_Harmony_M1", 
  name = "Tile_Harmony_M1_UMAP", 
  nNeighbors = 100, 
  dimsToUse = 1:100, 
  minDist = 0.3, 
  metric = "cosine",force = T
)
save.image("zhpf4to5.RData")
#
# create one pseudo seurat file of 4 hpf and 5 hpf for quick processing 
x.num <- c(rep(1,800000),rep(0,800000),rep(2,800000))
head(x.num)
x.num.sample <- cbind(sample(x.num,length(zhpf4to5$cellNames)),
                      sample(x.num,length(zhpf4to5$cellNames)),
                      sample(x.num,length(zhpf4to5$cellNames)),
                      sample(x.num,length(zhpf4to5$cellNames)),
                      sample(x.num,length(zhpf4to5$cellNames)))
head(x.num.sample)
colnames(x.num.sample) <- c("F1","F2","F3","F4","F5")
rownames(x.num.sample) <- zhpf4to5$cellNames
#
zhpf4to5.seurat <- CreateSeuratObject(counts = t(x.num.sample), assay = "ATAC", project = "zhpf4to5", min.cells = 0, min.features = 0)
zhpf4to5.seurat
zhpf4to5.seurat <- NormalizeData(zhpf4to5.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
zhpf4to5.seurat <- FindVariableFeatures(zhpf4to5.seurat, selection.method = "vst", nfeatures = 4)
zhpf4to5.seurat <- ScaleData(zhpf4to5.seurat, features = rownames(zhpf4to5.seurat))
variablegene <- VariableFeatures(object = zhpf4to5.seurat)
zhpf4to5.seurat <- RunPCA(zhpf4to5.seurat, features = variablegene,npcs =2,reduction.name = "UMAP",reduction.key = "UMAP_")
zhpf4to5.seurat
# coordinate transfer
zhpf4to5.seurat@reductions$UMAP@cell.embeddings[,1] <- zhpf4to5@embeddings$Tile_Harmony_M1_UMAP$df[,1]
zhpf4to5.seurat@reductions$UMAP@cell.embeddings[,2] <- zhpf4to5@embeddings$Tile_Harmony_M1_UMAP$df[,2]
#
DimPlot(zhpf4to5.seurat)
#
# loading all cell annotaion results of each stage
load("/home/sunkeyong/ZEPA/OBO/hpf4/zhpf4.seurat.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf5/zhpf5.seurat.RData")
# cell annotation transfer
zhpf4to5.merged.seurat <- merge(x=zhpf4.seurat,y=zhpf5.seurat)
table(zhpf4to5.merged.seurat$Sample)
zhpf4to5.merged.seurat.meta <- as.data.frame(zhpf4to5.merged.seurat@meta.data)
zhpf4.5.merged.seurat.meta <- zhpf4.5.merged.seurat.meta[Cells(zhpf4to5.seurat),]
#
zhpf4to5.seurat$ID <- zhpf4.5.merged.seurat.meta$ID
zhpf4to5.seurat$time <- unlist(strsplit(unlist(strsplit(rownames(zhpf4to5.seurat@meta.data),"#"))[seq(1,2*ncol(zhpf4to5.seurat),2)],".",fixed=T))[seq(1,3*ncol(zhpf4to5.seurat),3)]
#
table(zhpf4to5.seurat$time)
table(zhpf4to5.seurat$ID)
DimPlot(zhpf4to5.seurat,group.by = "ID",split.by = "time",label = T,repel = T)+NoLegend()
#
zhpf4to5.seurat@active.ident <- as.factor(zhpf4to5.seurat$time)
table(zhpf4to5.seurat@active.ident)
zhpf4to5.seurat.hpf4 <- subset(zhpf4to5.seurat,idents="4hpf")
zhpf4to5.seurat.hpf5 <- subset(zhpf4to5.seurat,idents="5hpf")
#
zhpf4to5.seurat.hpf4@active.ident <- as.factor(zhpf4to5.seurat.hpf4$ID)
zhpf4to5.seurat.hpf5@active.ident <- as.factor(zhpf4to5.seurat.hpf5$ID)
#
DimPlot(zhpf4to5.seurat.hpf4,label = T,repel = T,cols = c("#7BCAAF","#F8BDB7","#E0CA8A"))+NoLegend()
DimPlot(zhpf4to5.seurat.hpf5,label = T,repel = T,cols = c("#7BCAAF","#F8BDB7","#4BC5FF","#E0CA8A"))+NoLegend()
#
p <- DimPlot(zhpf4to5.seurat.hpf4,label = F,repel = T,cols = c("#7BCAAF","#F8BDB7","#E0CA8A"),pt.size = 0.001,
             shuffle=T,raster=F)+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/Tree/zhpf4to5/Coembedding.hpf4.png",plot=p,device="png",dpi=500,units = "cm",width=10,height = 10)
#
p <- DimPlot(zhpf4to5.seurat.hpf5,label = F,repel = T,cols = c("#7BCAAF","#F8BDB7","#4BC5FF","#E0CA8A"),pt.size = 0.001,
             shuffle=T,raster=F)+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/Tree/zhpf4to5/Coembedding.hpf5.png",plot=p,device="png",dpi=500,units = "cm",width=10,height = 10)
#
zhpf4to5.seurat@active.ident <- as.factor(zhpf4to5.seurat$ID)
table(zhpf4to5.seurat@active.ident)
p <- DimPlot(zhpf4to5.seurat,label = F,repel = T,cols = c("#7BCAAF","#F8BDB7","#E0CA8A","#7BCAAF","#F8BDB7","#4BC5FF","#E0CA8A"),pt.size = 0.001,
             shuffle=T,raster=F)+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/Tree/zhpf4to5/Coembedding.hpf5_and_hpf4.png",plot=p,device="png",dpi=500,units = "cm",width=10,height = 10)
#
zhpf4to5.seurat.meta <- as.data.frame(zhpf4to5.seurat@meta.data)
zhpf4to5.seurat.meta <- zhpf4to5.seurat.meta[,c(5,6)]
zhpf4to5.seurat.meta$UMAP1 <- zhpf4to5.seurat@reductions$UMAP@cell.embeddings[,1]
zhpf4to5.seurat.meta$UMAP2 <- zhpf4to5.seurat@reductions$UMAP@cell.embeddings[,2]
#
write.csv(zhpf4to5.seurat.meta,file = "/home/sunkeyong/ZEPA/SourceData/Fig3E.csv")
#
###########
# based on UMAP space
zhpf4to5.seurat.4hpf.UMAP <- as.data.frame(zhpf4to5.seurat.4hpf@reductions$UMAP@cell.embeddings)
zhpf4to5.seurat.5hpf.UMAP <- as.data.frame(zhpf4to5.seurat.5hpf@reductions$UMAP@cell.embeddings)
head(zhpf4to5.seurat.4hpf.UMAP)
#
dis1 <- c()
ED.matrix <- c()
for (i in 1:nrow(zhpf4to5.seurat.5hpf.UMAP)) {
  print(i)
  dis1 <- sqrt((zhpf4to5.seurat.4hpf.UMAP$UMAP_1-zhpf4to5.seurat.5hpf.UMAP[i,1])^2+(zhpf4to5.seurat.4hpf.UMAP$UMAP_2-zhpf4to5.seurat.5hpf.UMAP[i,2])^2)
  ED.matrix <- cbind(ED.matrix,dis1)
}
#
rownames(ED.matrix) <- rownames(zhpf4to5.seurat.4hpf@meta.data)
colnames(ED.matrix) <- rownames(zhpf4to5.seurat.5hpf@meta.data)
ED.matrix[1:2,1:2]
dim(ED.matrix)
#
#
zhpf4to5.seurat.4hpf.meta <- as.data.frame(zhpf4to5.seurat.4hpf@meta.data)
zhpf4to5.seurat.5hpf.meta <- as.data.frame(zhpf4to5.seurat.5hpf@meta.data)
zhpf4to5.seurat.4hpf.IDmeta <- as.data.frame(table(zhpf4to5.seurat.4hpf$ID))
zhpf4to5.seurat.5hpf.IDmeta <- as.data.frame(table(zhpf4to5.seurat.5hpf$ID))
#
zhpf4to5.seurat.4hpf.meta$ID <- as.character(zhpf4to5.seurat.4hpf.meta$ID)
zhpf4to5.seurat.5hpf.meta$ID <- as.character(zhpf4to5.seurat.5hpf.meta$ID)
zhpf4to5.seurat.4hpf.IDmeta$Var1 <- as.character(zhpf4to5.seurat.4hpf.IDmeta$Var1)
zhpf4to5.seurat.5hpf.IDmeta$Var1 <- as.character(zhpf4to5.seurat.5hpf.IDmeta$Var1)
#
n.length <- nrow(zhpf4to5.seurat.4hpf.IDmeta)
m.length <- nrow(zhpf4to5.seurat.5hpf.IDmeta)
#
cor.matrix <- c()
nnn <- 500 ###repeat number
for (hh in 1:nnn) {
  ED.matrix2 <- ED.matrix[sample(1:nrow(ED.matrix),nrow(ED.matrix)*0.8),sample(1:ncol(ED.matrix),ncol(ED.matrix)*0.8)]
  print(hh)
  #
  cor.matrix.tmp1 <- matrix(0, n.length, m.length) #n:行数,m:列数
  rownames(cor.matrix.tmp1) <- zhpf4to5.seurat.4hpf.IDmeta$Var1
  colnames(cor.matrix.tmp1) <- zhpf4to5.seurat.5hpf.IDmeta$Var1
  #
  for (ii in 1:m.length) {
    ht.ID.use <- zhpf4to5.seurat.5hpf.IDmeta$Var1[ii]
    cells.use <- intersect(rownames(zhpf4to5.seurat.5hpf.meta[which(zhpf4to5.seurat.5hpf.meta$ID==ht.ID.use),]),colnames(ED.matrix2))
    cell.collect <- c()
    for (jj in 1:length(cells.use)) {
      ED.matrix2.tmp <- ED.matrix2[,cells.use[jj]]
      cell.collect.x <- names(ED.matrix2.tmp[order(ED.matrix2.tmp)])[1:5]
      cell.collect <- c(cell.collect,cell.collect.x)
      Freq <- as.data.frame(table(zhpf4to5.seurat.4hpf.meta[cell.collect,]$ID))
      Freq$ratio <- Freq$Freq/sum(Freq$Freq)
      for (kk in 1:nrow(Freq)) {
        cor.matrix.tmp1[which(rownames(cor.matrix.tmp1)==Freq$Var1[kk]),ii] <- Freq$ratio[kk]
      }
    }
  }
  cor.matrix.tmp2 <- as.vector(cor.matrix.tmp1)
  cor.matrix <- cbind(cor.matrix,cor.matrix.tmp2)
}
head(cor.matrix)
#
# taking the median values of 500 permutations
cor.matrix.median <- apply(cor.matrix, 1, median)
# 
cor.matrix.Final <- c()
##  n.length, m.length
##  n:row number ,m:column number
for (i in 1:m.length) {
  cor.matrix.Final <- cbind(cor.matrix.Final,as.numeric(cor.matrix.median[((i-1)*n.length+1):(i*n.length)]))
}
cor.matrix.Final
#
rownames(cor.matrix.Final) <- rownames(cor.matrix.tmp1)
colnames(cor.matrix.Final) <- colnames(cor.matrix.tmp1)
#
pheatmap(cor.matrix.Final,cluster_rows = F,cluster_cols = F)
# save the edge weight 
# EW: edge weight 
zhpf4to5.UMAP.EW <- as.data.frame(cbind(rep(rownames(cor.matrix.Final),ncol(cor.matrix.Final)),
                                    rep(colnames(cor.matrix.Final),each=nrow(cor.matrix.Final)),
                                    as.vector(cor.matrix.Final)))
colnames(zhpf4to5.cor) <- c("T0","T1","UMAP.EW")
#
save(zhpf4to5.UMAP.EW,file="zhpf4to5.UMAP.EW.RData")
#
###### repeat the above procedures and save all edge weight between two neighbouring stages
# loading all edge weight matrix (19 neighbours)
#
load("/home/sunkeyong/ZEPA/Tree/zhpf4to5/zhpf4to5.UMAP.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf5to6/zhpf5to6.UMAP.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf6to7/zhpf6to7.UMAP.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf7to8/zhpf7to8.UMAP.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf8to9/zhpf8to9.UMAP.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf9to10/zhpf9to10.UMAP.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf10to11/zhpf10to11.UMAP.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf11to12/zhpf11to12.UMAP.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf12to14/zhpf12to14.UMAP.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf14to18/zhpf14to18.UMAP.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf18to20/zhpf18to20.UMAP.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf20to22/zhpf20to22.UMAP.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf22to24/zhpf22to24.UMAP.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf24to30/zhpf24to30.UMAP.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf30to34/zhpf30to34.UMAP.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf34to38/zhpf34to38.UMAP.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf38to42/zhpf38to42.UMAP.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf42to48/zhpf42to48.UMAP.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf48to72/zhpf48to72.UMAP.EW.RData")
#merge all edge weight 
UMAP.EW.merge <- rbind(zhpf4to5.UMAP.EW,zhpf5to6.UMAP.EW,zhpf6to7.UMAP.EW,zhpf7to8.UMAP.EW,
                       zhpf8to9.UMAP.EW,zhpf9to10.UMAP.EW,zhpf10to11.UMAP.EW,
                       zhpf11to12.UMAP.EW,zhpf12to14.UMAP.EW,zhpf14to18.UMAP.EW,
                       zhpf18to20.UMAP.EW,zhpf20to22.UMAP.EW,zhpf22to24.UMAP.EW,
                       zhpf24to30.UMAP.EW,zhpf30to34.UMAP.EW,zhpf34to38.UMAP.EW,
                       zhpf38to42.UMAP.EW,zhpf42to48.UMAP.EW,zhpf48to72.UMAP.EW)
#
dim(UMAP.EW.merge)
head(UMAP.EW.merge)
write.csv(UMAP.EW.merge,file = "/home/sunkeyong/ZEPA/Tree/UMAP.EW.merge.csv")
# used software Cytoscape to generate a directed acyclic graph 
#
bin100.values <- c()
bin100 <- 1:100
bin100 <- bin100/100
UMAP.EW.merge$UMAP.EW <- as.numeric(UMAP.EW.merge$UMAP.EW)
for (i in 1:100) {
  bin100.values <- c(bin100.values,dim(UMAP.EW.merge[which(UMAP.EW.merge$UMAP.EW <= bin100[i]),])[1])
}
#
bin100.values2 <-c()
for (i in 1:99) {
  bin100.values2 <- c(bin100.values2,c(bin100.values[i+1]-bin100.values[i]))
}
#
bin100.values3 <- c(bin100.values[1],bin100.values2)
#
barplot(bin100.values3)
vF <- log2(bin100.values3+1)
names(vF) <- 1:100
barplot(vF) ## for the Figure S3F
#
###################################### 
# based on LSI space
zhpf4to5.lsi <- as.data.frame(zhpf4to5@reducedDims$Tile_Harmony_M1$matDR)
dim(zhpf4to5.lsi)
head(zhpf4to5.lsi)
zhpf4to5.seurat.4hpf.lsi <- zhpf4to5.lsi[rownames(zhpf4to5.seurat.4hpf.UMAP),]
zhpf4to5.seurat.5hpf.lsi <- zhpf4to5.lsi[rownames(zhpf4to5.seurat.5hpf.UMAP),]
#
dis1 <- c()
ED.matrix <- c()
zhpf4to5.seurat.4hpf.lsi2 <- t(zhpf4to5.seurat.4hpf.lsi[,1:30])
zhpf4to5.seurat.5hpf.lsi2 <- t(zhpf4to5.seurat.5hpf.lsi[,1:30])
for (i in 1:nrow(zhpf4to5.seurat.5hpf.lsi)) {
  print(i)
  dis1 <- sqrt(colSums((zhpf4to5.seurat.4hpf.lsi2-zhpf4to5.seurat.5hpf.lsi2[,i])^2))
  ED.matrix <- cbind(ED.matrix,dis1)
}
dim(ED.matrix)
#
rownames(ED.matrix) <- rownames(zhpf4to5.seurat.4hpf@meta.data)
colnames(ED.matrix) <- rownames(zhpf4to5.seurat.5hpf@meta.data)
ED.matrix[1:2,1:2]
dim(ED.matrix)
#
#
zhpf4to5.seurat.4hpf.meta <- as.data.frame(zhpf4to5.seurat.4hpf@meta.data)
zhpf4to5.seurat.5hpf.meta <- as.data.frame(zhpf4to5.seurat.5hpf@meta.data)
zhpf4to5.seurat.4hpf.IDmeta <- as.data.frame(table(zhpf4to5.seurat.4hpf$ID))
zhpf4to5.seurat.5hpf.IDmeta <- as.data.frame(table(zhpf4to5.seurat.5hpf$ID))
#
zhpf4to5.seurat.4hpf.meta$ID <- as.character(zhpf4to5.seurat.4hpf.meta$ID)
zhpf4to5.seurat.5hpf.meta$ID <- as.character(zhpf4to5.seurat.5hpf.meta$ID)
zhpf4to5.seurat.4hpf.IDmeta$Var1 <- as.character(zhpf4to5.seurat.4hpf.IDmeta$Var1)
zhpf4to5.seurat.5hpf.IDmeta$Var1 <- as.character(zhpf4to5.seurat.5hpf.IDmeta$Var1)
#
n.length <- nrow(zhpf4to5.seurat.4hpf.IDmeta)
m.length <- nrow(zhpf4to5.seurat.5hpf.IDmeta)
#
cor.matrix <- c()
nnn <- 500 ###repeat times
for (hh in 1:nnn) {
  ED.matrix2 <- ED.matrix[sample(1:nrow(ED.matrix),nrow(ED.matrix)*0.8),sample(1:ncol(ED.matrix),ncol(ED.matrix)*0.8)]
  print(hh)
  #
  cor.matrix.tmp <- matrix(0, n.length, m.length) #n:行数,m:列数
  rownames(cor.matrix.tmp) <- zhpf4to5.seurat.4hpf.IDmeta$Var1
  colnames(cor.matrix.tmp) <- zhpf4to5.seurat.5hpf.IDmeta$Var1
  #
  for (ii in 1:m.length) {
    ht.ID.use <- zhpf4to5.seurat.5hpf.IDmeta$Var1[ii]
    cells.use <- intersect(rownames(zhpf4to5.seurat.5hpf.meta[which(zhpf4to5.seurat.5hpf.meta$ID==ht.ID.use),]),colnames(ED.matrix2))
    cell.collect <- c()
    for (jj in 1:length(cells.use)) {
      ED.matrix2.tmp <- ED.matrix2[,cells.use[jj]]
      cell.collect.x <- names(ED.matrix2.tmp[order(ED.matrix2.tmp)])[1:5]
      cell.collect <- c(cell.collect,cell.collect.x)
      Freq <- as.data.frame(table(zhpf4to5.seurat.4hpf.meta[cell.collect,]$ID))
      Freq$ratio <- Freq$Freq/sum(Freq$Freq)
      for (kk in 1:nrow(Freq)) {
        cor.matrix.tmp[which(rownames(cor.matrix.tmp)==Freq$Var1[kk]),ii] <- Freq$ratio[kk]
      }
    }
  }
  cor.matrix.tmp2 <- as.vector(cor.matrix.tmp)
  cor.matrix <- cbind(cor.matrix,cor.matrix.tmp2)
}
head(cor.matrix)
#
# taking the median values of 500 permutations based on LSI
cor.matrix.median <- apply(cor.matrix, 1, median)
# 
cor.matrix.Final <- c()
##  n.length, m.length
##  n:row number ,m:column number
for (i in 1:m.length) {
  cor.matrix.Final <- cbind(cor.matrix.Final,as.numeric(cor.matrix.median[((i-1)*n.length+1):(i*n.length)]))
}
cor.matrix.Final
#
rownames(cor.matrix.Final) <- rownames(cor.matrix.tmp1)
colnames(cor.matrix.Final) <- colnames(cor.matrix.tmp1)
#
pheatmap(cor.matrix.Final,cluster_rows = F,cluster_cols = F)
# save the edge weight 
# EW: edge weight 
zhpf4to5.LSI.EW <- as.data.frame(cbind(rep(rownames(cor.matrix.Final),ncol(cor.matrix.Final)),
                                        rep(colnames(cor.matrix.Final),each=nrow(cor.matrix.Final)),
                                        as.vector(cor.matrix.Final)))
colnames(zhpf4to5.cor) <- c("T0","T1","LSI.EW")
#
save(zhpf4to5.LSI.EW,file="zhpf4to5.LSI.EW.RData")
#
###### repeat the above procedures and save all edge weight between two neighbouring stages
# loading all edge weight matrix (19 neighbours)
#
load("/home/sunkeyong/ZEPA/Tree/zhpf4to5/zhpf4to5.LSI.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf5to6/zhpf5to6.LSI.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf6to7/zhpf6to7.LSI.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf7to8/zhpf7to8.LSI.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf8to9/zhpf8to9.LSI.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf9to10/zhpf9to10.LSI.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf10to11/zhpf10to11.LSI.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf11to12/zhpf11to12.LSI.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf12to14/zhpf12to14.LSI.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf14to18/zhpf14to18.LSI.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf18to20/zhpf18to20.LSI.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf20to22/zhpf20to22.LSI.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf22to24/zhpf22to24.LSI.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf24to30/zhpf24to30.LSI.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf30to34/zhpf30to34.LSI.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf34to38/zhpf34to38.LSI.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf38to42/zhpf38to42.LSI.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf42to48/zhpf42to48.LSI.EW.RData")
load("/home/sunkeyong/ZEPA/Tree/zhpf48to72/zhpf48to72.LSI.EW.RData")
#merge all edge weight 
LSI.EW.merge <- rbind(zhpf4to5.LSI.EW,zhpf5to6.LSI.EW,zhpf6to7.LSI.EW,zhpf7to8.LSI.EW,
                       zhpf8to9.LSI.EW,zhpf9to10.LSI.EW,zhpf10to11.LSI.EW,
                       zhpf11to12.LSI.EW,zhpf12to14.LSI.EW,zhpf14to18.LSI.EW,
                       zhpf18to20.LSI.EW,zhpf20to22.LSI.EW,zhpf22to24.LSI.EW,
                       zhpf24to30.LSI.EW,zhpf30to34.LSI.EW,zhpf34to38.LSI.EW,
                       zhpf38to42.LSI.EW,zhpf42to48.LSI.EW,zhpf48to72.LSI.EW)
#
dim(LSI.EW.merge)
head(LSI.EW.merge)
write.csv(LSI.EW.merge,file = "/home/sunkeyong/ZEPA_revision1/0628/Tree/LSI.EW.merge.csv")
# then used software Cytoscape to generate a directed acyclic graph 
# plot the distribution
bin100.values <- c()
bin100 <- 1:100
bin100 <- bin100/100
LSI.EW.merge$LSI.EW <- as.numeric(LSI.EW.merge$LSI.EW)
for (i in 1:100) {
  bin100.values <- c(bin100.values,dim(LSI.EW.merge[which(LSI.EW.merge$LSI.EW <= bin100[i]),])[1])
}
#
bin100.values2 <-c()
for (i in 1:99) {
  bin100.values2 <- c(bin100.values2,c(bin100.values[i+1]-bin100.values[i]))
}
#
bin100.values3 <- c(bin100.values[1],bin100.values2)
#
barplot(bin100.values3)
vF <- log2(bin100.values3+1)
names(vF) <- 1:100
barplot(vF) ## for the Figure S9A
#
# calculate the correlation of the edge weights between based on UMAP and LSI
cor(as.numeric(UMAP.EW.merge$UMAP.EW),as.numeric(LSI.EW.merge$LSI.EW),method = c("pearson"))
##
