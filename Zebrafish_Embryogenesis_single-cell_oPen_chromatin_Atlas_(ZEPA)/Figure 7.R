# Figure 7 
# focused on pigment and notochord development 
# 
# analyse the differentiation trajectory of pigment cells
setwd("/home/sunkeyong/ZEPA/Pigment")
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
addArchRThreads(threads = 1) 
getArchRChrPrefix()
#
v = c("/home/sunkeyong/ZEPA/Pigment/pigment.bed.gz")
names(v) = c("pigment")
#
addArchRThreads(threads = 1) 
getArchRChrPrefix()
#
ArrowFiles <- createArrowFiles(
  inputFiles = v,
  sampleNames = names(v),
  minTSS = 2,
  minFrags = 400,
  offsetPlus = 0,
  offsetMinus = 0,
  geneAnnotation = GRCz11.geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  addTileMat = T,
  addGeneScoreMat = T,
  excludeChr =c("MT",paste("chr",chr_rm$V1,sep = "")))
#
addArchRThreads(threads = 1) 
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)
#
pigment <- ArchRProject(ArrowFiles = ArrowFiles,
                        geneAnnotation = GRCz11.geneAnnotation,
                        genomeAnnotation = genomeAnnotation,
                        copyArrows = T)
pigment
pigment$batch <- unlist(strsplit(unlist(strsplit(pigment$cellNames,"#"))[seq(1,2*length(pigment$cellNames),2)],"hpf."))[seq(2,2*length(pigment$cellNames),2)]
pigment$time <- unlist(strsplit(unlist(strsplit(pigment$cellNames,"#"))[seq(1,2*length(pigment$cellNames),2)],".",fixed=T))[seq(1,3*length(pigment$cellNames),3)]
#
pigment <- addIterativeLSI(ArchRProj = pigment,useMatrix = "TileMatrix", name = "Tile_LSI_LM1", iterations = 3, 
                           clusterParams = list(resolution = c(1,2), sampleCells = 10000, n.start = 10,maxClusters = 200), 
                           sampleCellsPre = 10000,
                           UMAPParams = list(n_neighbors = 50, min_dist = 0.1, metric = "cosine", verbose =
                                               FALSE, fast_sgd = TRUE),
                           nPlot = 10000,
                           varFeatures = 100000, dimsToUse = 1:100,LSIMethod = 1,force=T)
pigment <- addHarmony(ArchRProj = pigment,reducedDims = "Tile_LSI_LM1",groupBy = "batch",
                      name = "Tile_LSI_LM1.batch",dimsToUse = 2:100,force = T)
pigment <- addUMAP(ArchRProj = pigment, reducedDims = "Tile_LSI_LM1.batch",
                   name = "Tile_LSI_LM1.batch.UMAP", nNeighbors = 100, 
                   dimsToUse = 1:100, minDist = 0.3, metric = "cosine",force = T)
pigment <- addClusters(input = pigment,name = "Tile_LSI_LM1.batch.C1",reducedDims = "Tile_LSI_LM1.batch",
                         resolution = 1,method = "Seurat", dimsToUse = 1:100,maxClusters= 300,force = T)
#
plotEmbedding(ArchRProj = pigment, colorBy = "cellColData", name = "time", embedding = "Tile_LSI_LM1.batch.UMAP")+NoLegend()
plotEmbedding(ArchRProj = pigment, colorBy = "cellColData", name = "batch", embedding = "Tile_LSI_LM1.batch.UMAP")+NoLegend()
#
x.num <- c(rep(1,800000),rep(0,800000),rep(2,800000))
head(x.num)
x.num.sample <- cbind(sample(x.num,length(pigment$cellNames)),
                      sample(x.num,length(pigment$cellNames)),
                      sample(x.num,length(pigment$cellNames)),
                      sample(x.num,length(pigment$cellNames)),
                      sample(x.num,length(pigment$cellNames)))
head(x.num.sample)
colnames(x.num.sample) <- c("F1","F2","F3","F4","F5")
rownames(x.num.sample) <- pigment$cellNames
#
pigment.seurat <- CreateSeuratObject(counts = t(x.num.sample), assay = "ATAC", project = "z72hpf", min.cells = 0, min.features = 0)
pigment.seurat
pigment.seurat <- NormalizeData(pigment.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
pigment.seurat <- FindVariableFeatures(pigment.seurat, selection.method = "vst", nfeatures = 4)
pigment.seurat <- ScaleData(pigment.seurat, features = rownames(pigment.seurat))
variablegene <- VariableFeatures(object = pigment.seurat)
pigment.seurat <- RunPCA(pigment.seurat, features = variablegene,npcs =2,reduction.name = "UMAP",reduction.key = "UMAP_")
pigment.seurat
#
pigment.seurat@reductions$UMAP@cell.embeddings[,1] <- pigment.seurat@embeddings$Tile_LSI_LM1.batch.UMAP$df[,1]
pigment.seurat@reductions$UMAP@cell.embeddings[,2] <- pigment.seurat@embeddings$Tile_LSI_LM1.batch.UMAP$df[,2]
#
DimPlot(pigment.seurat,reduction = "UMAP")+NoLegend()
#
pigment.seurat@active.ident <- as.factor(pigment.seurat$Tile_LSI_LM1.batch.C1)
pigment.seurat <- RenameIdents(pigment.seurat,"C1"="xanthophore")
pigment.seurat <- RenameIdents(pigment.seurat,"C2"="xanthophore")
pigment.seurat <- RenameIdents(pigment.seurat,"C3"="xanthophore")
pigment.seurat <- RenameIdents(pigment.seurat,"C4"="melanophore")
pigment.seurat <- RenameIdents(pigment.seurat,"C5"="melanophore")
pigment.seurat <- RenameIdents(pigment.seurat,"C6"="melanophore")
pigment.seurat <- RenameIdents(pigment.seurat,"C7"="iridophore")
pigment.seurat <- RenameIdents(pigment.seurat,"C8"="iridophore")
pigment.seurat <- RenameIdents(pigment.seurat,"C9"="iridophore")
pigment.seurat <- RenameIdents(pigment.seurat,"C10"="iridophore")
pigment.seurat <- RenameIdents(pigment.seurat,"C11"="precursors")
pigment.seurat <- RenameIdents(pigment.seurat,"C12"="precursors")
#
pigment.seurat$ID <- as.character(pigment.seurat@active.ident)
pigment$ID <- pigment.seurat$ID 
############### 
# Peak calling
addArchRThreads(threads = 1) 
pigment <- addGroupCoverages(ArchRProj = pigment, groupBy = "ID")
#
pigment <- addReproduciblePeakSet(
  ArchRProj = pigment, 
  groupBy = "ID",
  genomeSize = 1.5e9,
  excludeChr = c("MT",paste("chr",chr_rm$V1,sep = "")), 
  genomeAnnotation = genomeAnnotation,
  geneAnnotation = GRCz11.geneAnnotation,
  pathToMacs2 = "/home/sunkeyong/.local/bin/macs2", force=T
)
#
pigment.peakset <- getPeakSet(pigment)
#
pigment <- addPeakMatrix(ArchRProj = pigment, force =T )
# Motif annotation
# loading zebrafish motif collected from DANIO-CODE
load("/home/sunkeyong/ZEPA/Custom_danRer11_Motif/danRer11.Motif.RData")
#
pigment <- addMotifAnnotations(ArchRProj = pigment, motifPWMs=danRer11.Motif, name = "Motif")
#
pigment <- addBgdPeaks(pigment)
#
pigment <- addDeviationsMatrix(
  ArchRProj = pigment, 
  peakAnnotation = "Motif",
  force = TRUE
)
#
getAvailableMatrices(pigment)
#
plotVarDev <- getVarDeviations(pigment, name = "MotifMatrix", plot = TRUE)
plotVarDev
###
pigment.seurat@active.ident <- as.factor(pigment.seurat$ID)
DimPlot(pigment.seurat,reduction = "UMAP",cols = c("#C43530","#282D67","#44884A","#7E2F8A"),label = T)
p <- DimPlot(pigment.seurat,reduction = "UMAP",cols = c("#C43530","#282D67","#44884A","#7E2F8A"),
             label = F,pt.size =0.0001)+NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/Pigment/fig_save/pigment_ID.png",plot=p,device="png",dpi=300,units = "cm",width = 14,height = 14)
#
pigment.seurat$time <- pigment$time
##
DimPlot(pigment.seurat,reduction = "UMAP",cols = c("#CFE6C0","#A8D392","#0C9E4A","#006D32",
                                                 "#C6C4E0","#857BB8","#694799",
                                                 "#EFBF00","#B78510","#935913"),label = T)
p <- DimPlot(pigment.seurat,reduction = "UMAP",cols = c("#CFE6C0","#A8D392","#0C9E4A","#006D32",
                                                      "#C6C4E0","#857BB8","#694799",
                                                      "#EFBF00","#B78510","#935913"),
             label = F,pt.size =0.0001)+NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/Pigment/fig_save/pigment_time.png",plot=p,device="png",dpi=300,units = "cm",width = 14,height = 14)
#
#
trajectory12 <- c("C12", "C11","C10","C9","C7","C8")
trajectory13 <- c("C12", "C11","C4","C5","C6")
trajectory14 <- c("C12", "C11","C1","C3","C2")
#
pigment <- addTrajectory(ArchRProj = pigment, name = "CTt12", reducedDims ="Tile_LSI_LM1.batch",
                                   groupBy = "Tile_LSI_LM1.batch.C1",trajectory = trajectory12, 
                                   embedding = "Tile_LSI_LM1.batch.UMAP", force = TRUE)
pigment <- addTrajectory(ArchRProj = pigment, name = "CTt13", reducedDims ="Tile_LSI_LM1.batch",
                                   groupBy = "Tile_LSI_LM1.batch.C1",trajectory = trajectory13, 
                                   embedding = "Tile_LSI_LM1.batch.UMAP", force = TRUE)
pigment <- addTrajectory(ArchRProj = pigment, name = "CTt14", reducedDims ="Tile_LSI_LM1.batch",
                                   groupBy = "Tile_LSI_LM1.batch.C1",trajectory = trajectory14, 
                                   embedding = "Tile_LSI_LM1.batch.UMAP", force = TRUE)
#
p.CTt12 <- plotTrajectory(pigment, trajectory = "CTt12", colorBy = "cellColData", 
                          name = "CTt12",embedding = "Tile_LSI_LM1.batch.UMAP")
p.CTt12[[1]]
p.CTt13 <- plotTrajectory(pigment, trajectory = "CTt13", colorBy = "cellColData", 
                          name = "CTt13",embedding = "Tile_LSI_LM1.batch.UMAP")
p.CTt13[[1]]
p.CTt14 <- plotTrajectory(pigment, trajectory = "CTt14", colorBy = "cellColData", 
                          name = "CTt14",embedding = "Tile_LSI_LM1.batch.UMAP")
p.CTt14[[1]]
#
pigment.seurat$CTt12 <- pigment$CTt12
pigment.seurat$CTt13 <- pigment$CTt13
pigment.seurat$CTt14 <- pigment$CTt14
#
FeaturePlot(pigment.seurat,features = "CTt12")
FeaturePlot(pigment.seurat,features = "CTt13")
FeaturePlot(pigment.seurat,features = "CTt14")
#
features <- list(CTset1Score = intersect(c("paics","gch2","aox5","xdh"),GRCz11.geneAnnotation$genes$symbol),
                 CTset2Score = intersect(c("gstp1","atp6v0ca","pah","slc37a2","tspan36"),GRCz11.geneAnnotation$genes$symbol),
                 CTset3Score = intersect(c("ltk","tfec","apoda.1","pnp4a","gpnmb"),GRCz11.geneAnnotation$genes$symbol))
#
pigment <- addModuleScore(ArchRProj = pigment,useMatrix = "GeneScoreMatrix", features =features )
#
plotEmbedding(ArchRProj = pigment, colorBy = "cellColData",size = 0.0001, 
              name = "Module.CTset1Score", embedding = "Tile_LSI_LM1.batch.UMAP",imputeWeights = getImputeWeights(pigment))
plotEmbedding(ArchRProj = pigment, colorBy = "cellColData",size = 0.0001, 
              name = "Module.CTset2Score", embedding = "Tile_LSI_LM1.batch.UMAP",imputeWeights = getImputeWeights(pigment))
plotEmbedding(ArchRProj = pigment, colorBy = "cellColData",size = 0.0001, 
              name = "Module.CTset3Score", embedding = "Tile_LSI_LM1.batch.UMAP",imputeWeights = getImputeWeights(pigment))
#
p <- plotEmbedding(ArchRProj = pigment, size = 0.0005,plotAs="points",colorBy = "cellColData",
                   name = "Module.CTset3Score", embedding = "Tile_LSI_LM1.batch.UMAP",
                   imputeWeights = getImputeWeights(pigment))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/Pigment/fig_save/CTgenesetscore3.png",plot=p,device="png",
       dpi=300,units = "cm",width = 18,height = 18)
#
p <- plotEmbedding(ArchRProj = pigment, size = 0.0005,plotAs="points",colorBy = "cellColData",
                   name = "Module.CTset2Score", embedding = "Tile_LSI_LM1.batch.UMAP",
                   imputeWeights = getImputeWeights(pigment))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/Pigment/fig_save/CTgenesetscore2.png",plot=p,device="png",
       dpi=300,units = "cm",width = 18,height = 18)
#
p <- plotEmbedding(ArchRProj = pigment, size = 0.0005,plotAs="points",
                   colorBy = "cellColData",quantCut = c(0.00001, 0.81),
                   name = "Module.CTset1Score", embedding = "Tile_LSI_LM1.batch.UMAP",
                   imputeWeights = getImputeWeights(pigment))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/Pigment/fig_save/CTgenesetscore1.png",plot=p,device="png",
       dpi=300,units = "cm",width = 18,height = 18)
#
p <- plotEmbedding( ArchRProj = pigment,size = 0.08,plotAs="points", 
                    colorBy = "GeneScoreMatrix", name = "sox2",
                    pal =paletteContinuous("solarExtra"),
                    embedding = "Tile_LSI_LM1.batch.UMAP",
                    imputeWeights = getImputeWeights(pigment))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/Pigment/fig_save/CT.sox2.genescore.png",plot=p,device="png",
       dpi=300,units = "cm",width = 18,height = 18)
#
p <- plotBrowserTrack(ArchRProj = pigment, groupBy = "ID",pal = c("#C43530","#282D67","#44884A","#7E2F8A"),
                      geneSymbol = "apoda.1",ylim = c(0.1,1), upstream = 20000,downstream = 17000,verbose =F,plotSummary = c("bulkTrack", "geneTrack"))
grid::grid.draw(p$apoda.1)
plotPDF(plotList = p, name = "track-apoda.1", ArchRProj = pigment,addDOC = FALSE, width = 5, height = 5)
#
p <- plotBrowserTrack(ArchRProj = pigment, groupBy = "ID",pal = c("#C43530","#282D67","#44884A","#7E2F8A"),
                      geneSymbol = "gch2", upstream =25000,downstream = 20000,verbose =F,plotSummary = c("bulkTrack", "geneTrack"))
grid::grid.draw(p$gch2)
plotPDF(plotList = p, name = "track-gch2", ArchRProj = pigment,addDOC = FALSE, width = 5, height = 5)
#
p <- plotBrowserTrack(ArchRProj = pigment, groupBy = "ID",pal = c("#C43530","#282D67","#44884A","#7E2F8A"),
                      geneSymbol = "pmela", upstream = 30000,downstream = 30000,plotSummary = c("bulkTrack", "geneTrack"))
grid::grid.draw(p$pmela)
plotPDF(plotList = p, name = "track-pmela", ArchRProj = pigment,addDOC = FALSE, width = 5, height = 5)
##
p <- plotBrowserTrack(ArchRProj = pigment, groupBy = "ID",pal = c("#C43530","#282D67","#44884A","#7E2F8A"),
                      geneSymbol = "dct", upstream = 10000,downstream = 40000,verbose =F,plotSummary = c("bulkTrack", "geneTrack"))
grid::grid.draw(p$dct)
plotPDF(plotList = p, name = "track-dct", ArchRProj = pigment,addDOC = FALSE, width = 5, height = 5)
#
p <- plotBrowserTrack(ArchRProj = pigment, groupBy = "ID",pal = c("#C43530","#282D67","#44884A","#7E2F8A"),
                      geneSymbol = "tyr", upstream = 5000,downstream = 40000,verbose =F,plotSummary = c("bulkTrack", "geneTrack"))
grid::grid.draw(p$tyr)
plotPDF(plotList = p, name = "track-tyr", ArchRProj = pigment,addDOC = FALSE, width = 5, height = 5)
#
p <- plotBrowserTrack(ArchRProj = pigment, groupBy = "ID",pal = c("#C43530","#282D67","#44884A","#7E2F8A"),
                      geneSymbol = "sox2", ylim = c(0.1,1),upstream = 40000,downstream = 90000,verbose =F,plotSummary = c("bulkTrack", "geneTrack"))
grid::grid.draw(p$sox2)
plotPDF(plotList = p, name = "track-sox2", ArchRProj = pigment,addDOC = FALSE, width = 5, height = 5)
#
pigment <- addImputeWeights(pigment,reducedDims = "Tile_LSI_LM1.N1.batch")
plotEmbedding( ArchRProj = pigment, colorBy = "GeneScoreMatrix", name = "sox10", 
               quantCut = c(0.01, 0.4),
               embedding = "Tile_LSI_LM1.N1.batch.UMAP3",imputeWeights = getImputeWeights(pigment))
plotEmbedding( ArchRProj = pigment, colorBy = "GIM", name = "sox10",
               log2Norm =T,quantCut = c(0.01, 0.6), embedding = "Tile_LSI_LM1.N1.batch.UMAP3",imputeWeights = getImputeWeights(pigment))
#
# perform peak to gene correlation analysis
pigment <- addPeak2GeneLinks(
  ArchRProj = pigment,
  reducedDims = "Tile_LSI_LM1.batch",
  useMatrix = "GeneScoreMatrix",
  dimsToUse = 1:50,
  scaleDims = NULL)
#
p2g <- getPeak2GeneLinks(
  ArchRProj = pigment,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)
p2g
# track view with information about peak to gene links
# taking sox10 and pmela as examples
p <- plotBrowserTrack(ArchRProj = pigment, groupBy = "ID",loops = getPeak2GeneLinks(pigment),
                      pal = c("#C43530","#282D67","#44884A","#7E2F8A"),
                      geneSymbol = "pmela", upstream = 30000,downstream = 30000,verbose =F)
grid::grid.draw(p$pmela)
plotPDF(plotList = p, name = "track-p2g-pmela", ArchRProj = pigment,addDOC = FALSE, width = 5, height = 5)
##
p <- plotBrowserTrack(ArchRProj = pigment, groupBy = "IDx",loops = getPeak2GeneLinks(pigment),
                      pal = c("#C43530","#282D67","#44884A","#7E2F8A"),
                      geneSymbol = "sox10", upstream = 40000,downstream = 20000,verbose =F)
grid::grid.draw(p$sox10)
plotPDF(plotList = p, name = "track-p2g-sox10-2", ArchRProj = pigment,addDOC = FALSE, width = 5, height = 5)
#
# calculate  the peak To gene number 
hist(table(p2g$idxATAC)) 
median(table(p2g$idxATAC)) 
hist(table(p2g$idxRNA)) 
median(table(p2g$idxRNA))
#
p2geneDF <- metadata(pigment@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
p2geneDF
p2geneDF.tmp <- as.data.frame(p2geneDF)
p2geneDF.tmp
p2geneDF.tmp <- p2geneDF.tmp[which(p2geneDF.tmp$Correlation > 0.5),]
p2geneDF.tmp <- p2geneDF.tmp[which(p2geneDF.tmp$VarQATAC > 0.35),]
p2geneDF.tmp <- p2geneDF.tmp[which(p2geneDF.tmp$VarQRNA > 0.35),]
p2geneDF.tmp <- p2geneDF.tmp[which(p2geneDF.tmp$FDR < 1e-04),]
p2geneDF.tmp
dim(p2geneDF.tmp)
#
p <- plotPeak2GeneHeatmap(ArchRProj = pigment, groupBy = "ID",k = 15,
                          corCutOff = 0.5,
                          varCutOffATAC = 0.35,
                          varCutOffRNA = 0.35)
p
#
p2g.M <- plotPeak2GeneHeatmap(ArchRProj = pigment, groupBy = "ID",
                              k = 15,returnMatrices =T,
                              corCutOff = 0.5,
                              varCutOffATAC = 0.35,
                              varCutOffRNA = 0.35,
                              nPlot =48788)
p2g.M
#

library(ComplexHeatmap)
##
dim(p2g.M$ATAC$matrix)
dim(p2g.M$RNA$matrix)
#
p2g.M@listData$ATAC$kmeansId
#
pAm <- p2g.M$ATAC$matrix
head(pAm)
tail(pAm)
pAm.row.rank <- c()
for(i in 1:15){
  x <- which(p2g.M$ATAC$kmeansId==i)
  pAm.row.rank <- c(pAm.row.rank,x)
}
length(pAm.row.rank)
pAm <- pAm[pAm.row.rank,]
head(pAm)#
#
p2gkmeansId <- as.data.frame(table(p2g.M$ATAC$kmeansId))
p2gkmeansId$Freq[c(c(3,1,2,5,7,4,6),10,9,13,14,15,8,11,12)]
#
pAm.row.rank2 <- c(rownames(pAm)[1:6436],
                   sample(rownames(pAm)[6437:23465],17029),
                   rownames(pAm)[23466:48788])
pAm <- pAm[pAm.row.rank2,]
dim(pAm)
head(pAm)
#
col_fun = colorRamp2(seq(-127,128,1)/64, paletteContinuous("solarExtra"))
split.module <- rep(1:15,p2gkmeansId$Freq)
split.module <- factor(split.module,levels = 1:15)
#
zC1.length <- c(which(split.module=="10"),which(split.module=="12"),which(split.module=="13"))
zC1.length <- zC1.length[sample(length(zC1.length),length(zC1.length))]
zC2.length <- c(which(split.module=="8"))
zC2.length <- zC2.length[sample(length(zC2.length),length(zC2.length))]
zC3.length <- c(which(split.module=="4"),which(split.module=="5"),which(split.module=="6"))
zC3.length <- zC3.length[sample(length(zC3.length),length(zC3.length))]
zC4.length <- c(which(split.module=="1"),which(split.module=="2"),which(split.module=="3"))
zC4.length <- zC4.length[sample(length(zC4.length),length(zC4.length))]
#
split.module2 <- rep(1:4,c(length(zC1.length),length(zC2.length),length(zC3.length),length(zC4.length)))
split.module2 <- factor(split.module2,levels = 1:4)
pAm.x <- pAm[c(zC1.length,zC2.length,zC3.length,zC4.length),]
p2g.ATAC.plot <- Heatmap(pAm.x,
                         col = col_fun,
                         cluster_rows = F,
                         cluster_columns = F,
                         show_column_names = F,
                         show_row_names = F,
                         row_split = split.module2,
                         gap = unit(0.5, "mm"), 
                         border = "black",
                         use_raster = T,
                         #top_annotation=top_color ## 控制cluster的annotation
)
draw(p2g.ATAC.plot, heatmap_legend_side = "right", show_annotation_legend = T)
#
CzC1 <- 363:398
CzC2 <- 263:298
CzC3 <- 154:188
CzC4 <- 1:40
#
dim(pAm.x)
pAm.x2 <- pAm.x[,c(CzC1,CzC2,CzC3,CzC4)]
p2g.ATAC.plot <- Heatmap(pAm.x2,
                         col = col_fun,
                         cluster_rows = F,
                         cluster_columns = F,
                         show_column_names = F,
                         show_row_names = F,
                         row_split = split.module2,
                         gap = unit(0.5, "mm"), 
                         border = "black",
                         use_raster = T,
                         #top_annotation=top_color ## 控制cluster的annotation
)
draw(p2g.ATAC.plot, heatmap_legend_side = "right", show_annotation_legend = T)
#
pRm <- p2g.M$RNA$matrix
pRm <- pRm[rownames(pAm),]
head(pRm)
#
#
RNA_col_fun = colorRamp2(seq(-127,128,1)/64, paletteContinuous("blueYellow"))
#
p2g.RNA.plot <- Heatmap(pRm[rownames(pAm.x2),colnames(pAm.x2)],
                        col = RNA_col_fun,
                        cluster_rows = F,
                        cluster_columns = F,
                        show_column_names = F,
                        show_row_names = F,
                        row_split = split.module2,
                        gap = unit(0.5, "mm"), 
                        border = "black"
                        #top_annotation=top_color, 
                        #col = c("#FFFFFF","#4699C8")
)
draw(p2g.RNA.plot, heatmap_legend_side = "right", show_annotation_legend = T)
#
p2g.M.linkers.meta <- as.data.frame(p2g.M$Peak2GeneLinks)
p2g.M.linkers.meta2 <- p2g.M.linkers.meta[sort(unique(rownames(pAm.x2))),]
dim(p2g.M.linkers.meta2)
write.table(p2g.M.linkers.meta2,file = "/home/sunkeyong/ZEPA/Pigment/pigment.p2g.M.linkers.meta.csv",quote=F, sep = "\t",row.names = F,col.names = T)
#
hist(table(p2g.M.linkers.meta2$gene))
median(table(p2g.M.linkers.meta2$gene))
xx.peak <- table(table(p2g.M.linkers.meta2$peak))
barplot(xx.peak,space = 0)
median(table(p2g.M.linkers.meta2$peak))
#3 & 1
hist(table(p2g.M.linkers.meta2$gene))
median(table(p2g.M.linkers.meta2$gene))
hist(table(p2geneDF.tmp$peakName),breaks=150)
median(table(p2geneDF.tmp$peakName))
#
p2g.peak.num <- as.data.frame(table(p2g.M.linkers.meta2$gene))
p2g.peak.num[which(p2g.peak.num$Freq > 40 ),2] <- 40
max(p2g.peak.num$Freq)
hist(p2g.peak.num$Freq,breaks = 31)
hist(p2g.peak.num$Freq,freq = 1)
#
barplot(table(p2g.peak.num$Freq),space = 0)
#
p2g.gene.num <- as.data.frame(table(p2geneDF.tmp$peakName))
p2g.gene.num[which(p2g.gene.num$Freq > 7 ),2] <- 7
hist(p2g.gene.num$Freq,breaks = 8)
hist(p2g.gene.num$Freq,breaks = 1:8)
##
barplot(table(p2g.gene.num$Freq),space = 0)
#
split.module.dataframe <- as.data.frame(table(split.module))
split.module.dataframe 
#
p2g.M.Peak2GeneLinks <- as.data.frame(p2g.M$Peak2GeneLinks)
head(p2g.M.Peak2GeneLinks)
#
for (x in 1:15) {
  peak <- unique(p2g.M.Peak2GeneLinks[rownames(pAm)[(sum(split.module.dataframe$Freq[1:(x-1)])+1):sum(split.module.dataframe$Freq[1:x])],]$peak)
  x.p1 <- unlist(strsplit(peak,":"))[seq(1,2*length(peak),2)]
  x.p2 <- unlist(strsplit(peak,":"))[seq(2,2*length(peak),2)]
  x.p3 <- as.data.frame(cbind(x.p1,unlist(strsplit(x.p2,"-"))[seq(1,2*length(Cluster3.peak1),2)],unlist(strsplit(x.p2,"-"))[seq(2,2*length(Cluster3.peak1),2)]))
  write.table(x.p3,file = paste("/home/sunkeyong/ZEPA/Pigment/P2G/",x,".peak.csv",sep = ""),quote=F, sep = "\t",row.names = F,col.names = F)
}
#
write.csv(pAm.x2,file = "/home/sunkeyong/ZEPA/SourceData/Fig7F_left.csv")
write.csv(pRm[rownames(pAm.x2),colnames(pAm.x2)],file = "/home/sunkeyong/ZEPA/SourceData/Fig7F_right.csv")
#
pigment.seurat.meta <- as.data.frame(cbind(pigment.seurat.meta$time,
                                           pigment.seurat.meta$ID,
                                           pigment.seurat.meta$UMAP1,
                                           pigment.seurat.meta$UMAP2))
colnames(pigment.seurat.meta) <- c("time","celltype","UMAP1","UMAP2")
rownames(pigment.seurat.meta) <- rownames(pigment.seurat.meta)
head(pigment.seurat.meta)
#
write.csv(pigment.seurat.meta,file = "/home/sunkeyong/ZEPA/SourceData/Fig7_UMAP.csv")
#
write.csv(table(p2g.peak.num$Freq),file = "/home/sunkeyong/ZEPA/SourceData/Fig7E_Down.csv")
write.csv(table(p2g.gene.num$Freq),file = "/home/sunkeyong/ZEPA/SourceData/Fig7E_UP.csv")
#
######
############################################################################
#######################################
# Figure 7 
# analyse the differentiation trajectory of notochord
setwd("/home/sunkeyong/ZEPA/notochord")
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
addArchRThreads(threads = 1) 
getArchRChrPrefix()
#
v = c("/home/sunkeyong/ZEPA/notochord/notochord.bed.gz")
names(v) = c("notochord")
#
addArchRThreads(threads = 1) 
getArchRChrPrefix()
#
ArrowFiles <- createArrowFiles(
  inputFiles = v,
  sampleNames = names(v),
  minTSS = 2,
  minFrags = 400,
  offsetPlus = 0,
  offsetMinus = 0,
  geneAnnotation = GRCz11.geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  addTileMat = T,
  addGeneScoreMat = T,
  excludeChr =c("MT",paste("chr",chr_rm$V1,sep = "")))
#
addArchRThreads(threads = 1) 
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)
#
notochord <- ArchRProject(ArrowFiles = ArrowFiles,
                          geneAnnotation = GRCz11.geneAnnotation,
                          genomeAnnotation = genomeAnnotation,
                          copyArrows = T)
notochord
notochord$batch <- unlist(strsplit(unlist(strsplit(notochord$cellNames,"#"))[seq(1,2*length(notochord$cellNames),2)],"hpf."))[seq(2,2*length(notochord$cellNames),2)]
notochord$time <- unlist(strsplit(unlist(strsplit(notochord$cellNames,"#"))[seq(1,2*length(notochord$cellNames),2)],".",fixed=T))[seq(1,3*length(notochord$cellNames),3)]
#
notochord <- addIterativeLSI(ArchRProj = notochord,useMatrix = "TileMatrix", name = "Tile_LSI_LM1", iterations = 3, 
                             clusterParams = list(resolution = c(1,2), sampleCells = 10000, n.start = 10,maxClusters = 200), 
                             sampleCellsPre = 10000,
                             UMAPParams = list(n_neighbors = 50, min_dist = 0.1, metric = "cosine", verbose =
                                                 FALSE, fast_sgd = TRUE),
                             nPlot = 10000,
                             varFeatures = 100000, dimsToUse = 1:100,LSIMethod = 1,force=T)
notochord <- addHarmony(ArchRProj = notochord,reducedDims = "Tile_LSI_LM1",groupBy = "batch",
                        name = "Tile_LSI_LM1.batch",dimsToUse = 2:100,force = T)
notochord <- addUMAP(ArchRProj = notochord, reducedDims = "Tile_LSI_LM1.batch",
                     name = "Tile_LSI_LM1.batch.UMAP", nNeighbors = 100, 
                     dimsToUse = 1:100, minDist = 0.3, metric = "cosine",force = T)
notochord <- addClusters(input = notochord,name = "Tile_LSI_LM1.batch.C1",reducedDims = "Tile_LSI_LM1.batch",
                         resolution = 1,method = "Seurat", dimsToUse = 1:100,maxClusters= 300,force = T)
#
plotEmbedding(ArchRProj = notochord, colorBy = "cellColData", name = "time", embedding = "Tile_LSI_LM1.batch.UMAP")+NoLegend()
plotEmbedding(ArchRProj = notochord, colorBy = "cellColData", name = "batch", embedding = "Tile_LSI_LM1.batch.UMAP")+NoLegend()
#
x.num <- c(rep(1,800000),rep(0,800000),rep(2,800000))
head(x.num)
x.num.sample <- cbind(sample(x.num,length(notochord$cellNames)),
                      sample(x.num,length(notochord$cellNames)),
                      sample(x.num,length(notochord$cellNames)),
                      sample(x.num,length(notochord$cellNames)),
                      sample(x.num,length(notochord$cellNames)))
head(x.num.sample)
colnames(x.num.sample) <- c("F1","F2","F3","F4","F5")
rownames(x.num.sample) <- notochord$cellNames
#
notochord.seurat <- CreateSeuratObject(counts = t(x.num.sample), assay = "ATAC", project = "z72hpf", min.cells = 0, min.features = 0)
notochord.seurat
notochord.seurat <- NormalizeData(notochord.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
notochord.seurat <- FindVariableFeatures(notochord.seurat, selection.method = "vst", nfeatures = 4)
notochord.seurat <- ScaleData(notochord.seurat, features = rownames(notochord.seurat))
variablegene <- VariableFeatures(object = notochord.seurat)
notochord.seurat <- RunPCA(notochord.seurat, features = variablegene,npcs =2,reduction.name = "UMAP",reduction.key = "UMAP_")
notochord.seurat
#
notochord.seurat@reductions$UMAP@cell.embeddings[,1] <- notochord.seurat@embeddings$Tile_LSI_LM1.batch.UMAP$df[,1]
notochord.seurat@reductions$UMAP@cell.embeddings[,2] <- notochord.seurat@embeddings$Tile_LSI_LM1.batch.UMAP$df[,2]
#
DimPlot(notochord.seurat,reduction = "UMAP")+NoLegend()
#
notochord <- addImputeWeights(notochord,reducedDims = "Tile_LSI_LM1.batch")
plotEmbedding(ArchRProj = notochord, colorBy = "GeneScoreMatrix", name = "epcam", 
              embedding = "Tile_LSI_LM1.batch.UMAP",imputeWeights = getImputeWeights(notochord))
plotEmbedding(ArchRProj = notochord, colorBy = "GeneScoreMatrix", name = "noto", 
              embedding = "Tile_LSI_LM1.batch.UMAP",imputeWeights = getImputeWeights(notochord))
plotEmbedding(ArchRProj = notochord, colorBy = "GeneScoreMatrix", name = "col2a1a", 
              embedding = "Tile_LSI_LM1.batch.UMAP",imputeWeights = getImputeWeights(notochord))
plotEmbedding(ArchRProj = notochord, colorBy = "GeneScoreMatrix", name = "col11a2", 
              embedding = "Tile_LSI_LM1.batch.UMAP",imputeWeights = getImputeWeights(notochord))
plotEmbedding(ArchRProj = notochord, colorBy = "GeneScoreMatrix", name = "sox19a", 
              embedding = "Tile_LSI_LM1.batch.UMAP",imputeWeights = getImputeWeights(notochord))
plotEmbedding(ArchRProj = notochord, colorBy = "GeneScoreMatrix", name = "pacsin3", 
              embedding = "Tile_LSI_LM1.batch.UMAP",imputeWeights = getImputeWeights(notochord))
plotEmbedding(ArchRProj = notochord, colorBy = "GeneScoreMatrix", name = "ngs", 
              embedding = "Tile_LSI_LM1.batch.UMAP",imputeWeights = getImputeWeights(notochord))
plotEmbedding(ArchRProj = notochord, colorBy = "GeneScoreMatrix", name = "cnmd", 
              embedding = "Tile_LSI_LM1.batch.UMAP",imputeWeights = getImputeWeights(notochord))
plotEmbedding(ArchRProj = notochord, colorBy = "GeneScoreMatrix", name = "sox9a", 
              embedding = "Tile_LSI_LM1.batch.UMAP",imputeWeights = getImputeWeights(notochord))
plotEmbedding(ArchRProj = notochord, colorBy = "GeneScoreMatrix", name = "sox9b", 
              embedding = "Tile_LSI_LM1.batch.UMAP",imputeWeights = getImputeWeights(notochord))
#
p <- plotBrowserTrack(ArchRProj = notochord, groupBy = "Tile_LSI_LM1.batch.C1",
                      geneSymbol = "noto", upstream = 10000,downstream = 10000,
                      ylim = c(0.001,0.999999))
grid::grid.newpage()
grid::grid.draw(p$noto)
#
notochord.seurat$Tile_LSI_LM1.batch.C1 <- notochord$Tile_LSI_LM1.batch.C1
notochord.seurat$time <- notochord$time
notochord.seurat$batch <- notochord$batch
#
markersGS <- getMarkerFeatures(
  ArchRProj = notochord, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Tile_LSI_LM1.batch.C1",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
#
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$C1
#
notochord.seurat@active.ident <- as.factor(notochord.seurat$Tile_LSI_LM1.batch.C1)
DimPlot(notochord.seurat)
table(notochord.seurat@active.ident)
notochord.seurat <- RenameIdents(notochord.seurat,"C1"="zNotoC3")
notochord.seurat <- RenameIdents(notochord.seurat,"C2"="zNotoC4")
notochord.seurat <- RenameIdents(notochord.seurat,"C3"="zNotoC4")
notochord.seurat <- RenameIdents(notochord.seurat,"C4"="zNotoC4")
notochord.seurat <- RenameIdents(notochord.seurat,"C5"="zNotoC2")
notochord.seurat <- RenameIdents(notochord.seurat,"C7"="zNotoC2")
notochord.seurat <- RenameIdents(notochord.seurat,"C8"="zNotoC1")
notochord.seurat <- RenameIdents(notochord.seurat,"C9"="zNotoC1")
#
#"zNotoC1"="Notochord_progenitor"
#"zNotoC2"="Notochord_intermediate"
#"zNotoC3"="Notochord_ngs+"
#"zNotoC4"="Notochord_col2a1a+"
#
levels(notochord.seurat) <- c("zNotoC1", "zNotoC2", "zNotoC3", "zNotoC4")
DimPlot(notochord.seurat,label = T,cols = c("#FBB4AE","#00BFFF","#CC99FF","#48D1CC"))+NoLegend()
#
notochord.seurat$ID <- as.character(notochord.seurat@active.ident)
notochord$ID <- as.character(notochord.seurat@active.ident)
#
# Peak calling
addArchRThreads(threads = 1) 
notochord <- addGroupCoverages(ArchRProj = notochord, groupBy = "ID")
#
notochord <- addReproduciblePeakSet(
  ArchRProj = notochord, 
  groupBy = "ID",
  genomeSize = 1.5e9,
  excludeChr = c("MT",paste("chr",chr_rm$V1,sep = "")), 
  genomeAnnotation = genomeAnnotation,
  geneAnnotation = GRCz11.geneAnnotation,
  pathToMacs2 = "/home/sunkeyong/.local/bin/macs2", force=T
)
#
notochord.peakset <- getPeakSet(notochord)
#
notochord <- addPeakMatrix(ArchRProj = notochord, force =T )
# Motif annotation
# loading zebrafish motif collected from DANIO-CODE
load("/home/sunkeyong/ZEPA/Custom_danRer11_Motif/danRer11.Motif.RData")
#
notochord <- addMotifAnnotations(ArchRProj = notochord, motifPWMs=danRer11.Motif, name = "Motif")
#
notochord <- addBgdPeaks(notochord)
#
notochord <- addDeviationsMatrix(
  ArchRProj = notochord, 
  peakAnnotation = "Motif",
  force = TRUE
)
#
getAvailableMatrices(notochord)
#
plotVarDev <- getVarDeviations(notochord, name = "MotifMatrix", plot = TRUE)
plotVarDev
###
markersPeaks <- getMarkerFeatures(
  ArchRProj = notochord, 
  useMatrix = "PeakMatrix", 
  groupBy = "ID",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
#
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList$zNotoC1
markerList$zNotoC2
markerList$zNotoC3
markerList$zNotoC4
#
heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1",
  #transpose = TRUE
)
draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
#
enrichMotifs2 <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = notochord,
  peakAnnotation = "Motifraw",
  cutOff = "FDR <= 0.01 & Log2FC >= 1"
)
enrichMotifs2
#
heatmapEM2 <- plotEnrichHeatmap(enrichMotifs2, n = 37, transpose = TRUE)
heatmapEM2 <- plotEnrichHeatmap(enrichMotifs2, n = 37, transpose = F)
ComplexHeatmap::draw(heatmapEM2, heatmap_legend_side = "bot", annotation_legend_side = "bot")
#
heatmapEM2 <- plotEnrichHeatmap(enrichMotifs2, n = 100, transpose = F,returnMatrix = FALSE,cutOff = 1)
ComplexHeatmap::draw(heatmapEM2, heatmap_legend_side = "bot", annotation_legend_side = "bot")
#
heatmapEM3 <- plotEnrichHeatmap(enrichMotifs2, n = 100, transpose = F,returnMatrix = T,cutOff = 1)
heatmapEM3
dim(heatmapEM3)# 224
pheatmap(heatmapEM3,cluster_rows = F,cluster_cols = F)
heatmapEM4 <- heatmapEM3
pheatmap(heatmapEM4[112:224,],cluster_rows = F,cluster_cols = F)
pheatmap(t(heatmapEM4[112:224,]),cluster_rows = F,cluster_cols = F)
#
heatmapEM5 <- heatmapEM4[112:224,3:4]
pheatmap(heatmapEM5,cluster_rows = T,cluster_cols = F)
dim(heatmapEM5)
pheatmap(heatmapEM5[80:127,],cluster_rows = F,cluster_cols = F)
pheatmap(heatmapEM5[order(heatmapEM5[,2]),],cluster_rows = F,cluster_cols = F)
pheatmap(heatmapEM5[rev(order(heatmapEM5[,1])),],cluster_rows = F,cluster_cols = F)
pheatmap(t(heatmapEM5[rev(order(heatmapEM5[,1])),]),cluster_rows = F,cluster_cols = F,color = paletteContinuous(set = "comet", n = 100))
#
write.csv(heatmapEM5,"heatmapEM5.csv")
#
notochord <- addPeak2GeneLinks(
  ArchRProj = notochord,
  reducedDims = "Tile_LSI_LM1.batch",
  useMatrix = "GeneScoreMatrix",
  dimsToUse = 1:100,
  scaleDims = NULL,
  corCutOff = 0.45,
)
p2g <- getPeak2GeneLinks(
  ArchRProj = notochord,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)
p2g
#
Notochrod.ID.color <- c("#FBB4AE","#00BFFF","#CC99FF","#48D1CC")
#
Notochrod.ID.color.p2g <- paletteDiscrete(notochord$IDx)
Notochrod.ID.color.p2g[1] <- "#FBB4AE"
Notochrod.ID.color.p2g[2] <- "#00BFFF"
Notochrod.ID.color.p2g[3] <- "#CC99FF"
Notochrod.ID.color.p2g[4] <- "#48D1CC"
Notochrod.ID.color.p2g
#
p2g <- plotPeak2GeneHeatmap(ArchRProj = notochord, groupBy = "ID",
                            k = 5,palGroup =Notochrod.ID.color.p2g,
                            corCutOff = 0.7,
                            varCutOffATAC = 0.35,
                            varCutOffRNA = 0.35)
p2g
#
p2geneDF <- metadata(notochord@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
p2geneDF
p2geneDF.tmp <- as.data.frame(p2geneDF)
p2geneDF.tmp
p2geneDF.tmp <- p2geneDF.tmp[which(p2geneDF.tmp$Correlation > 0.7),]
p2geneDF.tmp <- p2geneDF.tmp[which(p2geneDF.tmp$VarQATAC > 0.35),]
p2geneDF.tmp <- p2geneDF.tmp[which(p2geneDF.tmp$VarQRNA > 0.35),]
p2geneDF.tmp
dim(p2geneDF.tmp)
write.table(p2geneDF.tmp,file = "/home/sunkeyong/ZEPA/notochord/notochord.p2geneDF.csv",quote=F, sep = "\t",row.names = F,col.names = T)
#
hist(table(p2geneDF.tmp$geneName))
median(table(p2geneDF.tmp$geneName))
hist(table(p2geneDF.tmp$peakName),breaks=50)
median(table(p2geneDF.tmp$peakName))
#
p2g.peak.num <- as.data.frame(table(p2geneDF.tmp$geneName))
p2g.peak.num[which(p2g.peak.num$Freq > 50 ),2] <- 50
max(p2g.peak.num$Freq)
hist(p2g.peak.num$Freq,breaks = 50)
#
p2g.gene.num <- as.data.frame(table(p2geneDF.tmp$peakName))
p2g.gene.num[which(p2g.gene.num$Freq > 10 ),2] <- 10
hist(p2g.gene.num$Freq,breaks = 20)
#
barplot(table(p2g.gene.num$Freq),space=0)
barplot(table(p2g.peak.num$Freq),space=0)
# 
p <- plotBrowserTrack(ArchRProj = notochord, groupBy = "ID",tileSize = 100,
                      loops = getPeak2GeneLinks(notochord),
                      geneSymbol = "sox9b",pal=c("#FBB4AE","#00BFFF","#CC99FF","#48D1CC"),
                      ylim =c(0.1,1), upstream = 12000,downstream = 28000)
grid::grid.newpage()
grid::grid.draw(p$sox9b)
plotPDF(plotList = p, name = "sox9b.pdf", ArchRProj = notochord, addDOC = FALSE, width = 5, height = 8)
#
p <- plotBrowserTrack(ArchRProj = notochord, groupBy = "ID",tileSize = 100,
                      loops = getPeak2GeneLinks(notochord),
                      geneSymbol = "ngs",
                      pal=c("#FBB4AE","#00BFFF","#CC99FF","#48D1CC"),ylim =c(0.1,1), 
                      upstream = 15000,downstream = 20000)
grid::grid.newpage()
grid::grid.draw(p$ngs)
plotPDF(plotList = p, name = "ngs.pdf", ArchRProj = notochord, addDOC = FALSE, width = 5, height = 8)
#

markerTest <- getMarkerFeatures(
  ArchRProj = notochrod, 
  useMatrix = "PeakMatrix",
  groupBy = "ID",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "zNotoC3",
  bgdGroups = "zNotoC4"
)
#
markerTest@metadata
#
markerList <- getMarkers(markerTest, cutOff = "FDR <= 0.001 & Log2FC >= 1")
markerList$zNotoC3
x.p2 <- as.data.frame(cbind(as.character(markerList$zNotoC3$seqnames),as.character(markerList$zNotoC3$start),as.character(markerList$zNotoC3$end)))
write.table(x.p2,file = "/home/sunkeyong/ZEPA/notochord/datasave/C3_UP_C4_down.csv",quote=F, sep = "\t",row.names = F,col.names = F)
#
markerList <- getMarkers(markerTest, cutOff = "FDR <= 0.001 & Log2FC <= -1")
markerList$zNotoC3
x.p2 <- as.data.frame(cbind(as.character(markerList$zNotoC3$seqnames),as.character(markerList$zNotoC3$start),as.character(markerList$zNotoC3$end)))
write.table(x.p2,file = "/home/sunkeyong/ZEPA/notochord/datasave/C3_Down_C4_up.csv",quote=F, sep = "\t",row.names = F,col.names = F)
#
#
pma <- markerPlot(seMarker = markerTest, name = "zNotoC3", cutOff = "FDR <= 0.001 & abs(Log2FC) >= 1", plotAs = "MA")
pma
#
pma.metadata <- as.data.frame(pma$data)
ggplot(pma.metadata[seq(1,nrow(pma.metadata),100),],aes(x=x,y=y,color=color))+geom_point(shape=19,size=3,alpha = 0.5)
pma.metadata2 <- pma.metadata[which(pma.metadata$x < 1.6),]
ggplot(pma.metadata2[seq(1,nrow(pma.metadata),10),],aes(x=x,y=y,color=color))+
  geom_point(shape=19,size=1.5,alpha = 1)+ 
  scale_color_manual(values=c("#3772C7", "#D3D3D3", "#BD362F"))
#
p <- ggplot(pma.metadata2,aes(x=x,y=y,color=color))+
  geom_point(shape=19,size=0.00001,alpha = 1)+ 
  scale_color_manual(values=c("#3772C7", "#D3D3D3", "#BD362F"))+
  theme_test()
p
ggsave("Notochord.PMA3.png",
       plot=p,device="png",dpi=400,units = "cm",width = 13,height = 9)
ggsave("Notochord.PMA4.png",
       plot=p,device="png",dpi=400,units = "cm",width = 10,height = 6.5)
#
#
notochord.seurat@active.ident <- as.factor(notochord.seurat$ID)
DimPlot(notochord.seurat,pt.size =0.0001,cols =Notochrod.ID.color)
p <- DimPlot(notochord.seurat,pt.size =0.01,cols =Notochrod.ID.color)+NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/notochord/UMAP/Notochord.ID.png",plot=p,device="png",dpi=300,units = "cm",width = 12,height = 12)
##
notochord.seurat@active.ident <- as.factor(notochord.seurat$time)
time.levels <- c("7hpf","8hpf","9hpf","10hpf",
                 "11hpf","12hpf","14hpf",
                 "18hpf","20hpf","22hpf","24hpf",
                 "30hpf","34hpf","38hpf",
                 "42hpf","48hpf","72hpf")
notochord.timeColor <- c("#AD3417","#F7C7D2","#EA66A1","#E50C84",
                         "#D2E2F2","#227DBD","#004695",
                         "#CFE6C0","#A8D392","#0C9E4A","#006D32",
                         "#C6C4E0","#857BB8","#694799",
                         "#EFBF00","#B78510","#935913")
table(notochord.seurat@active.ident)
levels(notochord.seurat) <- time.levels
#
DimPlot(notochord.seurat,pt.size =0.0001,cols =notochord.timeColor)
p <- DimPlot(notochord.seurat,pt.size =0.01,cols =notochord.timeColor)+NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/notochord/UMAP/Notochord.timepoint.png",plot=p,device="png",dpi=300,units = "cm",width = 12,height = 12)
##
p <- plotEmbedding(ArchRProj = notochord, colorBy = "GeneScoreMatrix", name = "noto", pal =paletteContinuous("solarExtra"),
                   embedding = "Tile_LSI_LM1.batch.UMAP",imputeWeights = getImputeWeights(notochord))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/notochord/UMAP/notochord.noto.genescore.png",plot=p,device="png",dpi=300,units = "cm",width = 15,height = 15)
p <- plotEmbedding(ArchRProj = notochord, colorBy = "GeneScoreMatrix", name = "col2a1a", pal =paletteContinuous("solarExtra"),
                   embedding = "Tile_LSI_LM1.batch.UMAP",imputeWeights = getImputeWeights(notochord))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/notochord/UMAP/notochord.col2a1a.genescore.png",plot=p,device="png",dpi=300,units = "cm",width = 15,height = 15)
p <- plotEmbedding(ArchRProj = notochord, colorBy = "GeneScoreMatrix", name = "col11a2", pal =paletteContinuous("solarExtra"),
                   embedding = "Tile_LSI_LM1.batch.UMAP",imputeWeights = getImputeWeights(notochord))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/notochord/UMAP/notochord.col11a2.genescore.png",plot=p,device="png",dpi=300,units = "cm",width = 15,height = 15)
p <- plotEmbedding(ArchRProj = notochord, colorBy = "GeneScoreMatrix", name = "sox19a", pal =paletteContinuous("solarExtra"),
                   embedding = "Tile_LSI_LM1.batch.UMAP",imputeWeights = getImputeWeights(notochord))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/notochord/UMAP/notochord.sox19a.genescore.png",plot=p,device="png",dpi=300,units = "cm",width = 15,height = 15)
p <- plotEmbedding(ArchRProj = notochord, colorBy = "GeneScoreMatrix", name = "pacsin3", pal =paletteContinuous("solarExtra"),
                   embedding = "Tile_LSI_LM1.batch.UMAP",imputeWeights = getImputeWeights(notochord))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/notochord/UMAP/notochord.pacsin3.genescore.png",plot=p,device="png",dpi=300,units = "cm",width = 15,height = 15)
p <- plotEmbedding(ArchRProj = notochord, colorBy = "GeneScoreMatrix", name = "ngs", pal =paletteContinuous("solarExtra"),
                   embedding = "Tile_LSI_LM1.batch.UMAP",imputeWeights = getImputeWeights(notochord))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/notochord/UMAP/notochord.ngs.genescore.png",plot=p,device="png",dpi=300,units = "cm",width = 15,height = 15)
p <- plotEmbedding(ArchRProj = notochord, colorBy = "GeneScoreMatrix", name = "cnmd", pal =paletteContinuous("solarExtra"),
                   embedding = "Tile_LSI_LM1.batch.UMAP",imputeWeights = getImputeWeights(notochord))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/notochord/UMAP/notochord.cnmd.genescore.png",plot=p,device="png",dpi=300,units = "cm",width = 15,height = 15)
p <- plotEmbedding(ArchRProj = notochord, colorBy = "GeneScoreMatrix", name = "sox9a", pal =paletteContinuous("solarExtra"),
                   embedding = "Tile_LSI_LM1.batch.UMAP",imputeWeights = getImputeWeights(notochord))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/notochord/UMAP/notochord.sox9a.genescore.png",plot=p,device="png",dpi=300,units = "cm",width = 15,height = 15)
p <- plotEmbedding(ArchRProj = notochord, colorBy = "GeneScoreMatrix", name = "sox9b", pal =paletteContinuous("solarExtra"),
                   embedding = "Tile_LSI_LM1.batch.UMAP",imputeWeights = getImputeWeights(notochord))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/notochord/UMAP/notochord.sox9b.genescore.png",plot=p,device="png",dpi=300,units = "cm",width = 15,height = 15)
### for source data
pma.metadata2
colnames(pma.metadata2) <- c("Log2Accessiblity","Log2FC","color")
head(pma.metadata2)
write.csv(pma.metadata2,file = "/home/sunkeyong/ZEPA/SourceData/Fig7K.csv")
#
notochord.seurat.meta <- cbind(notochord.seurat$time,
                               as.character(notochord.seurat$ID),
                               notochord.seurat@reductions$UMAP@cell.embeddings[,1],
                               notochord.seurat@reductions$UMAP@cell.embeddings[,2])
head(notochord.seurat.meta)
colnames(notochord.seurat.meta) <- c("time","cell type","UMAP1","UMAP2")
rownames(notochord.seurat.meta) <- rownames(notochord.seurat.meta)
write.csv(notochord.seurat.meta,file = "/home/sunkeyong/ZEPA/SourceData/Fig7H_I.csv")
#
write.csv(table(p2g.gene.num$Freq),file = "/home/sunkeyong/ZEPA/SourceData/FigS11H_Left.csv")
write.csv(table(p2g.peak.num$Freq),file = "/home/sunkeyong/ZEPA/SourceData/FigS11H_Right.csv")
#
p2g.M <- plotPeak2GeneHeatmap(ArchRProj = notochord, 
                              groupBy = "ID",
                              k = 5,
                              corCutOff = 0.7,
                              returnMatrices =T,
                              varCutOffATAC = 0.35,
                              varCutOffRNA = 0.35)
p2g.M$ATAC$matrix

#
write.csv(p2g.M$ATAC$matrix,file = "/home/sunkeyong/ZEPA/SourceData/FigS11G_Left.csv")
write.csv(p2g.M$RNA$matrix,file = "/home/sunkeyong/ZEPA/SourceData/FigS11G_Right.csv")
#








notochrod.ID.color <- c("#FBB4AE","#00BFFF","#CC99FF","#48D1CC")
notochord.seurat@active.ident <- as.factor(notochord.seurat$ID)
DimPlot(notochord.seurat,reduction = "UMAP",cols = notochrod.ID.color,label = T)
p <- DimPlot(notochord.seurat,reduction = "UMAP",cols = notochrod.ID.color,
             label = F,pt.size =0.0001)+NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/notochord/fig_save/notochord_ID.png",plot=p,device="png",dpi=300,units = "cm",width = 14,height = 14)
#
##
DimPlot(notochord.seurat,reduction = "UMAP",cols = c("#CFE6C0","#A8D392","#0C9E4A","#006D32",
                                                     "#C6C4E0","#857BB8","#694799",
                                                     "#EFBF00","#B78510","#935913"),label = T)
p <- DimPlot(notochord.seurat,reduction = "UMAP",cols = c("#CFE6C0","#A8D392","#0C9E4A","#006D32",
                                                          "#C6C4E0","#857BB8","#694799",
                                                          "#EFBF00","#B78510","#935913"),
             label = F,pt.size =0.0001)+NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+  
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("/home/sunkeyong/ZEPA/notochord/fig_save/notochord_time.png",plot=p,device="png",dpi=300,units = "cm",width = 14,height = 14)
#