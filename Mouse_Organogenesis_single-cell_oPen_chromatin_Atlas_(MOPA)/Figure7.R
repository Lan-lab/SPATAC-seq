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
library(ComplexHeatmap)
library("RColorBrewer")
library("circlize")
#
addArchRThreads(threads = 30) 
#
## load all arrow files
fetal.ArrowFiles1 <- c("/home/sunkeyong/MOPA_project/fragment_total/E1A.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E1B.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E1C.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E1D.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E1E.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E1F.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E1G.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E1H.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E2A.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E2B.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E2C.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E2D.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E2E.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E2F.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E2G.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E2H.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E7A.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E7B.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E7C.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E7D.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E7E.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E7F.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E7G.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E7H.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E8A.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E8B.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E8C.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E8D.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E8E.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E8F.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E8G.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E8H.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E11A.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E11B.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E11C.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E11D.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E11E.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E11F.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E11G.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E11H.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E17A.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E17B.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E17C.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E17D.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E17E.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E17F.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E17G.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E17H.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E18A.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E18B.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E18C.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E18D.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E18E.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E18F.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E18G.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E18H.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E19A.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E19B.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E19C.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E19D.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E19E.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E19F.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E19G.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E19H.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E153.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E154.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E155.arrow",
                       "/home/sunkeyong/MOPA_project/fragment_total/E156.arrow")
ArrowFiles.adult  <- c("/home/sunkeyong/MAE/arrow.file/ArchROutput/ArrowFiles/BoneMarrow_62016.arrow",
                       "/home/sunkeyong/MAE/arrow.file/ArchROutput/ArrowFiles/BoneMarrow_62216.arrow",
                       "/home/sunkeyong/MAE/arrow.file/ArchROutput/ArrowFiles/Cerebellum_62216.arrow",
                       "/home/sunkeyong/MAE/arrow.file/ArchROutput/ArrowFiles/HeartA_62816.arrow",
                       "/home/sunkeyong/MAE/arrow.file/ArchROutput/ArrowFiles/Kidney_62016.arrow",
                       "/home/sunkeyong/MAE/arrow.file/ArchROutput/ArrowFiles/LargeIntestineA_62816.arrow",
                       "/home/sunkeyong/MAE/arrow.file/ArchROutput/ArrowFiles/LargeIntestineB_62816.arrow",
                       "/home/sunkeyong/MAE/arrow.file/ArchROutput/ArrowFiles/Liver_62016.arrow",
                       "/home/sunkeyong/MAE/arrow.file/ArchROutput/ArrowFiles/Lung1_62216.arrow",
                       "/home/sunkeyong/MAE/arrow.file/ArchROutput/ArrowFiles/Lung2_62216.arrow",
                       "/home/sunkeyong/MAE/arrow.file/ArchROutput/ArrowFiles/PreFrontalCortex_62216.arrow",
                       "/home/sunkeyong/MAE/arrow.file/ArchROutput/ArrowFiles/SmallIntestine_62816.arrow",
                       "/home/sunkeyong/MAE/arrow.file/ArchROutput/ArrowFiles/Spleen_62016.arrow",
                       "/home/sunkeyong/MAE/arrow.file/ArchROutput/ArrowFiles/Testes_62016.arrow",
                       "/home/sunkeyong/MAE/arrow.file/ArchROutput/ArrowFiles/Thymus_62016.arrow",
                       "/home/sunkeyong/MAE/arrow.file/ArchROutput/ArrowFiles/WholeBrainA_62216.arrow",
                       "/home/sunkeyong/MAE/arrow.file/ArchROutput/ArrowFiles/WholeBrainA_62816.arrow")
ArrowFiles <- c(ArrowFiles1,ArrowFiles.adult)
#
## Create ArchR object
Fetal.Adult.ArchROBO <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  copyArrows = TRUE)
#
load("Merge.cluster.RData") # Anotation
head(merge.cluster)
rownames(merge.cluster)  <- merge.cluster$cell
Stage <- c(rep("Fetal",357171),rep("Adult",76174))
merge.cluster$stage <- Stage
#
merge.cluste.F <- merge.cluster[Fetal.Adult.ArchROBO$cellNames,]
#
Fetal.Adult.ArchROBO$merge.cluster <- merge.cluste.F$merge.cluster
Fetal.Adult.ArchROBO$stage <- merge.cluste.F$stage
table(Fetal.Adult.ArchROBO$stage)
##
# Create Peak Matrix
load("Fetal.Adult.peakset.RData")
Fetal.Adult.ArchROBO <- addPeakSet(Fetal.Adult.ArchROBO,peakSet = Fetal.Adult.peakset)
Fetal.Adult.ArchROBO <- addPeakMatrix(Fetal.Adult.ArchROBO,force = T)
#
Fetal.Adult.ArchROBO <- addIterativeLSI(
  ArchRProj = Fetal.Adult.ArchROBO,
  useMatrix = "PeakMatrix", 
  name = "Peak_IterativeLSI_M1", 
  iterations = 3, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(3,4), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 100000, 
  dimsToUse = 1:100,
  LSIMethod = 1
)
#
# Correct batch effect between fetal and adult dataset
# Using Harmony
# using 2:100 dims
Fetal.Adult.ArchROBO <- addHarmony(
  ArchRProj = Fetal.Adult.ArchROBO,
  reducedDims = "Peak_IterativeLSI_M1",
  name = "Peak_Harmony_M1",
  groupBy = "fetal_or_adult", #######
  dimsToUse = 2:100
)
#
Fetal.Adult.ArchROBO <- addUMAP(
  ArchRProj = Fetal.Adult.ArchROBO, 
  reducedDims = "Peak_Harmony_M1", 
  name = "Peak_Harmony_M1_UMAP", 
  nNeighbors = 100, 
  dimsToUse = 1:100, 
  minDist = 0.4, 
  metric = "cosine",force = T
)
Fetal.Adult.ArchROBO <- addUMAP(
  ArchRProj = Fetal.Adult.ArchROBO, 
  reducedDims = "Peak_IterativeLSI_M1", 
  name = "Peak_IterativeLSI_M1_UMAP", 
  nNeighbors = 100, 
  dimsToUse = 1:100, 
  minDist = 0.4, 
  metric = "cosine",force = T
)
# Plot some marker genes
p <- plotEmbedding(ArchRProj = Fetal.Adult.ArchROBO, colorBy = "GeneScoreMatrix", name = "Epcam",pal = paletteContinuous("solarExtra"), quantCut = c(0.1, 0.99), embedding = "Peak_Harmony.re1_UMAP1",imputeWeights = getImputeWeights(projHeme2))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
ggsave("featureplot/Epcam.genescore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
#
# Create pseudo-Seurat object for rapid processing metadata
#
x.num <- c(rep(1,800000),rep(0,800000),rep(2,800000))
head(x.num)
x.num.sample <- cbind(sample(x.num,length(Fetal.Adult.ArchROBO$cellNames)),
                      sample(x.num,length(Fetal.Adult.ArchROBO$cellNames)),
                      sample(x.num,length(Fetal.Adult.ArchROBO$cellNames)),
                      sample(x.num,length(Fetal.Adult.ArchROBO$cellNames)),
                      sample(x.num,length(Fetal.Adult.ArchROBO$cellNames)))
head(x.num.sample)
colnames(x.num.sample) <- c("F1","F2","F3","F4","F5")
rownames(x.num.sample) <- z.Notochordx$cellNames
#
Fetal.Adult.seuratOBO <- CreateSeuratObject(counts = t(x.num.sample), assay = "ATAC", project = "Fetal.Adult", min.cells = 0, min.features = 0)
Fetal.Adult.seuratOBO
Fetal.Adult.seuratOBO <- NormalizeData(Fetal.Adult.seuratOBO, normalization.method = "LogNormalize", scale.factor = 10000)
Fetal.Adult.seuratOBO <- FindVariableFeatures(Fetal.Adult.seuratOBO, selection.method = "vst", nfeatures = 4)
Fetal.Adult.seuratOBO <- ScaleData(Fetal.Adult.seuratOBO, features = rownames(Fetal.Adult.seuratOBO))
variablegene <- VariableFeatures(object = Fetal.Adult.seuratOBO)
Fetal.Adult.seuratOBO <- RunPCA(Fetal.Adult.seuratOBO, features = variablegene,npcs =2,reduction.name = "UMAP",reduction.key = "UMAP_")
Fetal.Adult.seuratOBO
DimPlot(Fetal.Adult.seuratOBO)+NoLegend()
#
x <- Fetal.Adult.ArchROBO@embeddings$Peak_Harmony_M1_UMAP
y <- x$df
s.UMAP <- cbind(y$`Peak_Harmony_M1#UMAP_Dimension_1`,
                y$`Peak_Harmony_M1#UMAP_Dimension_2`)
#
head(s.UMAP )
rownames(s.UMAP ) <- rownames(y)
head(s.UMAP )
colnames(s.UMAP ) <- c("UMAP_1","UMAP_2")
head(s.UMAP )
class(s.UMAP )
#s.UMAP <- as.data.frame(s.UMAP)
Fetal.Adult.seuratOBO@reductions$UMAP@cell.embeddings <- s.UMAP 
#
DimPlot(Fetal.Adult.seuratOBO,reduction = "UMAP")+NoLegend()
#
Fetal.Adult.seuratOBO$stage <- metadatax.F$stage
Fetal.Adult.seuratOBO$merge.cluster <- metadatax.F$merge.cluster
DimPlot(Fetal.Adult.seuratOBO,group.by = "stage")
DimPlot(Fetal.Adult.seuratOBO,group.by = "merge.cluster",rep=T,label = T)+NoLegend()
DimPlot(Fetal.Adult.seuratOBO,group.by = "merge.cluster",split.by = "stage",rep=T,label = T)+NoLegend()
#
##
Fetal.Adult.seuratOBO.tmp <- Fetal.Adult.seuratOBO
Fetal.Adult.seuratOBO@active.ident <- as.factor(Fetal.Adult.seuratOBO$merge.cluster)
#
Fetal.Adult.seuratOBO.tmp@active.ident <- as.factor(Fetal.Adult.seuratOBO.tmp$stage)
Idents(Fetal.Adult.seuratOBO,cells=Cells(subset(Fetal.Adult.seuratOBO.tmp,idents="Adult"))) <- "Adult"
Fetal.Adult.seuratOBO$plot.for.fetal <- Fetal.Adult.seuratOBO@active.ident
#
Fetal.Adult.seuratOBO@active.ident <- as.factor(Fetal.Adult.seuratOBO$merge.cluster)
Idents(Fetal.Adult.seuratOBO,cells=Cells(subset(Fetal.Adult.seuratOBO.tmp,idents="Fetal"))) <- "Fetal"
Fetal.Adult.seuratOBO$plot.for.Adult <- Fetal.Adult.seuratOBO@active.ident
#
#
Fetal.Adult.seuratOBO@active.ident <- as.factor(Fetal.Adult.seuratOBO$plot.for.Adult)
levels(Fetal.Adult.seuratOBO) <- c("Fetal","Unknown","T cells-1","Hepatocytes","B cells","Ex or in neurons",
                  "Enterocytes","Cardiomyocytes","Cerebellar granule cells","Endothelial cells-1","Hematopoietic progenitors",
                  "Proximal tubule","T or NK cells","Erythroblasts","Sperm","Inhibitory neurons","B or Mac or Microglia",
                  "Mono or Mac or DC","DCT or CD or Loop of henle","Astrocytes","Type I pneumocytes","Oligodendrocytes","Endothelial cells-2",
                  "Endothelial cells-3","Monocytes","Endothelial cells-4","T cells-2","Purkinje cells","Immature B cells","Ex. neurons SCPN",
                  "Type II pneumocytes")
color.for.adult <- c("gray","#FDCDAC","#FFFF99","#B3CDE3","#6495ED","#00FFFF", ## 1--5 
                     "#FBB4AE","#F4CAE4","#ABDDA4","#B15928","#D6604D", ## 6--10
                     "#999999","#00BFFF","#CC99FF","#48D1CC","#FFED6F", ## 11--15 
                     "#9E0142","#0099FF","#B3E2CD","#276419","#66C2A5", ## 16--20
                     "#B3EE3A","#6A3D9A","#87CEFA","#ABD9E9","#B2ABD2", ## 21--25 
                     "#0000FF","#FF33CC","#66BD63","#E08214","#DE77AE")
DimPlot(Fetal.Adult.seuratOBO,repel = T,label = T,cols = color.for.adult)+NoLegend()
#
Fetal.Adult.seuratOBO@active.ident <- as.factor(Fetal.Adult.seuratOBO$plot.for.fetal)
levels(Fetal.Adult.seuratOBO) <- c("Adult","Neural Tube","Notochord cells","Neural progenitor cells", "Radial glia","Postmitotic premature neurons",
                  "Sensory neurons","Oligodendrocyte Progenitors","Premature oligodendrocyte","Schwann cell precursor",
                  "Inhibitory neuron progenitors","Inhibitory interneurons","Inhibitory neurons","Excitatory neurons",
                  "Cholinergic neurons","Granule neurons","Ependymal cell","Isthmic organizer cells","Connective tissue progenitors",
                  "Chondrocytes & osteoblasts","Chondroctye progenitors","Intermediate Mesoderm","Limb mesenchyme","Jaw and tooth progenitors",
                  "Early mesenchyme","Myocytes","Cardiac muscle lineages","Epithelial cells","Melanocytes","Hepatocytes",
                  "Endothelial cells","White blood cells","Definitive erythroid lineage","Primitive erythroid lineage")
color.for.fetal <- c("gray","#FDCDAC","#FFFF99","#B3CDE3","#6495ED","#00FFFF", ## 1--5 
                     "#FBB4AE","#F4CAE4","#ABDDA4","#B15928","#D6604D", ## 6--10
                     "#999999","#00BFFF","#CC99FF","#48D1CC","#FFED6F", ## 11--15 
                     "#9E0142","#0099FF","#B3E2CD","#276419","#66C2A5", ## 16--20
                     "#B3EE3A","#6A3D9A","#87CEFA","#ABD9E9","#B2ABD2", ## 21--25 
                     "#0000FF","#FF33CC","#66BD63","#E08214","#DE77AE", ## 26--30
                     "#9970AB","#EE82EE","#99CCCC")
DimPlot(Fetal.Adult.seuratOBO,repel = T,label = T,cols = color.for.fetal,order = order.plot,shuffle =F,
        raster= T,raster.dpi=c(1000,1000))+NoLegend()
#
Fetal.Adult.seuratOBO@active.ident <- as.factor(Fetal.Adult.seuratOBO$stage)
p <- DimPlot(Fetal.Adult.seuratOBO,repel = T, label = F, pt.size = 0.000001)+NoLegend()+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 1))+
  theme(panel.border = element_blank())
p
ggsave("Fetal.Adult.stage.coembedding.png",plot=p,device="png",dpi=400,units = "cm",width = 15,height = 15)
#
#
#
# Create dowmsample data
#
table(Fetal.Adult.seuratOBO@active.ident)
Fetal.Adult.seuratOBO@active.ident <- Fetal.Adult.seuratOBO$merge.cluster
Fetal.Adult.seuratOBO.downsample <- subset(Fetal.Adult.seuratOBO,downsample=600)
Fetal.Adult.seuratOBO.downsample
DimPlot(Fetal.Adult.seuratOBO.downsample,repel = T,label = T)+NoLegend()                         
#
projHeme2 <- Fetal.Adult.ArchROBO[Cells(Fetal.Adult.seuratOBO.downsample), ]
projHeme2
#
###
getAvailableMatrices(projHeme2)
## get peak matrix
Get.Peak.mat <- getMatrixFromProject(
  ArchRProj = projHeme2,
  useMatrix = "PeakMatrix")
#
Get.Peak.mat.matrix <- Get.Peak.mat@assays$data$PeakMatrix
dim(Get.Peak.mat.matrix)
####
peak_chr_name <- rep(as.character(Get.Peak.mat@rowRanges@seqnames@values),Get.Peak.mat@rowRanges@seqnames@lengths)
peak_start_site <- Get.Peak.mat@rowRanges@ranges@start
peak_end_site <- peak_start_site+500
peak.names <- paste(peak_chr_name,peak_start_site,peak_end_site,sep = "-")
head(peak.names)
###
rownames(Get.Peak.mat.matrix) <- peak.names
head(Get.Peak.mat.matrix)
dim(Get.Peak.mat.matrix)
#
# Create peak-martrix of downsample dataset
Peak.matrix.suerat <- CreateSeuratObject(counts = Get.Peak.mat.matrix, assay = "ATAC", project = "fetal_adult", min.cells = 0, min.features = 0)
Peak.matrix.suerat
#
#
Peak.matrix.suerat$ID <- Fetal.Adult.seuratOBO.downsample$merge.cluster
table(Peak.matrix.suerat$ID)
#
Peak.matrix.suerat <- NormalizeData(Peak.matrix.suerat, normalization.method = "RC",scale.factor = 1000000)
#
# for calculating the mean values of peak-cell type accessibility matrix
Peak.matrix.suerat@active.ident <- as.factor(Peak.matrix.suerat$ID)
tx <- as.data.frame(table(Peak.matrix.suerat@active.ident))
as.character(tx$Var1[1])
ax <- c()
axsum <- c()
i=1
for(i in 1:length(tx$Var1)){
  ax <- as.matrix(rowSums(as.matrix(subset(Peak.matrix.suerat,idents = as.character(tx$Var1[i]))@assays$ATAC@data))/tx$Freq[i])
  axsum <- cbind(axsum,ax)
}
colnames(axsum) <- as.character(tx$Var1)
#
dim(axsum)
axsum <- as.data.frame(axsum)
#
Peak.matrix.High.accessible.peaks <- names(which(Matrix::rowSums(Peak.matrix.suerat@assays$ATAC@counts) > 100))
length(Peak.matrix.High.accessible.peaks)
axsum2 <- as.data.frame(axsum[Peak.matrix.High.accessible.peaks,])
#
t.cor.log.spearman <- cor(log2(axsum2+1),method = c("spearman"))
#
#
color.set.for.heatmap <- colorRampPalette(c("#2166AC","#4393C3","#92C5DE","#D1E5F0","white","#FDDBC7","#F4A582","#D6604D","#B2182B"))(30)
pheatmap::pheatmap(t.cor.log.spearman[31:63,1:30],border_color = NA,
                   color = color.set.for.heatmap)
#
save.image("Mouse.Fetal.Adult.RData")
#