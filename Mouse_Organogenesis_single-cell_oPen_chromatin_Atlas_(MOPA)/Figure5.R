setwd("/home/sunkeyong/MOPA_revision1/Epithelial")
#
library(Seurat)
library(ArchR)
addArchRGenome("mm10")
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
addArchRThreads(threads = 13) 
library(enrichplot)
library(RColorBrewer)
library(clusterProfiler)
library(circlize)
library(ComplexHeatmap)
library(DOSE)
library(org.Mm.eg.db)
# load the annotation results from the figure 4
load("MOPA.epithelial.RData")
#
# Identification of Positive TF-Regulators
getAvailableMatrices(MOPA.epithelial)
#
seGroupMotif <- getGroupSE(ArchRProj = MOPA.epithelial, useMatrix = "MotifMatrix", groupBy = "ID.F3")
seGroupMotif
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs
corGSM_MM <- correlateMatrices(
  ArchRProj = MOPA.epithelial,
  useMatrix1 = "GeneScoreMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "PeakMatrix_LSI"
)
corGSM_MM
corGIM_MM <- correlateMatrices(
  ArchRProj = MOPA.epithelial,
  useMatrix1 = "GIM_PeakX4",
  useMatrix2 = "MotifMatrix",
  reducedDims = "PeakMatrix_LSI"
)
corGIM_MM
#
corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
corGIM_MM$maxDelta <- rowData(seZ)[match(corGIM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]
####
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
corGIM_MM <- corGIM_MM[order(abs(corGIM_MM$cor), decreasing = TRUE), ]
corGIM_MM <- corGIM_MM[which(!duplicated(gsub("\\-.*","",corGIM_MM[,"MotifMatrix_name"]))), ]
corGIM_MM$TFRegulator <- "NO"
corGIM_MM$TFRegulator[which(corGIM_MM$cor > 0.5 & corGIM_MM$padj < 0.01 & corGIM_MM$maxDelta > quantile(corGIM_MM$maxDelta, 0.75))] <- "YES"
sort(corGIM_MM[corGIM_MM$TFRegulator=="YES",1])
corGIM_MM
#
View(data.frame(corGIM_MM))
x <- data.frame(corGIM_MM)
dim(x)
table(x$TFRegulator)
hist(x$cor,breaks = 60)
#
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
##
FeaturePlot(TOME_V2_Epithelial,features = c("Nr5a2"),min.cutoff = 0)
####
# Nr5a2 PAX8 Grhl1
# Trp63 top rank
motifs <- c("Nkx")
markerMotifs <- getFeatures(MOPA.epithelial, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
markerMotifs
#
corGIM_MM$MotifMatrix_name[1:50]
#
plotEmbedding(
  ArchRProj = MOPA.epithelial, 
  colorBy = "MotifMatrix", 
  name = "z:Nkx21_399", 
  embedding = "TileMatrix_UMAP_1",
  imputeWeights = getImputeWeights(MOPA.epithelial))
########
##
#
# Motif Footprinting
motifPositions <- getPositions(MOPA.epithelial)
motifPositions
#
motifs <- c("")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs
##
addArchRThreads(threads = 20) 
#
seFoot <- getFootprints(
  ArchRProj = MOPA.epithelial, 
  positions = motifPositions[markerMotifs[751:884]], 
  groupBy = "ID.F3"
)
#
plotFootprints(
  seFoot = seFoot,
  ArchRProj = MOPA.epithelial, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5)
####
# import the motif deviation matrix
Get.motif.mat <- getMatrixFromProject(ArchRProj = MOPA.epithelial,useMatrix = "MotifMatrix")
rownames(Get.motif.mat@assays@data$deviations) <- Get.motif.mat@elementMetadata$name
Get.deviations.matrix <- Get.motif.mat@assays@data$deviations
head(Get.deviations.matrix)
dim(Get.deviations.matrix)
hist(as.numeric(Get.deviations.matrix))
rownames(Get.motif.mat@assays@data$z) <- Get.motif.mat@elementMetadata$name
Get.z.matrix <- Get.motif.mat@assays@data$z
head(Get.z.matrix)
dim(Get.z.matrix)
hist(as.numeric(Get.z.matrix))
save(Get.motif.mat,file = "motif/Get.motif.mat.RData")
save(Get.deviations.matrix,file = "motif/Get.deviations.matrix.RData")
save(Get.z.matrix,file = "motif/Get.z.matrix.RData")
#
plotVarDev <- getVarDeviations(MOPA.epithelial, name = "MotifMatrix", plot = TRUE)
plotVarDev
##
Get.deviations <- CreateSeuratObject(counts = Get.deviations.matrix, assay = "deviations", project = "z.deviations", min.cells = 0, min.features = 0)
hist(as.numeric(Get.deviations@assays$deviations@counts))
Get.deviations <- CreateSeuratObject(counts = Get.z.matrix, assay = "deviations", project = "z.deviations", min.cells = 0, min.features = 0)
hist(as.numeric(Get.deviations@assays$deviations@counts))
Get.deviations
##
MOPA.epithelial.metadata <- as.data.frame(MOPA.epithelial@cellColData)
MOPA.epithelial.metadata <- MOPA.epithelial.metadata[rownames(Get.deviations@meta.data),]
Get.deviations$ID.F3 <- MOPA.epithelial.metadata$ID.F3
#
Get.deviations@active.ident <- as.factor(Get.deviations$ID.F3)
downsample.num <- 170
Get.deviations.sub <- subset(Get.deviations, downsample = downsample.num)
Get.deviations.sub
Get.deviations.sub <- ScaleData(Get.deviations.sub, features = rownames(Get.deviations.sub))
hist(as.numeric(Get.deviations.sub@assays$deviations@counts))
#
# laod the PseudoCell function
library(Seurat)
library(purrr)
#
PseudoCell.num <- 170
MOPA.epithelial.motif <- PseudoCell(Get.deviations.sub, "deviations","data","ID.F3",PseudoCell.num)
MOPA.epithelial.motif
rownames(MOPA.epithelial.motif@meta.data)
MOPA.epithelial.motif$ID <- unlist(strsplit(rownames(MOPA.epithelial.motif@meta.data),"_Cell"))[seq(1,2*length(rownames(MOPA.epithelial.motif@meta.data)),2)]
MOPA.epithelial.motif@active.ident <- as.factor(MOPA.epithelial.motif$ID)
MOPA.epithelial.motif <- ScaleData(MOPA.epithelial.motif, features = rownames(MOPA.epithelial.motif))
#
DoHeatmap(MOPA.epithelial.motif, features = "Zfp263-877") + NoLegend()
#
library(dplyr)
MOPA.epithelial.motif.markers <- FindAllMarkers(MOPA.epithelial.motif, only.pos = TRUE, min.pct = 0.01, logfc.threshold = 0.01)
MOPA.epithelial.motif.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top10 <- MOPA.epithelial.motif.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DoHeatmap(MOPA.epithelial.motif, features = top10$gene) + NoLegend()
#
x.dev <- MOPA.epithelial.motif@assays$RNA@scale.data
rownames(x.dev) <- paste(unlist(strsplit(rownames(x.dev),"-"))[seq(1,2*dim(x.dev)[1],2)],
                         unlist(strsplit(rownames(x.dev),"-"))[seq(2,2*dim(x.dev)[1],2)],sep = "_")

dim(x.dev)
View(x.dev)
head(x.dev)
hist(x.dev)
x.dev[1:3,1:3]
#
x.dev.num <-as.numeric(x.dev)
hist(x.dev.num)
#
x.dev2 <- x.dev
x.dev2[x.dev2 < -2]<- -2
x.dev2[x.dev2 > 4]<- 4
x.dev2 <- as.matrix(x.dev2)
x.dev2.num <-as.numeric(x.dev2)
hist(x.dev2.num)
#
library(pheatmap)
library(ComplexHeatmap)
library("RColorBrewer")
library("circlize")
#col_fun = colorRamp2(c(0.7, 1, 1.3), c("#1F77B4","#FFFFFF","#D62728"))
col_fun = colorRamp2(c(-4, 0, 4), c("#1F77B4","#FFFFFF","#D62728"))
col_fun = colorRamp2(c(-2, 0, 2), c("Blue","#FFFFFF","yellow"))
#
pheatmap(x.dev)
pheatmap(x.dev2)
##
Heatmap(x.dev,
        cluster_rows = T,
        col = col_fun,
        cluster_columns = T,
        show_column_names = F,
        show_row_names = F)
#
sample_info = data.frame(CellType = sort(rep(as.character(as.data.frame(table(MOPA.epithelial.motif@active.ident))$Var1),downsample.num/PseudoCell.num)))
col.epi.27 <- c ("#32CD32", "#008B8B", "#FFE4B5", "#40E0D0", "#FF6347", # 1-5
                 "#FF1493", "#8B008B", "#FFA500", "#6A5ACD", "#DEB887", # 7-11
                 "#0000CD", "#9400D3", "#87CEEB", "#20B2AA", "#E9967A", # 12-16
                 "#800080", "#00FFFF", "#F08080", "#228B22", "#1E90FF", # 17-18-19-22-23
                 "#8B4513", "#FF00FF")
col.use.1.tmp <- col.epi.27
names(col.use.1.tmp) <- as.data.frame(table(MOPA.epithelial.motif@active.ident))$Var1
col.use.1.tmp
class(col.use.1.tmp)
col.use.1.tmp2 <- list("CellType"=col.use.1.tmp)
col.use.1.tmp2
top_color <- HeatmapAnnotation(df=sample_info,col=col.use.1.tmp2, annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"))
#
colnames(x.dev) <- sort(rep(as.character(as.data.frame(table(MOPA.epithelial.motif@active.ident))$Var1),downsample.num/PseudoCell.num))
motif.plot1 <- Heatmap(x.dev,
                       row_km=0,
                       col = col_fun,
                       cluster_rows = F,
                       cluster_columns = F,
                       show_column_names = F,
                       show_row_names = F,
                       heatmap_legend_param = list(
                         title = "Modified Expression",
                         title_position = "leftcenter-rot", 
                         #at=c(-1,0,1),
                         legend_height = unit(2, "cm")
                       ),
                       #row_title = "motif (1082)",
                       row_title_side = c("left"),
                       column_title = "CellType",
                       column_title_side = c("top"),
                       #top_annotation=top_color,
                       show_column_dend=F,
                       show_row_dend=F, row_dend_reorder = TRUE)
#
draw(motif.plot1, heatmap_legend_side = "right", show_annotation_legend = T)
x.dev.cor <- cor(x.dev)
Heatmap(x.dev.cor,
        cluster_rows = T,
        cluster_columns = T,
        top_annotation=top_color,
        show_column_dend=F,
        show_row_dend=F,
        show_column_names = F,
        show_row_names = F)
#
#
dim(x.dev)
x.dev[1:3,1:3]
x.dev[x.dev < -2]<- -2
x.dev[x.dev > 2]<- 2
#
length(which (plotVarDev$data$combinedVars > 1.6)) ## 127
length(which (plotVarDev$data$combinedVars > 1.5)) ## 163
# 
x.dev.var <- x.dev[plotVarDev$data$idx[1:163],] 
dim(x.dev.var)
x.dev.var[1:3,1:3]
#
motif.plot2 <- Heatmap(x.dev.var,
                       row_km=0,
                       col = col_fun,
                       cluster_rows = T,
                       cluster_columns = T,
                       show_column_names = T,
                       show_row_names = F,
                       heatmap_legend_param = list(
                         title = "Modified Expression",
                         title_position = "leftcenter-rot", 
                         #at=c(-1,0,1), 
                         legend_height = unit(2, "cm") 
                       ),
                       #row_title = "motif (1082)",
                       row_title_side = c("left"),
                       column_title = "CellType",
                       column_title_side = c("top"),
                       top_annotation=top_color,
                       show_column_dend=F,
                       show_row_dend=F, row_dend_reorder = TRUE)
#
draw(motif.plot2, heatmap_legend_side = "right", show_annotation_legend = F)
##
x.dev.var2 <- x.dev.var[rownames(x.dev.var)[row_order(motif.plot2)],]
rownames(x.dev.var2)
marker_TFs <- c("Lef1-734")
marker_TFs <- as.character(plotVarDev$data$name[81:120])
#marker_TFs <- subset(cor.metadata.merge,maxDelta > maxDelta_cutoff & cor > cor_cutoff)$MotifMatrix_name
marker_TFs <- intersect(rownames(x.dev.var2),unique(marker_TFs))
marker_TFs
TF_pos <- which(rownames(x.dev.var2) %in% marker_TFs)
row_anno <-  rowAnnotation(marker_TFs = anno_mark(at = TF_pos, 
                                                  labels = marker_TFs))
motif.plot3 <- Heatmap(x.dev.var2,
                       row_km=0,
                       #col = colorRamp2(seq(-127,128,1)/64, paletteContinuous("horizonExtra")),
                       col = colorRamp2(seq(-127,128,1)/64, paletteContinuous("solarExtra")),
                       #col = colorRamp2(seq(-127,128,1)/64, paletteContinuous("blueYellow")),
                       #col = colorRamp2(seq(-127,128,1)/64, colorRampPalette(c("#2166AC","#4393C3","white","#F4A582","#D6604D","#B2182B"))(256)),
                       #col = c("#01665e", "white", "red"),
                       cluster_rows = F,
                       cluster_columns = T,
                       show_column_names = T,
                       show_row_names = F,
                       heatmap_legend_param = list(
                         title = "Modified Expression",
                         title_position = "leftcenter-rot", 
                         #at=c(-1,0,1), 
                         legend_height = unit(2, "cm")
                       ),
                       #row_title = "motif (1082)",
                       row_title_side = c("left"),
                       column_title = "CellType",
                       column_title_side = c("top"),
                       top_annotation=top_color,
                       show_column_dend=F,
                       show_row_dend=F, 
                       row_dend_reorder = TRUE,
                       right_annotation = row_anno)
#
draw(motif.plot3, heatmap_legend_side = "right", show_annotation_legend = F)
##
#
dim(x.dev.var2)
x.dev.var3 <- x.dev.var2[,c("C10","C22","C14","C25","C17","C19","C23","C12",
                            "C11","C16","C8","C2","C24","C4","C3",
                            "C7","C15","C18","C9","C5","C1","C13")]
dim(x.dev.var3)
#
colnames(x.dev.var3) <- c("epi_C9","epi_C19","epi_C13","epi_C22","epi_C16","epi_C18",
                          "epi_C20","epi_C11","epi_C10","epi_C15",
                          "epi_C7","epi_C2","epi_C21","epi_C4","epi_C3","epi_C6",
                          "epi_C14","epi_C17","epi_C8","epi_C5","epi_C1","epi_C12")
#
motif.plot4 <- Heatmap(x.dev.var3,
                       row_km=0,
                       #col = colorRamp2(seq(-127,128,1)/64, paletteContinuous("horizonExtra")),
                       col = colorRamp2(seq(-127,128,1)/64, paletteContinuous("solarExtra")),
                       #col = colorRamp2(seq(-127,128,1)/64, paletteContinuous("blueYellow")),
                       #col = colorRamp2(seq(-127,128,1)/64, colorRampPalette(c("#2166AC","#4393C3","white","#F4A582","#D6604D","#B2182B"))(256)),
                       #col = c("#01665e", "white", "red"),
                       cluster_rows = F,
                       cluster_columns = F,
                       show_column_names = T,
                       show_row_names = F,
                       heatmap_legend_param = list(
                         title = "Modified Expression",
                         title_position = "leftcenter-rot", 
                         #at=c(-1,0,1),
                         legend_height = unit(2, "cm")
                       ),
                       #row_title = "motif (1082)",
                       row_title_side = c("left"),
                       column_title = "CellType",
                       column_title_side = c("top"),
                       #top_annotation=top_color,
                       show_column_dend=F,
                       show_row_dend=F, 
                       row_dend_reorder = TRUE,
                       right_annotation = row_anno)
#
draw(motif.plot4, heatmap_legend_side = "right", show_annotation_legend = F)
##
head(x.dev.var3)
write.csv(x.dev.var3,file="Mepi.motif.heatmap.csv",col.names =0)
#
#######
#######
# plot the correlation between TF deviation and gene expression for the individual TF
x.dev.cor <- x.dev
hist(x.dev.cor)
head(x.dev.cor)
View(x.dev.cor)
colnames(x.dev.cor) <- unlist(strsplit(rownames(MOPA.epithelial.motif@meta.data),"_Cell"))[seq(1,2*length(rownames(MOPA.epithelial.motif@meta.data)),2)]
#
table(epi.co.seurat$ID.F)
table(epi.co.seurat$id)
table(epi.co.seurat$ID.batch)
epi.co.seurat2 <- subset(epi.co.seurat, cells = Cells(Get.deviations.sub))
epi.co.seurat2@active.ident <- as.factor(as.character(epi.co.seurat2$ID.F))
table(epi.co.seurat2@active.ident)
####
epi.co.seurat2
epi.co.seurat2 <- ScaleData(epi.co.seurat2, features = rownames(epi.co.seurat2))
#
PseudoCell.num <- 170
epi.co.seurat2$ID.F3 <- epi.co.seurat2@active.ident
MOPA.epithelial.motif <- PseudoCell(epi.co.seurat2, "RNA","data","ID.F3",PseudoCell.num)
MOPA.epithelial.motif
rownames(MOPA.epithelial.motif@meta.data)
MOPA.epithelial.motif$ID <- unlist(strsplit(rownames(MOPA.epithelial.motif@meta.data),"_Cell"))[seq(1,2*length(rownames(MOPA.epithelial.motif@meta.data)),2)]
MOPA.epithelial.motif@active.ident <- as.factor(MOPA.epithelial.motif$ID)
MOPA.epithelial.motif <- ScaleData(MOPA.epithelial.motif, features = rownames(MOPA.epithelial.motif))
#
DoHeatmap(MOPA.epithelial.motif, features = "Nkx2-1") + NoLegend()
#
x.RNA.dev <- MOPA.epithelial.motif@assays$RNA@scale.data
dim(x.RNA.dev)
x.RNA.dev <- MOPA.epithelial.motif@assays$RNA@data
#
head(x.RNA.dev)
colnames(x.RNA.dev) <- unlist(strsplit(rownames(MOPA.epithelial.motif@meta.data),"_Cell"))[seq(1,2*length(rownames(MOPA.epithelial.motif@meta.data)),2)]
dim(x.RNA.dev)
###
head(corGSM_MM$GeneScoreMatrix_name)
head(corGSM_MM$MotifMatrix_name)
motif.name <- as.data.frame(rownames(x.dev.cor))
View(motif.name)
plot(x.RNA.dev["Nkx2-1",],x.dev.cor["Nkx21_399",][1:22])
plot(x.RNA.dev[1,])
plot(x.RNA.dev[2,])
cor.test(x.RNA.dev["Nkx2-1",],as.numeric(x.dev.cor["Nkx21_399",][1:22]))
#
#############
##############
length(unlist(strsplit(rownames(x.dev),"_"))[seq(1,2*dim(x.dev)[1],2)])
length(unique(unlist(strsplit(rownames(x.dev),"_"))[seq(1,2*dim(x.dev)[1],2)]))
#
dim(x.dev)
x.dev.cor <- x.dev
x.dev.cor[1:3,1:3]
colnames(x.dev.cor) <- unlist(strsplit(rownames(MOPA.epithelial.motif@meta.data),"_Cell"))[seq(1,2*length(rownames(MOPA.epithelial.motif@meta.data)),2)]
class(x.dev.cor)
head(x.dev.cor)
x.dev.cor <- as.data.frame(x.dev.cor)
dim(x.dev.cor)
x.dev.cor$name <- unlist(strsplit(rownames(x.dev.cor),"_"))[seq(1,2*dim(x.dev.cor)[1],2)]
head(x.dev.cor)
#
x.dev.cor1 <- x.dev.cor[!duplicated(x.dev.cor[,"name"]),]
head(x.dev.cor1)
dim(x.dev.cor1)
#
x.dev.cor2 <- x.dev.cor1
rownames(x.dev.cor2) <- x.dev.cor2$name
head(x.dev.cor2)
x.dev.cor2 <- x.dev.cor2[,-23]
dim(x.dev.cor2) ## 797
#
intersect(rownames(x.dev.cor2),rownames(x.RNA.dev)) ### 691
#
co.TF.name <- intersect(rownames(x.dev.cor2),rownames(x.RNA.dev)) 
x.dev.cor3 <- x.dev.cor2[co.TF.name,]
x.RNA.dev3 <- x.RNA.dev[co.TF.name,]
#
dim(x.dev.cor3) 
dim(x.RNA.dev3)
x.dev.cor3[1:3,1:3]
x.RNA.dev3[1:3,1:3]
#
plot(x.RNA.dev3[1,],x.dev.cor3[1,])
plot(x.RNA.dev3[2,],x.dev.cor3[2,])
plot(x.RNA.dev3[3,],x.dev.cor3[3,])
plot(x.RNA.dev3[4,],x.dev.cor3[4,])
plot(x.RNA.dev3[5,],x.dev.cor3[5,])
##
#
#
#
TF.correlation.color <- as.data.frame(cbind(c("Keratinocyte_1","Epidermal stem cells_2","Branchial arch ectodermal cells_3",
                                              "Epidermal progenitors_4","Keratinocyte_5","Otic vesicle epithelium_6",
                                              "Otic sensory epithelium_7","Pericardium_8","Intestinal epithelium_9",
                                              "Second branchial arch_10","Olfactory epithelium_11","Hair follicle stem cell_12",
                                              "Intestinal stem cells_13","Renal epithelium_14","Retina epithelium_15",
                                              "First branchial arch epithelium_16","Utricle and saccule epithelium_17",
                                              "Lung epithelium_18","Endocrine cells_19",
                                              "AER_20","Urothelium_21","Stomach epithelium_22"),
                                            c("#32CD32", "#008B8B", "#FFE4B5", "#40E0D0", "#FF6347","#FF1493", "#8B008B", "#FFA500", "#6A5ACD", "#DEB887",
                                              "#0000CD", "#9400D3", "#87CEEB", "#20B2AA", "#E9967A","#800080", "#00FFFF", "#F08080", "#228B22", "#1E90FF",
                                              "#8B4513", "#FF00FF","#9370DB", "#98FB98", "#FFFF00","#F0E68C", "red")))
rownames(TF.correlation.color) <- TF.correlation.color$V1
TF.correlation.color.rank <- TF.correlation.color[colnames(x.RNA.dev3),]
#
plot(x.RNA.dev3["Trp63",],x.dev.cor3["Trp63",],main = "Trp63",pch= 19,cex=2, col=TF.correlation.color.rank$V2,xlab="gene expression",ylab="TF Deviation score")
plot(x.RNA.dev3["Trp73",],x.dev.cor3["Trp73",],main = "Trp73",pch= 19,cex=2, col=TF.correlation.color.rank$V2,xlab="gene expression",ylab="TF Deviation score")
plot(x.RNA.dev3["Grhl1",],x.dev.cor3["Grhl1",],main = "Grhl1",pch= 19,cex=2, col=TF.correlation.color.rank$V2,xlab="gene expression",ylab="TF Deviation score")
#
#
getAvailableMatrices(MOPA.epithelial)
MOPA.epithelial <- addPeak2GeneLinks(
  ArchRProj = MOPA.epithelial,
  useMatrix = "GIM_M1_LSI.M2",
  reducedDims = "PeakMatrix_LSI.M2",
  dimsToUse = 1:100,
  corCutOff = 0.2,
  predictionCutoff = 0.2,
  maxDist = 500000,
  overlapCutoff = 0.9
)
#
p2g <- getPeak2GeneLinks(
  ArchRProj = MOPA.epithelial,
  corCutOff = 0.4,
  resolution = 1000,
  returnLoops = FALSE,
  varCutOffATAC = 0.15,
  varCutOffRNA = 0.15
)
p2g ### 73104
#
RNA.peak.num <- as.data.frame(table(p2g$idxRNA))
head(RNA.peak.num)
RNA.peak.num <- as.numeric(RNA.peak.num$Freq)
head(RNA.peak.num)
median(RNA.peak.num) ## value= 4
mean(RNA.peak.num)   ## value= 6.317864
boxplot(RNA.peak.num,outline =F)
#
RNA.peak.num[which(RNA.peak.num > 30)] <- 30
barplot(table(RNA.peak.num),space = 0,main="idxRNA",col="#BFB5D5")
#
Peak.RNA.num <- as.data.frame(table(p2g$idxATAC))
head(Peak.RNA.num)
Peak.RNA.num <- as.numeric(Peak.RNA.num$Freq)
head(Peak.RNA.num)
median(Peak.RNA.num) ## value= 1
mean(Peak.RNA.num)  ## value= 1.571217
boxplot(Peak.RNA.num,outline =F)
#
Peak.RNA.num[which(Peak.RNA.num > 6)] <- 6
barplot(table(Peak.RNA.num),space = 0,main="idxATAC",col="#AECDE0")
#
#
length(unique(p2g$idxRNA))## 11571 unique genes
length(p2g$idxRNA) #73104
#
length(unique(p2g$idxATAC)) ## 46527 unique cCREs
length(p2g$idxATAC) #73104
#
p2geneDF <- metadata(MOPA.epithelial@peakSet)$Peak2GeneLinks
p2geneDF$geneName <- mcols(metadata(p2geneDF)$geneSet)$name[p2geneDF$idxRNA]
p2geneDF$peakName <- (metadata(p2geneDF)$peakSet %>% {paste0(seqnames(.), "_", start(.), "_", end(.))})[p2geneDF$idxATAC]
p2geneDF
p2geneDF.tmp <- as.data.frame(p2geneDF)
p2geneDF.tmp
p2geneDF.tmp <- p2geneDF.tmp[which(p2geneDF.tmp$Correlation > 0.4),]
p2geneDF.tmp
p2geneDF.tmp <- p2geneDF.tmp[which(p2geneDF.tmp$VarQATAC > 0.15),]
p2geneDF.tmp
p2geneDF.tmp <- p2geneDF.tmp[which(p2geneDF.tmp$VarQRNA > 0.15),]
p2geneDF.tmp
dim(p2geneDF.tmp)
head(p2geneDF.tmp)
#
length(p2geneDF.tmp$peakName) # 73104
length(unique(p2geneDF.tmp$peakName)) # 46527
#
length(p2geneDF.tmp$geneName) # 73104
length(unique(p2geneDF.tmp$geneName)) #11571
#
#
peakset20230412 <- getPeakSet(MOPA.epithelial)
peakset20230412.metadata <- as.data.frame(rep(as.character(peakset20230412@seqnames@values),peakset20230412@seqnames@lengths))
colnames(peakset20230412.metadata) <- "chr"
peakset20230412.metadata$start <- peakset20230412@ranges@start
peakset20230412.metadata$end<- peakset20230412@ranges@start+500
peakset20230412.metadata$idx <- peakset20230412$idx
peakset20230412.metadata$peakType <- peakset20230412$peakType
peakset20230412.metadata$distToGeneStart <- peakset20230412$distToGeneStart
peakset20230412.metadata$distToTSS <- peakset20230412$distToTSS
peakset20230412.metadata$nearestGene <- peakset20230412$nearestGene
peakset20230412.metadata$N <- peakset20230412$N
peakset20230412.metadata$name <- paste(peakset20230412.metadata$chr,peakset20230412.metadata$start,peakset20230412.metadata$end,sep = "_")
rownames(peakset20230412.metadata) <- peakset20230412.metadata$name
#
head(peakset20230412.metadata)
tail(peakset20230412.metadata)
dim(peakset20230412.metadata)
#
mouseTSS <- read.table("/home/sunkeyong/softeware/10Xcellranger/refdata-cellranger-atac-mm10-1.2.0/regions/tss.bed")
head(mouseTSS)
dim(mouseTSS)
mouseTSS <- mouseTSS[!duplicated(mouseTSS[,"V4"]),]
dim(mouseTSS)
rownames(mouseTSS) <- mouseTSS$V4
#
length(intersect(mouseTSS$V4,unique(p2geneDF.tmp$geneName)))
#
mouseTSS <- mouseTSS[unique(p2geneDF.tmp$geneName),]
dim(mouseTSS)
head(mouseTSS)
tail(mouseTSS)
length(unique(rownames(mouseTSS)))
#
#
head(p2geneDF.tmp)
length(intersect(p2geneDF.tmp$geneName,mouseTSS$V4))
#
dis1 <- c()
dis2 <- c()
for (i in 1:nrow(p2geneDF.tmp)) {
  dis1 <- mouseTSS[p2geneDF.tmp[i,7],]$V2
  dis2 <- c(dis2,dis1)
}
#
p2geneDF.tmp2 <- p2geneDF.tmp
p2geneDF.tmp2$peaktoTTSS <- dis2
p2geneDF.tmp2
#
library(tidyverse)
p2geneDF.tmp2 <- p2geneDF.tmp2 %>% drop_na(peaktoTTSS)
dim(p2geneDF.tmp2)
p2geneDF.tmp2$peakcenter <- as.numeric(unlist(strsplit(p2geneDF.tmp2$peakName,"_"))[seq(2,3*70022,3)])+250
p2geneDF.tmp2$peaktoTTSSdis <- abs(p2geneDF.tmp2$peakcenter-p2geneDF.tmp2$peaktoTTSS)
#
p2geneDF.tmp2 <- p2geneDF.tmp2[which(p2geneDF.tmp2$peaktoTTSSdis < 500000),]
boxplot(p2geneDF.tmp2$peaktoTTSSdis,outline =F)
mean(p2geneDF.tmp2$peaktoTTSSdis) # 198913.8
median(p2geneDF.tmp2$peaktoTTSSdis) #176809
hist(p2geneDF.tmp2$peaktoTTSSdis)
#
head(p2geneDF.tmp2)
#
###
#####################################################################################
## Figure 5H
linkgene <- table(p2geneDF.tmp$geneName)
linkgene <- as.data.frame(linkgene)
linkgene
colnames(linkgene) <- c("genename","Freq")
rownames(linkgene) <- linkgene$genename
head(linkgene)
#
load("co_embed/TOME_M.V2.Epithelial.pseudo.Final.m1.averages.all.RData")
load("co_embed/Get.genescore.mat.matrix.suerat.averages.all.RData")
load("co_embed/TOME_M.V2.Epithelial.pseudo.Final.m1.markers.RData")
#
ATAC.cluster.averages.data2 <- Get.genescore.mat.matrix.suerat.averages.all[as.character(linkgene$genename),]
dim(ATAC.cluster.averages.data2)
head(ATAC.cluster.averages.data2)
epi.sci.pseudo.averages.data2 <- TOME_M.V2.Epithelial.pseudo.Final.m1.averages.all[as.character(linkgene$genename),]
dim(epi.sci.pseudo.averages.data2)
head(epi.sci.pseudo.averages.data2)
#
cor.values <- c()
cor.pvalues <- c()
cor.va <- c()
for (i in 1:nrow(epi.sci.pseudo.averages.data2)) {
  cor.va <- cor.test(epi.sci.pseudo.averages.data2[i,],ATAC.cluster.averages.data2[i,],method = c("spearman"))
  cor.values <- c(cor.values,cor.va$estimate)
  cor.pvalues <- c(cor.pvalues,cor.va$p.value)
}
cor.values.meta <- as.data.frame(cor.values)
head(cor.values.meta)
cor.values.meta$genename <- rownames(epi.sci.pseudo.averages.data2)
head(cor.values.meta)
rownames(cor.values.meta) <- cor.values.meta$genename
cor.values.meta$pvalues <- cor.pvalues
head(cor.values.meta)
dim(cor.values.meta)
#
cor.values.meta2 <- cor.values.meta[linkgene$genename,]
dim(cor.values.meta2)
cor.values.meta2$Freq <- linkgene$Freq
head(cor.values.meta2)
View(cor.values.meta2)
#
plot(cor.values.meta2$Freq,cor.values.meta2$cor.values)
cor.test(cor.values.meta2$Freq,cor.values.meta2$cor.values,method = c("spearman"))
fit <- lm(Freq~cor.values,data = cor.values.meta2)
fit
abline(fit,col="red") 
#
boxplot(subset(cor.values.meta2,cor.values.meta2$Freq < 5)$cor.values,
        subset(cor.values.meta2,cor.values.meta2$Freq >= 5 & cor.values.meta2$Freq < 10)$cor.values,
        subset(cor.values.meta2,cor.values.meta2$Freq >= 10 & cor.values.meta2$Freq < 15)$cor.values,
        subset(cor.values.meta2,cor.values.meta2$Freq >= 15 & cor.values.meta2$Freq < 20)$cor.values,
        subset(cor.values.meta2,cor.values.meta2$Freq >= 20 & cor.values.meta2$Freq < 25)$cor.values,
        subset(cor.values.meta2,cor.values.meta2$Freq >= 25 & cor.values.meta2$Freq < 30)$cor.values,
        subset(cor.values.meta2,cor.values.meta2$Freq >= 30)$cor.values,outline =F,
        col=brewer.pal(9,"Reds")[2:9]) # YlOrRd
#
boxplot(subset(cor.values.meta2,cor.values.meta2$Freq <= 5)$cor.values,
        subset(cor.values.meta2,cor.values.meta2$Freq > 5)$cor.values,
        subset(cor.values.meta2,cor.values.meta2$Freq > 10)$cor.values,
        subset(cor.values.meta2,cor.values.meta2$Freq > 15)$cor.values,
        subset(cor.values.meta2,cor.values.meta2$Freq > 20)$cor.values,
        subset(cor.values.meta2,cor.values.meta2$Freq > 25)$cor.values,
        subset(cor.values.meta2,cor.values.meta2$Freq > 30)$cor.values,
        outline =F,
        col=brewer.pal(9,"Reds")[2:9]) # YlOrRd

#
View(TOME_M.V2.Epithelial.pseudo.Final.m1.markers)
top500 <- TOME_M.V2.Epithelial.pseudo.Final.m1.markers %>% group_by(cluster) %>% top_n(n = 31, wt = avg_log2FC)
top500.gene <- unique(top500$gene)
#
Endocrine.sp.gene <- intersect(subset(TOME_M.V2.Epithelial.pseudo.Final.m1.markers,cluster== "Rank6.22")$gene,top500.gene)
Renal.sp.gene <- intersect(subset(TOME_M.V2.Epithelial.pseudo.Final.m1.markers,cluster== "Rank6.15")$gene,top500.gene)
Stomach.sp.gene <- intersect(subset(TOME_M.V2.Epithelial.pseudo.Final.m1.markers,cluster== "Rank6.25")$gene,top500.gene)
Lung.sp.gene <- intersect(subset(TOME_M.V2.Epithelial.pseudo.Final.m1.markers,cluster== "Rank6.19")$gene,top500.gene)
Intestinal.sp.gene <- intersect(subset(TOME_M.V2.Epithelial.pseudo.Final.m1.markers,cluster== "Rank6.10")$gene,top500.gene)
Intestinalstem.sp.gene <- intersect(subset(TOME_M.V2.Epithelial.pseudo.Final.m1.markers,cluster== "Rank6.14")$gene,top500.gene)
Pericardium.sp.gene <- intersect(subset(TOME_M.V2.Epithelial.pseudo.Final.m1.markers,cluster== "Rank6.9")$gene,top500.gene)
Urothelium.sp.gene <- intersect(subset(TOME_M.V2.Epithelial.pseudo.Final.m1.markers,cluster== "Rank6.24")$gene,top500.gene)
Oticvesicle.sp.gene <- intersect(subset(TOME_M.V2.Epithelial.pseudo.Final.m1.markers,cluster== "Rank6.7")$gene,top500.gene)
Oticsensory.sp.gene <- intersect(subset(TOME_M.V2.Epithelial.pseudo.Final.m1.markers,cluster== "Rank6.8")$gene,top500.gene)
Hair.sp.gene <- intersect(subset(TOME_M.V2.Epithelial.pseudo.Final.m1.markers,cluster== "Rank6.13")$gene,top500.gene)
Utricle.sp.gene <- intersect(subset(TOME_M.V2.Epithelial.pseudo.Final.m1.markers,cluster== "Rank6.18")$gene,top500.gene)
#
#
dim(p2geneDF.tmp)
head(p2geneDF.tmp)
#
endocrine.index1 <- c()
endocrine.index2 <- c()
for (i in 1:length(Endocrine.sp.gene)) {
  endocrine.index1 <- which(p2geneDF.tmp$geneName==Endocrine.sp.gene[i])
  endocrine.index2 <- c(endocrine.index2,endocrine.index1)
}
Renal.index1 <- c()
Renal.index2 <- c()
for (i in 1:length(Renal.sp.gene)) {
  Renal.index1 <- which(p2geneDF.tmp$geneName==Renal.sp.gene[i])
  Renal.index2 <- c(Renal.index2,Renal.index1)
}
Stomach.index1 <- c()
Stomach.index2 <- c()
for (i in 1:length(Stomach.sp.gene)) {
  Stomach.index1 <- which(p2geneDF.tmp$geneName==Stomach.sp.gene[i])
  Stomach.index2 <- c(Stomach.index2,Stomach.index1)
}
Lung.index1 <- c()
Lung.index2 <- c()
for (i in 1:length(Lung.sp.gene)) {
  Lung.index1 <- which(p2geneDF.tmp$geneName==Lung.sp.gene[i])
  Lung.index2 <- c(Lung.index2,Lung.index1)
}
Intestinal.index1 <- c()
Intestinal.index2 <- c()
for (i in 1:length(Intestinal.sp.gene)) {
  Intestinal.index1 <- which(p2geneDF.tmp$geneName==Intestinal.sp.gene[i])
  Intestinal.index2 <- c(Intestinal.index2,Intestinal.index1)
}
Intestinalstem.index1 <- c()
Intestinalstem.index2 <- c()
for (i in 1:length(Intestinalstem.sp.gene)) {
  Intestinalstem.index1 <- which(p2geneDF.tmp$geneName==Intestinalstem.sp.gene[i])
  Intestinalstem.index2 <- c(Intestinalstem.index2,Intestinalstem.index1)
}
Pericardium.index1 <- c()
Pericardium.index2 <- c()
for (i in 1:length(Pericardium.sp.gene)) {
  Pericardium.index1 <- which(p2geneDF.tmp$geneName==Pericardium.sp.gene[i])
  Pericardium.index2 <- c(Pericardium.index2,Pericardium.index1)
}
Urothelium.index1 <- c()
Urothelium.index2 <- c()
for (i in 1:length(Urothelium.sp.gene)) {
  Urothelium.index1 <- which(p2geneDF.tmp$geneName==Urothelium.sp.gene[i])
  Urothelium.index2 <- c(Urothelium.index2,Urothelium.index1)
}
Oticvesicle.index1 <- c()
Oticvesicle.index2 <- c()
for (i in 1:length(Oticvesicle.sp.gene)) {
  Oticvesicle.index1 <- which(p2geneDF.tmp$geneName==Oticvesicle.sp.gene[i])
  Oticvesicle.index2 <- c(Oticvesicle.index2,Oticvesicle.index1)
}
Oticsensory.index1 <- c()
Oticsensory.index2 <- c()
for (i in 1:length(Oticsensory.sp.gene)) {
  Oticsensory.index1 <- which(p2geneDF.tmp$geneName==Oticsensory.sp.gene[i])
  Oticsensory.index2 <- c(Oticsensory.index2,Oticsensory.index1)
}
Hair.index1 <- c()
Hair.index2 <- c()
for (i in 1:length(Hair.sp.gene)) {
  Hair.index1 <- which(p2geneDF.tmp$geneName==Hair.sp.gene[i])
  Hair.index2 <- c(Hair.index2,Hair.index1)
}
Utricle.index1 <- c()
Utricle.index2 <- c()
for (i in 1:length(Utricle.sp.gene)) {
  Utricle.index1 <- which(p2geneDF.tmp$geneName==Utricle.sp.gene[i])
  Utricle.index2 <- c(Utricle.index2,Utricle.index1)
}
#
Endocrine.sp.gene.peak <- p2geneDF.tmp[endocrine.index2,]
Renal.sp.gene.peak <- p2geneDF.tmp[Renal.index2,]
Stomach.sp.gene.peak <- p2geneDF.tmp[Stomach.index2,]
Lung.sp.gene.peak <- p2geneDF.tmp[Lung.index2,]
Intestinal.sp.gene.peak <- p2geneDF.tmp[Intestinal.index2,]
Intestinalstem.sp.gene.peak <- p2geneDF.tmp[Intestinalstem.index2,]
Pericardium.sp.gene.peak <- p2geneDF.tmp[Pericardium.index2,]
Urothelium.sp.gene.peak <- p2geneDF.tmp[Urothelium.index2,]
Oticvesicle.sp.gene.peak <- p2geneDF.tmp[Oticvesicle.index2,]
Oticsensory.sp.gene.peak <- p2geneDF.tmp[Oticsensory.index2,]
Hair.sp.gene.peak <- p2geneDF.tmp[Hair.index2,]
Utricle.sp.gene.peak <- p2geneDF.tmp[Utricle.index2,]
#
#
dim(Endocrine.sp.gene.peak)
dim(Renal.sp.gene.peak)
dim(Stomach.sp.gene.peak)
#
x.p1 <- unlist(strsplit(Endocrine.sp.gene.peak$peakName,"_"))
x.p2 <- as.data.frame(cbind(x.p1[seq(1,length(x.p1),3)],x.p1[seq(2,length(x.p1),3)],x.p1[seq(3,length(x.p1),3)]))
write.table(x.p2,file = "p2g.cluster.homer/Endocrine.forhomer.csv",quote=F, sep = "\t",row.names = F,col.names = F)
x.p1 <- unlist(strsplit(Renal.sp.gene.peak$peakName,"_"))
x.p2 <- as.data.frame(cbind(x.p1[seq(1,length(x.p1),3)],x.p1[seq(2,length(x.p1),3)],x.p1[seq(3,length(x.p1),3)]))
write.table(x.p2,file = "p2g.cluster.homer/Renal.forhomer.csv",quote=F, sep = "\t",row.names = F,col.names = F)
x.p1 <- unlist(strsplit(Stomach.sp.gene.peak$peakName,"_"));x.p2 <- as.data.frame(cbind(x.p1[seq(1,length(x.p1),3)],x.p1[seq(2,length(x.p1),3)],x.p1[seq(3,length(x.p1),3)]))
write.table(x.p2,file = "p2g.cluster.homer/Stomach.forhomer.csv",quote=F, sep = "\t",row.names = F,col.names = F)
#
x.p1 <- unlist(strsplit(Lung.sp.gene.peak$peakName,"_"));x.p2 <- as.data.frame(cbind(x.p1[seq(1,length(x.p1),3)],x.p1[seq(2,length(x.p1),3)],x.p1[seq(3,length(x.p1),3)]))
write.table(x.p2,file = "p2g.cluster.homer/Lung.forhomer.csv",quote=F, sep = "\t",row.names = F,col.names = F)
x.p1 <- unlist(strsplit(Intestinal.sp.gene.peak$peakName,"_"));x.p2 <- as.data.frame(cbind(x.p1[seq(1,length(x.p1),3)],x.p1[seq(2,length(x.p1),3)],x.p1[seq(3,length(x.p1),3)]))
write.table(x.p2,file = "p2g.cluster.homer/Intestinal.forhomer.csv",quote=F, sep = "\t",row.names = F,col.names = F)
x.p1 <- unlist(strsplit(Intestinalstem.sp.gene.peak$peakName,"_"));x.p2 <- as.data.frame(cbind(x.p1[seq(1,length(x.p1),3)],x.p1[seq(2,length(x.p1),3)],x.p1[seq(3,length(x.p1),3)]))
write.table(x.p2,file = "p2g.cluster.homer/Intestinalstem.forhomer.csv",quote=F, sep = "\t",row.names = F,col.names = F)
x.p1 <- unlist(strsplit(Pericardium.sp.gene.peak$peakName,"_"));x.p2 <- as.data.frame(cbind(x.p1[seq(1,length(x.p1),3)],x.p1[seq(2,length(x.p1),3)],x.p1[seq(3,length(x.p1),3)]))
write.table(x.p2,file = "p2g.cluster.homer/Pericardium.forhomer.csv",quote=F, sep = "\t",row.names = F,col.names = F)
x.p1 <- unlist(strsplit(Urothelium.sp.gene.peak$peakName,"_"));x.p2 <- as.data.frame(cbind(x.p1[seq(1,length(x.p1),3)],x.p1[seq(2,length(x.p1),3)],x.p1[seq(3,length(x.p1),3)]))
write.table(x.p2,file = "p2g.cluster.homer/Urothelium.forhomer.csv",quote=F, sep = "\t",row.names = F,col.names = F)
x.p1 <- unlist(strsplit(Oticvesicle.sp.gene.peak$peakName,"_"));x.p2 <- as.data.frame(cbind(x.p1[seq(1,length(x.p1),3)],x.p1[seq(2,length(x.p1),3)],x.p1[seq(3,length(x.p1),3)]))
write.table(x.p2,file = "p2g.cluster.homer/Oticvesicle.forhomer.csv",quote=F, sep = "\t",row.names = F,col.names = F)
x.p1 <- unlist(strsplit(Oticsensory.sp.gene.peak$peakName,"_"));x.p2 <- as.data.frame(cbind(x.p1[seq(1,length(x.p1),3)],x.p1[seq(2,length(x.p1),3)],x.p1[seq(3,length(x.p1),3)]))
write.table(x.p2,file = "p2g.cluster.homer/Oticsensory.forhomer.csv",quote=F, sep = "\t",row.names = F,col.names = F)
x.p1 <- unlist(strsplit(Hair.sp.gene.peak$peakName,"_"));x.p2 <- as.data.frame(cbind(x.p1[seq(1,length(x.p1),3)],x.p1[seq(2,length(x.p1),3)],x.p1[seq(3,length(x.p1),3)]))
write.table(x.p2,file = "p2g.cluster.homer/Hair.forhomer.csv",quote=F, sep = "\t",row.names = F,col.names = F)
x.p1 <- unlist(strsplit(Utricle.sp.gene.peak$peakName,"_"));x.p2 <- as.data.frame(cbind(x.p1[seq(1,length(x.p1),3)],x.p1[seq(2,length(x.p1),3)],x.p1[seq(3,length(x.p1),3)]))
write.table(x.p2,file = "p2g.cluster.homer/Utricle.forhomer.csv",quote=F, sep = "\t",row.names = F,col.names = F)
#
#
##
#####
#############################################
# Peak2GeneHeatmap
pp <- plotPeak2GeneHeatmap(ArchRProj = MOPA.epithelial, 
                           corCutOff = 0.4,
                           varCutOffATAC = 0.15,
                           varCutOffRNA = 0.15,
                           FDRCutOff = 1e-04,
                           groupBy = "Rank.ID",
                           k = 25,
                           palGroup = col.epi.27,
                           nPlot = 25000,
                           returnMatrices = F)
pp
ComplexHeatmap::draw(pp, heatmap_legend_side = "bot", annotation_legend_side = "bot")
#
#
p <- plotPeak2GeneHeatmap(ArchRProj = MOPA.epithelial, 
                          corCutOff = 0.4,
                          varCutOffATAC = 0.15,
                          varCutOffRNA = 0.15,
                          FDRCutOff = 1e-04,
                          groupBy = "Rank.ID",
                          k = 25,
                          palGroup = col.epi.27,
                          nPlot = 25000,
                          returnMatrices = T)
p # 73104
##
dim(p$ATAC$matrix)
dim(p$RNA$matrix)
#
p@listData$ATAC$kmeansId
#
#
## 20230412 
### plot peak to gene hetmap
# fot ATAC-seq
library(ComplexHeatmap)
pAm <- p$ATAC$matrix
head(pAm)
pAm.row.rank <- c()
for(i in 1:25){
  x <- which(p$ATAC$kmeansId==i)
  pAm.row.rank <- c(pAm.row.rank,x)
}
pAm <- pAm[pAm.row.rank,]
head(pAm)
dim(pAm)
colnames(pAm)
#
pAm.col.rank <- read.table("p2grank.csv")
pAm <- pAm[,pAm.col.rank$V1]
head(pAm)
#
#write.csv(p$RNA$colData$groupBy,"P2G.Rank.csv")
#pheatmap(pAm[seq(1,24989,25),seq(1,500,5)],cluster_rows = F,cluster_cols = F)
#ComplexHeatmap::Heatmap(pAm[seq(1,24989,25),seq(1,500,5)],cluster_rows = F,cluster_columns = F)
#ComplexHeatmap::Heatmap(pAm,cluster_rows = F,cluster_columns = F)
#
#
#
col_fun = colorRamp2(seq(-127,128,1)/64, paletteContinuous("solarExtra"))
split.module <- sort(p$ATAC$kmeansId)
split.module  <- factor(split.module,levels = 1:25)
#
p2g.cluster.rank <- c("T.C20","R.C18","V.C22","S.C19","I.C9",
                      "M.C13","N.C14","F.C6","G.C7","Q.C17","K.C11",
                      "O.C15","C.C3","P.C16","H.C8","U.C21","E.C5",
                      "J.C10","D.C4",
                      "A.C1","B.C2","L.C12")
p2g.cluster.rank.tmp1 <- as.data.frame(table(p$RNA$colData$groupBy))
p2g.cluster.rank.tmp1$color <- c ("#32CD32", "#008B8B", "#FFE4B5", "#40E0D0", "#FF6347", # 1-5
                                  "#FF1493", "#8B008B", "#FFA500", "#6A5ACD", "#DEB887", # 7-11
                                  "#0000CD", "#9400D3", "#87CEEB", "#20B2AA", "#E9967A", # 12-16
                                  "#800080", "#00FFFF", "#F08080", "#228B22", "#1E90FF", # 17-18-19-22-23
                                  "#8B4513", "#FF00FF")
rownames(p2g.cluster.rank.tmp1) <- p2g.cluster.rank.tmp1$Var1
p2g.cluster.rank.tmp1 <- p2g.cluster.rank.tmp1[p2g.cluster.rank,]
p2g.cluster.rank.tmp1
#
col.use.p2g.col <- p2g.cluster.rank.tmp1$color
names(col.use.p2g.col) <- p2g.cluster.rank
col.use.p2g.col <- list("CellType"=col.use.p2g.col)
#
p2g.cluster.rank2 <- c()
p2g.cluster.rank2  <- as.character(p2g.cluster.rank2)
x <- c()
for(i in 1:22){
  x <- rep(p2g.cluster.rank.tmp1$Var1[i],p2g.cluster.rank.tmp1$Freq[i])
  p2g.cluster.rank2 <- as.character(c(p2g.cluster.rank2,as.character(x)))
}
sample_info.z = data.frame(CellType = p2g.cluster.rank2 )
top_color <- HeatmapAnnotation(df=sample_info.z,col = col.use.p2g.col,annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"))
#
library(ComplexHeatmap)
dim(pAm)
p2g.ATAC.plot <- Heatmap(pAm,
                         col = col_fun,
                         cluster_rows = F,
                         cluster_columns = F,
                         show_column_names = F,
                         show_row_names = F,
                         row_split = split.module,
                         gap = unit(0.5, "mm"), 
                         border = "black",
                         use_raster = T,
                         top_annotation=top_color ## 控制cluster的annotation
)
draw(p2g.ATAC.plot, heatmap_legend_side = "right", show_annotation_legend = T)
#
# 
# for RNA-seq
pRm <- p$RNA$matrix
head(pRm)
pRm.row.rank <- c()
for(i in 1:25){
  x <- which(p$RNA$kmeansId==i)
  pRm.row.rank <- c(pRm.row.rank,x)
}
pRm <- pRm[pRm.row.rank,]
head(pRm)
#
pRm.col.rank <- read.table("p2grank.csv")
pRm <- pRm[,pRm.col.rank$V1]
head(pRm)
#
col_fun = colorRamp2(seq(-127,128,1)/64, paletteContinuous("blueYellow"))
split.module <- sort(p$ATAC$kmeansId)
split.module  <- factor(split.module,levels = 1:25)
#
p2g.cluster.rank <- c("T.C20","R.C18","V.C22","S.C19","I.C9",
                      "M.C13","N.C14","F.C6","G.C7","Q.C17","K.C11",
                      "O.C15","C.C3","P.C16","H.C8","U.C21","E.C5",
                      "J.C10","D.C4",
                      "A.C1","B.C2","L.C12")
p2g.cluster.rank.tmp1 <- as.data.frame(table(p$RNA$colData$groupBy))
p2g.cluster.rank.tmp1$color <- c ("#32CD32", "#008B8B", "#FFE4B5", "#40E0D0", "#FF6347", # 1-5
                                  "#FF1493", "#8B008B", "#FFA500", "#6A5ACD", "#DEB887", # 7-11
                                  "#0000CD", "#9400D3", "#87CEEB", "#20B2AA", "#E9967A", # 12-16
                                  "#800080", "#00FFFF", "#F08080", "#228B22", "#1E90FF", # 17-18-19-22-23
                                  "#8B4513", "#FF00FF")
rownames(p2g.cluster.rank.tmp1) <- p2g.cluster.rank.tmp1$Var1
p2g.cluster.rank.tmp1 <- p2g.cluster.rank.tmp1[p2g.cluster.rank,]
p2g.cluster.rank.tmp1
#
col.use.p2g.col <- p2g.cluster.rank.tmp1$color
names(col.use.p2g.col) <- p2g.cluster.rank
col.use.p2g.col <- list("CellType"=col.use.p2g.col)
#
p2g.cluster.rank2 <- c()
p2g.cluster.rank2  <- as.character(p2g.cluster.rank2)
x <- c()
for(i in 1:22){
  x <- rep(p2g.cluster.rank.tmp1$Var1[i],p2g.cluster.rank.tmp1$Freq[i])
  p2g.cluster.rank2 <- as.character(c(p2g.cluster.rank2,as.character(x)))
}
sample_info.z = data.frame(CellType = p2g.cluster.rank2 )
top_color <- HeatmapAnnotation(df=sample_info.z,col = col.use.p2g.col,annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"))
#
p2g.RNA.plot <- Heatmap(pRm,
                        col = col_fun,
                        cluster_rows = F,
                        cluster_columns = F,
                        show_column_names = F,
                        show_row_names = F,
                        row_split = split.module,
                        gap = unit(0.5, "mm"), 
                        border = "black",
                        top_annotation=top_color
)
draw(p2g.RNA.plot, heatmap_legend_side = "right", show_annotation_legend = T)
#
write.csv(pAm,file = "SourceData/Fig5_P2G_ATAC_Heatmap.csv")
write.csv(pRm,file = "SourceData/Fig5_P2G_RNA_Heatmap.csv")
#
#
#
MOPA.epithelial <- addMotifAnnotations(ArchRProj = MOPA.epithelial, motifSet = "cisbp", name = "Motif",force = TRUE)
MOPA.epithelial <- addBgdPeaks(MOPA.epithelial,force=T)
MOPA.epithelial <- addDeviationsMatrix(ArchRProj = MOPA.epithelial, peakAnnotation = "Motif",force = T)
#
#
plotEmbedding(ArchRProj = MOPA.epithelial, colorBy = "MotifMatrix", embedding = "PeakMatrix_LSI.M2_UMAP1",
              name = "deviations:Nr5a2_675", imputeWeights = getImputeWeights(MOPA.epithelial))
p
#
#
p <- plotBrowserTrack(
  ArchRProj = MOPA.epithelial, 
  groupBy = "ID", 
  geneSymbol = c("Nr5a2"), 
  upstream = 27000,
  downstream = 50000,
  ylim = c(0.1,1),
  loops = getPeak2GeneLinks(MOPA.epithelial)
)
grid::grid.newpage()
grid::grid.draw(p$Nr5a2)
#
#
p <- plotBrowserTrack(ArchRProj = MOPA.epithelial, groupBy = "Rank.ID", geneSymbol = c("Nr5a2"), 
                      upstream = 20000,downstream = 45000,ylim = c(0.1,1),pal =col.epi.27.plot,
                      loops = getPeak2GeneLinks(MOPA.epithelial,corCutOff = 0.4,FDRCutOff = 1e-04,
                                                varCutOffATAC = 0.15,varCutOffRNA = 0.15,))
grid::grid.newpage()
grid::grid.draw(p$Nr5a2)
#
p <- plotBrowserTrack(ArchRProj = MOPA.epithelial, groupBy = "Rank.ID", geneSymbol = c("Nkx2-1"), 
                      upstream = 60000,downstream = 15000,ylim = c(0.1,1),pal =col.epi.27.plot,
                      loops = getPeak2GeneLinks(MOPA.epithelial,corCutOff = 0.4,FDRCutOff = 1e-04,
                                                varCutOffATAC = 0.15,varCutOffRNA = 0.15,))
grid::grid.newpage()
grid::grid.draw(p$`Nkx2-1`)
#
p <- plotBrowserTrack(ArchRProj = MOPA.epithelial, groupBy = "Rank.ID", geneSymbol = c("Grhl1"), 
                      upstream = 40000,downstream = 27000,ylim = c(0.1,1),pal =col.epi.27.plot,
                      loops = getPeak2GeneLinks(MOPA.epithelial,corCutOff = 0.4,FDRCutOff = 1e-04,
                                                varCutOffATAC = 0.15,varCutOffRNA = 0.15,))
grid::grid.newpage()
grid::grid.draw(p$Grhl1)
#
#
#
#
#
#
#### Source Data in 20231015
# p2geneDF.tmp2
write.csv(p2geneDF.tmp2,file = "SourceData/Epithelial_cCRE_to_Gene_list.csv")
write.csv(RNA.peak.num,file = "SourceData/Fig5E.csv")
write.csv(Peak.RNA.num,file = "SourceData/Fig5F.csv")
write.csv(p2geneDF.tmp2$peaktoTTSSdis,file = "SourceData/Fig5G.csv")
#
boxplot(subset(cor.values.meta2,cor.values.meta2$Freq <= 5)$cor.values,
        subset(cor.values.meta2,cor.values.meta2$Freq > 5)$cor.values,
        subset(cor.values.meta2,cor.values.meta2$Freq > 10)$cor.values,
        subset(cor.values.meta2,cor.values.meta2$Freq > 15)$cor.values,
        subset(cor.values.meta2,cor.values.meta2$Freq > 20)$cor.values,
        subset(cor.values.meta2,cor.values.meta2$Freq > 25)$cor.values,
        subset(cor.values.meta2,cor.values.meta2$Freq > 30)$cor.values,
        outline =F,
        col=brewer.pal(9,"Reds")[2:9]) # YlOrRd
#
write.csv(subset(cor.values.meta2,cor.values.meta2$Freq <= 5)$cor.values,
          file = "SourceData/Fig5H_bar1.csv")
write.csv(subset(cor.values.meta2,cor.values.meta2$Freq > 5)$cor.values,
          file = "SourceData/Fig5H_bar2.csv")
write.csv(subset(cor.values.meta2,cor.values.meta2$Freq > 10)$cor.values,
          file = "SourceData/Fig5H_bar3.csv")
write.csv(subset(cor.values.meta2,cor.values.meta2$Freq > 15)$cor.values,
          file = "SourceData/Fig5H_bar4.csv")
write.csv(subset(cor.values.meta2,cor.values.meta2$Freq > 20)$cor.values,
          file = "SourceData/Fig5H_bar5.csv")
write.csv(subset(cor.values.meta2,cor.values.meta2$Freq > 25)$cor.values,
          file = "SourceData/Fig5H_bar6.csv")
write.csv(subset(cor.values.meta2,cor.values.meta2$Freq > 30)$cor.values,
          file = "SourceData/Fig5H_bar7.csv")
#
MOPA.epithelial 
rownames(Fiugre4E_S4G) <- colnames(immune.combined)
head(Fiugre4E_S4G)
Fiugre4E_S4G <- Fiugre4E_S4G[MOPA.epithelial$cellNames,]
Fiugre4E_S4G <- as.data.frame(Fiugre4E_S4G)
MOPA.epithelial@embeddings$PeakMatrix_LSI.M1_UMAP1$df$`PeakMatrix_LSI.M1#UMAP_Dimension_1` <- Fiugre4E_S4G[,3]
MOPA.epithelial@embeddings$PeakMatrix_LSI.M1_UMAP1$df$`PeakMatrix_LSI.M1#UMAP_Dimension_2` <- Fiugre4E_S4G[,4]
#
TFlist <- c("deviations:Nr5a2_675","deviations:Nkx21_399","deviations:Grhl1_390","deviations:Trp63_853",
            "deviations:Gata4_386","deviations:Gata6_382","deviations:Nfix_722","deviations:Hnf4g_664",
            "deviations:Pax8_854")
for ( i in 1:8) {
  p <- plotEmbedding(ArchRProj = MOPA.epithelial, colorBy = "MotifMatrix",size = 0.0001, name = TFlist[i], 
                     embedding = "PeakMatrix_LSI.M1_UMAP1",imputeWeights = getImputeWeights(MOPA.epithelial))+NoLegend()+
    theme(plot.title = element_blank())+
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
    theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
    theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
    theme(panel.background = element_rect( colour = "black", size = 0))+
    theme(panel.border = element_blank())
  ggsave(paste("Motif_",TFlist[i],".png",sep = ""),plot=p,device="png",
         dpi=300,units = "cm",width = 7,height = 7)
}
#
TFlist2 <- c("Nr5a2","Nkx2-1","Grhl1","Trp63",
             "Gata4","Gata6","Nfix","Hnf4g","Pax8")
for ( i in 1:8) {
  p <- plotEmbedding(ArchRProj = MOPA.epithelial, colorBy = "GIM_M1_LSI.M2",size = 0.0001, name = TFlist2[i], pal =c("lightgrey", "blue"),
                     embedding = "PeakMatrix_LSI.M1_UMAP1",imputeWeights = getImputeWeights(MOPA.epithelial))+NoLegend()+
    theme(plot.title = element_blank())+
    theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
    theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
    theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
    theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
    theme(panel.background = element_rect( colour = "black", size = 0))+
    theme(panel.border = element_blank())
  ggsave(paste("GIM_",TFlist2[i],".png",sep = ""),plot=p,device="png",
         dpi=300,units = "cm",width = 7,height = 7)
}
#
p <- plotEmbedding(ArchRProj = MOPA.epithelial, colorBy = "GIM_M1_LSI.M2",size = 0.0001, name = TFlist2[i], 
                   quantCut = c(0.4, 0.99),pal =c("lightgrey", "blue"),
                   embedding = "PeakMatrix_LSI.M1_UMAP1",imputeWeights = getImputeWeights(MOPA.epithelial))+NoLegend()+
  theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
ggsave(paste("GIM_",TFlist2[i],".png",sep = ""),plot=p,device="png",
       dpi=300,units = "cm",width = 7,height = 7)
#
#
epi.peak.rev <- peakset20230412.metadata
epi.peak.rev <- epi.peak.rev[,c(1,2,3)]
head(epi.peak.rev)
write.table(epi.peak.rev,file = "p2g.cluster.homer/epi.peak.rev.csv",sep = "\t",quote =F,row.names =F,col.names = F)
#
#
write.csv(cor.values.meta2,file = "SourceData/Figure.Epi_cor_list_tableS8.csv")
#
cor.values.meta2x <- cor.values.meta2[unique(top3$gene)[1:500],]
cor.values.meta2x[which(cor.values.meta2x$Freq>5),]$genename
write.csv(cor.values.meta2x[which(cor.values.meta2x$Freq>5),]$genename,file = "SourceData/Figure.Epi_GO.csv")
#
boxplot(cor.values.meta2[unique(top3$gene)[1:500],]$Freq,
        cor.values.meta2[setdiff(cor.values.meta2$genename,unique(top3$gene)[1:500]),]$Freq,
        outline =F)
#
wilcox.test(cor.values.meta2[unique(top3$gene)[1:500],]$Freq,
            cor.values.meta2[setdiff(cor.values.meta2$genename,unique(top3$gene)[1:500]),]$Freq)
#
write.csv(cor.values.meta2[unique(top3$gene)[1:500],]$Freq, file = "SourceData/Figure.New.5L_left.csv")
write.csv(cor.values.meta2[setdiff(cor.values.meta2$genename,unique(top3$gene)[1:500]),]$Freq,file = "SourceData/Figure.New.5L_right.csv")
#
#
save.image("MOPA.epithelial.RData")


