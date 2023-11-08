setwd("/home/sunkeyong/MOPA_project/")
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
addArchRGenome("mm10")
#
# Based on the results of Figure 1
load("MOPA.RData")
# Here first performing peak calling baesd on 296 cell subtypes.
MOPA <- addGroupCoverages(ArchRProj = MOPA, 
                          groupBy = "MOPA.296_subclusters",
                          minCells = 40,
                          maxCells = 500,
                          minReplicates = 2,
                          maxReplicates = 5,
                          sampleRatio = 0.8)
MOPA <- addReproduciblePeakSet(ArchRProj = MOPA, 
                               groupBy = "MOPA.296_subclusters",
                               pathToMacs2 = "/home/sunkeyong/.local/bin/macs2",
                               reproducibility = "2",
                               minCells = 40,
                               shift = -75,
                               extsize = 150,
                               additionalParams = "--nomodel --nolambda",
                               extendSummits = 250,
                               cutOff = 0.1)
#
MOPA <- addPeakMatrix(ArchRProj = MOPA, force =T)
#
MOPA_peakset <- getPeakSet(MOPA)
MOPA_peakset <- as.data.frame(MOPA_peakset)
table(MOPA_peakset$peakType) # for Figure 2A
#
#
MOPA$BW.ID <- paste(MOPA$ID,"_in_",MOPA$Timepoint,sep = "")
table(MOPA$BW.ID)
#
#
MOPA.BW.ID.matrix <- getGroupSE(ArchRProj = MOPA,
  useMatrix = "PeakMatrix",
  groupBy = "BW.ID",
  scaleTo = 1e6)
#
axsum.filter <- MOPA.BW.ID.matrix@assays@data$PeakMatrix
rownames(axsum.filter) <- paste(MOPA.BW.ID.matrix@elementMetadata$seqnames,
                                MOPA.BW.ID.matrix@elementMetadata$start,
                                MOPA.BW.ID.matrix@elementMetadata$end,sep="-")
axsum.filter.row.max <- apply(axsum.filter , 1, max)
km_result2 <- kmeans(df,150,nstart = 24)
#
dd <- cbind(axsum.filter,cluster=km_result2$cluster)
head(dd)
dd <- as.data.frame(dd)
min(table(dd$cluster))
#
# pick one cluster for shoing
df.1 <- df[rownames(subset(dd,dd$cluster== "20")),]
df.1.sub <- df.1[sample(dim(df.1)[1],150),]
pheatmap(df.1.sub,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=F,show_colnames=F)
#
pheatmap(df.1.sub,cluster_rows = T,cluster_cols = F,annotation_names_row = T,show_rownames=F,show_colnames=T)
#
#
# showing 400 peaks in each cluster
df.x.sub <- c()  
for(i in 1:length(unique(dd$cluster))){
  df.x <- df[rownames(subset(dd,dd$cluster== i)),]
  df.x.sub <- c(df.x.sub,rownames(df.x)[sample(dim(df.x)[1],400)])
}
df.x.sub.m <- df[df.x.sub,]
dim(df.x.sub.m)
pheatmap(df.x.sub.m,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=F,show_colnames=F)
#
#
###
# calculating the mean aceessibility of each cluster of K-means
df.x.m <- c()
df.m <- c()
for(i in 1:length(unique(dd$cluster))){
  df.x.m <- colMeans(df[rownames(subset(dd,dd$cluster== i)),])
  df.m <- rbind(df.m,df.x.m)
}
dim(df.m)
head(df.m)
pheatmap(df.m,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=F,show_colnames=T)
pheatmap(df.m,cluster_rows = T,cluster_cols = F,annotation_names_row = F,show_rownames=F,show_colnames=T)
pheatmap(df.m,cluster_rows = T,cluster_cols = T,annotation_names_row = F,show_rownames=F,show_colnames=T)
#
pheatmap(t(df.m),cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=T,show_colnames=T)
#
df.m.1 <- df.m
rownames(df.m.1) <- paste("M",seq(1:150),sep = "")
#
pheatmap(df.m.1,cluster_rows = T,cluster_cols = T,annotation_names_row = F,show_rownames=T,show_colnames=T)
#
head(df.m.1)
Heatmap(df.m.1,row_dend_reorder = T)
pheatmap(cor(t(df.m.1)))
#
# merge similar module
head(df.m.1)
head(dd)
module.merge1 <- read.csv("peakmodule.csv")
head(module.merge1)
module.merge1 <- as.data.frame(module.merge1)
head(module.merge1)
dd.mm1 <- dd
dd.mm1$peak <- rownames(dd.mm1)
dd.mm1 <- merge(dd.mm1,module.merge1,by="cluster")
head(dd.mm1)
rownames(dd.mm1) <- dd.mm1$peak
head(dd.mm1)
##
###
#
# re-calculating the mean aceessibility of each cluster of K-means
dd.mm1.x.m <- c()
dd.mm1.m <- c()
length(unique(dd.mm1$cluster.mt9))
for(i in 1:length(unique(dd.mm1$cluster.mt9))){
  dd.mm1.x.m  <- colMeans(df[rownames(subset(dd.mm1,dd.mm1$cluster.mt9== i)),])
  dd.mm1.m <- rbind(dd.mm1.m,dd.mm1.x.m )
}
dim(dd.mm1.m)
head(dd.mm1.m)
#
rownames(dd.mm1.m) <- paste("M",seq(1:89),sep = "")
#
pheatmap(dd.mm1.m,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=F,show_colnames=T)
#
pheatmap(t(dd.mm1.m),cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=F,show_colnames=T)
pheatmap(t(dd.mm1.m),cluster_rows = F,cluster_cols = T,annotation_names_row = F,show_rownames=T,show_colnames=T)
#
#
pheatmap(dd.mm1.m,cluster_rows = T,cluster_cols = F,annotation_names_row = F,show_rownames=F,show_colnames=T)
pheatmap(dd.mm1.m,cluster_rows = T,cluster_cols = T,annotation_names_row = F,show_rownames=F,show_colnames=T)
#
#
cluster.rank <- c("Neural Tube_E10.5","Neural Tube_E11.5","Neural Tube_E12.5","Neural Tube_E13.5",
                  "Notochord cells_E10.5","Notochord cells_E11.5","Notochord cells_E12.5","Notochord cells_E13.5",
                  "Neural progenitor cells_E10.5","Neural progenitor cells_E11.5","Neural progenitor cells_E12.5","Neural progenitor cells_E13.5",
                  "Radial glia_E10.5","Radial glia_E11.5","Radial glia_E12.5","Radial glia_E13.5",
                  "Postmitotic premature neurons_E10.5","Postmitotic premature neurons_E11.5","Postmitotic premature neurons_E12.5","Postmitotic premature neurons_E13.5",                  "Sensory neurons_E10.5","Sensory neurons_E11.5","Sensory neurons_E12.5","Sensory neurons_E13.5",
                  "Oligodendrocyte Progenitors_E10.5","Oligodendrocyte Progenitors_E11.5","Oligodendrocyte Progenitors_E12.5","Oligodendrocyte Progenitors_E13.5",
                  "Premature oligodendrocyte_E10.5","Premature oligodendrocyte_E11.5","Premature oligodendrocyte_E12.5","Premature oligodendrocyte_E13.5",
                  "Schwann cell precursor_E10.5","Schwann cell precursor_E11.5","Schwann cell precursor_E12.5","Schwann cell precursor_E13.5",
                  "Inhibitory neuron progenitors_E10.5","Inhibitory neuron progenitors_E11.5","Inhibitory neuron progenitors_E12.5","Inhibitory neuron progenitors_E13.5",
                  "Inhibitory interneurons_E10.5","Inhibitory interneurons_E11.5","Inhibitory interneurons_E12.5","Inhibitory interneurons_E13.5",
                  "Inhibitory neurons_E10.5","Inhibitory neurons_E11.5","Inhibitory neurons_E12.5","Inhibitory neurons_E13.5",
                  "Excitatory neurons_E10.5","Excitatory neurons_E11.5","Excitatory neurons_E12.5","Excitatory neurons_E13.5",
                  "Cholinergic neurons_E10.5","Cholinergic neurons_E11.5","Cholinergic neurons_E12.5","Cholinergic neurons_E13.5",
                  "Granule neurons_E10.5","Granule neurons_E11.5","Granule neurons_E12.5","Granule neurons_E13.5",
                  "Ependymal cell_E10.5","Ependymal cell_E11.5","Ependymal cell_E12.5","Ependymal cell_E13.5",
                  "Isthmic organizer cells_E10.5","Isthmic organizer cells_E11.5","Isthmic organizer cells_E12.5","Isthmic organizer cells_E13.5",
                  "Connective tissue progenitors_E10.5","Connective tissue progenitors_E11.5","Connective tissue progenitors_E12.5","Connective tissue progenitors_E13.5",
                  "Chondrocytes & osteoblasts_E10.5","Chondrocytes & osteoblasts_E11.5","Chondrocytes & osteoblasts_E12.5","Chondrocytes & osteoblasts_E13.5",
                  "Chondroctye progenitors_E10.5","Chondroctye progenitors_E11.5","Chondroctye progenitors_E12.5","Chondroctye progenitors_E13.5",
                  "Intermediate Mesoderm_E10.5","Intermediate Mesoderm_E11.5","Intermediate Mesoderm_E12.5","Intermediate Mesoderm_E13.5",
                  "Limb mesenchyme_E10.5","Limb mesenchyme_E11.5","Limb mesenchyme_E12.5","Limb mesenchyme_E13.5",
                  "Jaw and tooth progenitors_E10.5","Jaw and tooth progenitors_E11.5","Jaw and tooth progenitors_E12.5","Jaw and tooth progenitors_E13.5",
                  "Early mesenchyme_E10.5","Early mesenchyme_E11.5","Early mesenchyme_E12.5","Early mesenchyme_E13.5",
                  "Myocytes_E10.5","Myocytes_E11.5","Myocytes_E12.5","Myocytes_E13.5",
                  "Cardiac muscle lineages_E10.5","Cardiac muscle lineages_E11.5","Cardiac muscle lineages_E12.5","Cardiac muscle lineages_E13.5",
                  "Epithelial cells_E10.5","Epithelial cells_E11.5","Epithelial cells_E12.5","Epithelial cells_E13.5",
                  "Melanocytes_E10.5","Melanocytes_E11.5","Melanocytes_E12.5","Melanocytes_E13.5",
                  "Hepatocytes_E10.5","Hepatocytes_E11.5","Hepatocytes_E12.5","Hepatocytes_E13.5",
                  "Endothelial cells_E10.5","Endothelial cells_E11.5","Endothelial cells_E12.5","Endothelial cells_E13.5",
                  "White blood cells_E10.5","White blood cells_E11.5","White blood cells_E12.5","White blood cells_E13.5",
                  "Definitive erythroid lineage_E10.5","Definitive erythroid lineage_E11.5","Definitive erythroid lineage_E12.5","Definitive erythroid lineage_E13.5",
                  "Primitive erythroid lineage_E10.5","Primitive erythroid lineage_E11.5","Primitive erythroid lineage_E12.5","Primitive erythroid lineage_E13.5")
length(unique(cluster.rank))
length(cluster.rank)
#
module.rank <- c("M8","M14","M74","M73","M68","M29","M47","M46","M62","M42","M78",
                 "M89","M83","M37","M41","M36","M26","M67","M38","M53","M4","M28","M71","M10",
                 "M55","M65","M15","M60","M39","M18","M33","M66","M69","M3","M5","M85","M20","M86",
                 "M64","M34","M87","M9","M63","M6","M48","M31","M58","M30","M40",
                 "M50","M19","M44","M52","M72","M16","M51","M17","M23", "M7","M56",
                 "M76","M45","M75","M11","M57","M1","M13","M84","M35","M61",
                 "M88","M2","M79","M54","M24","M70","M82","M81","M80","M49",
                 "M27","M12","M43","M21","M77","M32","M59","M22","M25")
length(unique(module.rank ))
length(module.rank )
dd.mm1.m.1 <- dd.mm1.m
#
dd.mm1.m.2 <- dd.mm1.m.1
dd.mm1.m.2 <- dd.mm1.m.2[,cluster.rank]
#
dd.mm1.m.2 <- dd.mm1.m.2[module.rank,]
#
pheatmap(dd.mm1.m.2,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=T,show_colnames=T)
#
peak.plot1 <- Heatmap(dd.mm1.m.2,
                      cluster_rows = F,
                      cluster_columns = F,
                      show_column_names = F,
                      show_row_names = F,
                      col = paletteContinuous("solarExtra")[85:256])
draw(peak.plot1, heatmap_legend_side = "right", show_annotation_legend = T)
#
color.use.em <- c("#FDCDAC","#FFFF99","#B3CDE3","#6495ED","#00FFFF", ## 1--5 
                  "#FBB4AE","#F4CAE4","#ABDDA4","#B15928","#D6604D", ## 6--10
                  "#999999","#00BFFF","#CC99FF","#48D1CC","#FFED6F", ## 11--15 
                  "#9E0142","#0099FF","#B3E2CD","#276419","#66C2A5", ## 16--20
                  "#B3EE3A","#6A3D9A","#87CEFA","#ABD9E9","#B2ABD2", ## 21--25 
                  "#0000FF","#FF33CC","#66BD63","#E08214","#DE77AE", ## 26--30
                  "#9970AB","#EE82EE","#99CCCC") 
color.use.1  <- as.character(color.use.em)
color.use.1  <-  rep(color.use.1,each =4) 
#
sample_info.z = data.frame(CellType = cluster.rank )
col.use.1.tmp <- color.use.1
names(col.use.1.tmp) <- cluster.rank
col.use.1.tmp
class(col.use.1.tmp)
#col.use.1.tmp2 <- list("CellType"=col.use.1.tmp)
#col.use.1.tmp2
## 
timepoint <- c("E10.5","E11.5","E12.5","E13.5")
timepoint.color <-  rep(as.character(color.for.zebrafish.F1[1:4]),33)
names(timepoint.color) <- rep(c("E10.5","E11.5","E12.5","E13.5"),33)
#
anot_df <- data.frame(celltype=cluster.rank,timepoint=  rep(timepoint,33) )
anot_col <- list(celltype= rep(as.character(color.for.zebrafish.F1[1:33]),each =4),
                 timepoint= rep(as.character(color.for.zebrafish.F1[1:4]),33))
ha <- HeatmapAnnotation(anot_df,col=anot_col)
#
timepoint <- c("E10.5","E11.5","E12.5","E13.5")
ha <- HeatmapAnnotation(celltype=cluster.rank,timepoint= timepoint,
                        col=list(celltype= col.use.1.tmp,
                                 timepoint=timepoint.color))
#
#
library(RColorBrewer)
peak.plot1 <- Heatmap(dd.mm1.m.2,
                      cluster_rows = F,
                      cluster_columns = F,
                      #row_split = split.module,
                      show_column_names = T,
                      show_row_names = F,
                      gap = unit(0.5, "mm"), 
                      border = "black",
                      top_annotation=ha, 
                      col = paletteContinuous("solarExtra")[85:256])
draw(peak.plot1, heatmap_legend_side = "right", show_annotation_legend = F)
#
module.rank.s <- c("8","14","20","74","73","68","29","62","46","47","42","78",
                   "89","83","37","41","36","26","67","38","53","4","28","71","10",
                   "65","15","60","55","39","18","33","66","69","3","5","85","86",
                   "64","34","87","9","63","6","48","31","58","30","40",
                   "50","19","44","52","72","16","51","17","23", "7","56",
                   "76","45","75","11","57","1","13","84","35","61",
                   "88","2","79","54","24","70","82","12","43","21","77","81","80","49",
                   "27","32","59","22","25")
##
cluster.mt10 <-  as.character(1:89)
module.rank.s.re <- as.data.frame(cbind(module.rank.s,cluster.mt10))
#
colnames(module.rank.s.re) <- c("cluster.mt9","cluster.mt10")
module.merge2 <- merge(module.merge1, module.rank.s.re,by="cluster.mt9")
head(module.merge2)
##
dd.mm1 <- dd
dd.mm1$peak <- rownames(dd.mm1)
dd.mm1 <- merge(dd.mm1,module.merge2,by="cluster")
head(dd.mm1)
rownames(dd.mm1) <- dd.mm1$peak
head(dd.mm1)
dim(dd.mm1)
##
###
#
# re-calculating the mean aceessibility of each cluster of K-means after ranking the modules
dd.mm1.x.m <- c()
dd.mm1.m <- c()
length(unique(dd.mm1$cluster.mt10))
for(i in 1:length(unique(dd.mm1$cluster.mt10))){
  dd.mm1.x.m  <- colMeans(df[rownames(subset(dd.mm1,dd.mm1$cluster.mt10== i)),])
  dd.mm1.m <- rbind(dd.mm1.m,dd.mm1.x.m )
}
dim(dd.mm1.m)
head(dd.mm1.m)
#
rownames(dd.mm1.m) <- paste("M",seq(1:89),sep = "")
#
pheatmap(dd.mm1.m,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=F,show_colnames=T)
pheatmap(dd.mm1.m,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=T,show_colnames=T)
#
cluster.rank <- c("Neural Tube_E10.5","Neural Tube_E11.5","Neural Tube_E12.5","Neural Tube_E13.5",
                  "Notochord cells_E10.5","Notochord cells_E11.5","Notochord cells_E12.5","Notochord cells_E13.5",
                  "Neural progenitor cells_E10.5","Neural progenitor cells_E11.5","Neural progenitor cells_E12.5","Neural progenitor cells_E13.5",
                  "Radial glia_E10.5","Radial glia_E11.5","Radial glia_E12.5","Radial glia_E13.5",
                  "Postmitotic premature neurons_E10.5","Postmitotic premature neurons_E11.5","Postmitotic premature neurons_E12.5","Postmitotic premature neurons_E13.5",                  "Sensory neurons_E10.5","Sensory neurons_E11.5","Sensory neurons_E12.5","Sensory neurons_E13.5",
                  "Oligodendrocyte Progenitors_E10.5","Oligodendrocyte Progenitors_E11.5","Oligodendrocyte Progenitors_E12.5","Oligodendrocyte Progenitors_E13.5",
                  "Premature oligodendrocyte_E10.5","Premature oligodendrocyte_E11.5","Premature oligodendrocyte_E12.5","Premature oligodendrocyte_E13.5",
                  "Schwann cell precursor_E10.5","Schwann cell precursor_E11.5","Schwann cell precursor_E12.5","Schwann cell precursor_E13.5",
                  "Inhibitory neuron progenitors_E10.5","Inhibitory neuron progenitors_E11.5","Inhibitory neuron progenitors_E12.5","Inhibitory neuron progenitors_E13.5",
                  "Inhibitory interneurons_E10.5","Inhibitory interneurons_E11.5","Inhibitory interneurons_E12.5","Inhibitory interneurons_E13.5",
                  "Inhibitory neurons_E10.5","Inhibitory neurons_E11.5","Inhibitory neurons_E12.5","Inhibitory neurons_E13.5",
                  "Excitatory neurons_E10.5","Excitatory neurons_E11.5","Excitatory neurons_E12.5","Excitatory neurons_E13.5",
                  "Cholinergic neurons_E10.5","Cholinergic neurons_E11.5","Cholinergic neurons_E12.5","Cholinergic neurons_E13.5",
                  "Granule neurons_E10.5","Granule neurons_E11.5","Granule neurons_E12.5","Granule neurons_E13.5",
                  "Ependymal cell_E10.5","Ependymal cell_E11.5","Ependymal cell_E12.5","Ependymal cell_E13.5",
                  "Isthmic organizer cells_E10.5","Isthmic organizer cells_E11.5","Isthmic organizer cells_E12.5","Isthmic organizer cells_E13.5",
                  "Connective tissue progenitors_E10.5","Connective tissue progenitors_E11.5","Connective tissue progenitors_E12.5","Connective tissue progenitors_E13.5",
                  "Chondrocytes & osteoblasts_E10.5","Chondrocytes & osteoblasts_E11.5","Chondrocytes & osteoblasts_E12.5","Chondrocytes & osteoblasts_E13.5",
                  "Chondroctye progenitors_E10.5","Chondroctye progenitors_E11.5","Chondroctye progenitors_E12.5","Chondroctye progenitors_E13.5",
                  "Intermediate Mesoderm_E10.5","Intermediate Mesoderm_E11.5","Intermediate Mesoderm_E12.5","Intermediate Mesoderm_E13.5",
                  "Limb mesenchyme_E10.5","Limb mesenchyme_E11.5","Limb mesenchyme_E12.5","Limb mesenchyme_E13.5",
                  "Jaw and tooth progenitors_E10.5","Jaw and tooth progenitors_E11.5","Jaw and tooth progenitors_E12.5","Jaw and tooth progenitors_E13.5",
                  "Early mesenchyme_E10.5","Early mesenchyme_E11.5","Early mesenchyme_E12.5","Early mesenchyme_E13.5",
                  "Myocytes_E10.5","Myocytes_E11.5","Myocytes_E12.5","Myocytes_E13.5",
                  "Cardiac muscle lineages_E10.5","Cardiac muscle lineages_E11.5","Cardiac muscle lineages_E12.5","Cardiac muscle lineages_E13.5",
                  "Epithelial cells_E10.5","Epithelial cells_E11.5","Epithelial cells_E12.5","Epithelial cells_E13.5",
                  "Melanocytes_E10.5","Melanocytes_E11.5","Melanocytes_E12.5","Melanocytes_E13.5",
                  "Hepatocytes_E10.5","Hepatocytes_E11.5","Hepatocytes_E12.5","Hepatocytes_E13.5",
                  "Endothelial cells_E10.5","Endothelial cells_E11.5","Endothelial cells_E12.5","Endothelial cells_E13.5",
                  "White blood cells_E10.5","White blood cells_E11.5","White blood cells_E12.5","White blood cells_E13.5",
                  "Definitive erythroid lineage_E10.5","Definitive erythroid lineage_E11.5","Definitive erythroid lineage_E12.5","Definitive erythroid lineage_E13.5",
                  "Primitive erythroid lineage_E10.5","Primitive erythroid lineage_E11.5","Primitive erythroid lineage_E12.5","Primitive erythroid lineage_E13.5")
length(unique(cluster.rank))
length(cluster.rank)
#
module.rank.new.F <- paste("M",as.character(1:89),sep = "")
dd.mm1.m.1 <- dd.mm1.m
dd.mm1.m.2 <- dd.mm1.m.1
dd.mm1.m.2 <- dd.mm1.m.2[,cluster.rank]
dd.mm1.m.2 <- dd.mm1.m.2[module.rank.new.F,]
pheatmap(dd.mm1.m.2,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=T,show_colnames=T)
pheatmap(dd.mm1.m.2,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=T,show_colnames=F)
#
color.use.em <- c("#FDCDAC","#FFFF99","#B3CDE3","#6495ED","#00FFFF", ## 1--5 
                  "#FBB4AE","#F4CAE4","#ABDDA4","#B15928","#D6604D", ## 6--10
                  "#999999","#00BFFF","#CC99FF","#48D1CC","#FFED6F", ## 11--15 
                  "#9E0142","#0099FF","#B3E2CD","#276419","#66C2A5", ## 16--20
                  "#B3EE3A","#6A3D9A","#87CEFA","#ABD9E9","#B2ABD2", ## 21--25 
                  "#0000FF","#FF33CC","#66BD63","#E08214","#DE77AE", ## 26--30
                  "#9970AB","#EE82EE","#99CCCC") 
color.use.1  <- color.use.em
color.use.1  <-  rep(color.use.1,each =4) 
#
#
cluster.rank.new <- c("Neural Tube","Notochord cells","Neural progenitor cells",
                      "Radial glia","Postmitotic premature neurons","Sensory neurons",
                      "Oligodendrocyte Progenitors","Premature oligodendrocyte",
                      "Schwann cell precursor","Inhibitory neuron progenitors",
                      "Inhibitory interneurons","Inhibitory neurons","Excitatory neurons",
                      "Cholinergic neurons","Granule neurons","Ependymal cell","Isthmic organizer cells",
                      "Connective tissue progenitors","Chondrocytes & osteoblasts",
                      "Chondroctye progenitors","Intermediate Mesoderm","Limb mesenchyme","Jaw and tooth progenitors",
                      "Early mesenchyme","Myocytes","Cardiac muscle lineages","Epithelial cells","Melanocytes",
                      "Hepatocytes","Endothelial cells","White blood cells",
                      "Definitive erythroid lineage","Primitive erythroid lineage")
cluster.rank.new.color <-  color.use.1
names(cluster.rank.new.color) <- rep(cluster.rank.new,each= 4)
#
timepoint <- c("E10.5","E11.5","E12.5","E13.5")
Timepoint.color <- c("#98304E","#DBCCE0","#48D1CC","#6861A9")
timepoint.color <-  rep(Timepoint.color,33)
names(timepoint.color) <- rep(c("E10.5","E11.5","E12.5","E13.5"),33)
#
#
ha <- HeatmapAnnotation(celltype=rep(cluster.rank.new,each= 4),timepoint= rep(timepoint,33),
                        col=list(celltype= cluster.rank.new.color,
                                 timepoint=timepoint.color))

#
peak.plot1 <- Heatmap(dd.mm1.m.2,
                      cluster_rows = F,
                      cluster_columns = F,
                      #row_split = split.module,
                      show_column_names = F,
                      #row_split= split.module,
                      show_row_names = T,
                      gap = unit(0.5, "mm"), 
                      border = "black",
                      top_annotation=ha, 
                      use_raster=T,
                      col = paletteContinuous("solarExtra")[85:256])
draw(peak.plot1, heatmap_legend_side = "right", show_annotation_legend = F)
#
write.csv(dd.mm1.m.2,file = "./Source_Data/Fig2F.csv")
#
module.rank.new.F <- paste("M",as.character(1:89),sep = "")
split.module  <- factor(module.rank.new.F)
split.module  <- factor(split.module,levels =module.rank.new.F )
#
######## mergr module from one common cluster fro rGREAT and homer analysis
peak.num.m.0613 <- read.csv("peakmergemodule0613.csv",sep=",")
head(peak.num.m.0613)
head(dd.mm1)
dim(dd.mm1)
dd.mm2 <- merge(dd.mm1,peak.num.m.0613,by="cluster.mt10")
head(dd.mm2)
rownames(dd.mm2) <- dd.mm2$peak
head(dd.mm2)
dim(dd.mm1)
dim(dd.mm2)
###
for(i in 1:length(unique(dd.mm2$cluster.mt11))){
  x.p1 <- unlist(strsplit(subset(dd.mm2,dd.mm2$cluster.mt11== i)$peak,"-"))
  x.p2 <- as.data.frame(cbind(x.p1[seq(1,length(x.p1),3)],x.p1[seq(2,length(x.p1),3)],x.p1[seq(3,length(x.p1),3)]))
  write.table(x.p2,file = paste("peakmodule0613/","M",i,".csv",sep = ""),quote=F, sep = "\t",row.names = F,col.names = F)
}
#
#
#
pheatmap(cbind(as.data.frame(peak.num.m.0613$cluster.mt12),as.data.frame(peak.num.m.0613$cluster.mt12)),
         cluster_cols = F,cluster_rows = F,scale = "column",show_rownames = T,show_colnames = F)
#
peaksave <- dd.mm2[,c("cluster.mt10","cluster.mt11")]
head(peaksave)
dim(peaksave)
#
MOPA_peakset
#
peak_chr_name <- rep(as.character(MOPA_peakset@seqnames@values),MOPA_peakset@seqnames@lengths)
peak_start_site <- MOPA_peakset@ranges@start
peak_end_site <- MOPA_peakset@ranges@start+500
peak.names <- paste(peak_chr_name,peak_start_site,peak_end_site,sep = "-")
###
MOPA_peakset.metadata <- as.data.frame(peak.names)
head(MOPA_peakset.metadata)
colnames(MOPA_peakset.metadata) <- "peakname"
MOPA_peakset.metadata$chr <- peak_chr_name
MOPA_peakset.metadata$start <- peak_start_site
MOPA_peakset.metadata$end <- peak_end_site
MOPA_peakset.metadata$nearestTSS <- MOPA_peakset$nearestTSS
MOPA_peakset.metadata$distToGeneStart <- MOPA_peakset$distToGeneStart
MOPA_peakset.metadata$peakType <- MOPA_peakset$peakType
MOPA_peakset.metadata$distToTSS <- MOPA_peakset$distToTSS
MOPA_peakset.metadata$GC <- MOPA_peakset$GC
#
head(MOPA_peakset.metadata)
rownames(MOPA_peakset.metadata) <- MOPA_peakset.metadata$peakname
MOPA_peakset.metadata.tmp <- MOPA_peakset.metadata[rownames(peaksave),]
dim(MOPA_peakset.metadata.tmp)
#
#
#
colnames(peaksave) <- c("Module_89","mergedModule_for_GREATandMotif")
head(peaksave)
peaksave$Peak_name <- rownames(peaksave)
peaksave$peaktype <- MOPA_peakset.metadata.tmp$peakType
peaksave$nearestTSS <- MOPA_peakset.metadata.tmp$nearestTSS
peaksave$distToGeneStart <-MOPA_peakset.metadata.tmp$distToGeneStart
peaksave$distToTSS <- MOPA_peakset.metadata.tmp$distToTSS
head(peaksave)
peaksave2 <- peaksave[,c("Peak_name","peaktype","nearestTSS","distToGeneStart","distToTSS","Module_89","mergedModule_for_GREATandMotif")]
head(peaksave2)
write.csv(peaksave2,file="MOPA.allpeak.csv",row.names = F)
#
save.image("MOPA.Figure2.RData")
#