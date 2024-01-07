# Figure 1
# the species mixing experiment
setwd("/home/sunkeyong/SPATAC/Figure1/K562/")
library(Seurat)
library(ArchR)
addArchRGenome("hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
addArchRThreads(threads = 25) 
###
Bed_files= c("./human/SPATACseqSL1.sort.bam.fragments.filter.sort.bed.gz",
             "./human/SPATACseqSL2.sort.bam.fragments.filter.sort.bed.gz",
             "./human/SPATACseqSL3.sort.bam.fragments.filter.sort.bed.gz",
             "./human/SPATACseqSL4.sort.bam.fragments.filter.sort.bed.gz",
             "./human/SPATACseqSL5.sort.bam.fragments.filter.sort.bed.gz",
             "./human/SPATACseqSL6.sort.bam.fragments.filter.sort.bed.gz",
             "./human/SPATACseqSL7.sort.bam.fragments.filter.sort.bed.gz",
             "./human/SPATACseqSL8.sort.bam.fragments.filter.sort.bed.gz",
             "./human/Jason.sort.bed.gz",
             "./human/T_C1.sort.bed.gz",
             "./human/T_plate.sort.bed.gz")
names(Bed_files) = c("SPATACseqSL1","SPATACseqSL2","SPATACseqSL3","SPATACseqSL4",
                     "SPATACseqSL5","SPATACseqSL6","SPATACseqSL7","SPATACseqSL8",
                     "Jason","T_C1","T_plate")
#
ArrowFiles <- createArrowFiles(inputFiles = Bed_files,
                                sampleNames = names(Bed_files),
                                minTSS = 4,
                                minFrags = 500,
                                addTileMat = TRUE,
                                addGeneScoreMat = TRUE,
                                force = T)
#
K562_ArchR <- ArchRProject(ArrowFiles = ArrowFiles,
                          copyArrows = TRUE)
K562_ArchR
#
table(K562_ArchR$Sample)
#
K562_ArchR$ReadsInTSS.ratio <- K562_ArchR$ReadsInTSS/K562_ArchR$nFrags
K562_ArchR$ReadsInPromoter.ratio <- K562_ArchR$ReadsInPromoter/K562_ArchR$nFrags
K562_ArchR$nFrags.log2 <- log2(K562_ArchR$nFrags)
####### Quality plot 
plotGroups(ArchRProj = K562_ArchR, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment",plotAs = "violin",alpha = 0.4,addBoxPlot = TRUE)
plotGroups(ArchRProj = K562_ArchR, groupBy = "Sample", colorBy = "cellColData", name = "nFrags.log2",plotAs = "violin",alpha = 0.4,addBoxPlot = TRUE)
plotGroups(ArchRProj = K562_ArchR, groupBy = "Sample", colorBy = "cellColData", name = "ReadsInTSS.ratio",plotAs = "violin",alpha = 0.4,addBoxPlot = TRUE)
plotGroups(ArchRProj = K562_ArchR, groupBy = "Sample", colorBy = "cellColData", name = "ReadsInPromoter.ratio",plotAs = "violin",alpha = 0.4,addBoxPlot = TRUE)
#
#### plot DNA fragment distribution and TSS enrichment
#
p1 <- plotFragmentSizes(ArchRProj = K562_ArchR,returnDF =F,maxSize = 1000)
p1
p2 <- plotTSSEnrichment(ArchRProj = K562_ArchR)
p2
idxPass <- which(K562_ArchR$Sample=="SPATACseqSL8" | K562_ArchR$Sample=="Jason" | K562_ArchR$Sample=="T_C1" | K562_ArchR$Sample=="T_plate")
cellsPass <- K562_ArchR$cellNames[idxPass]
K562_ArchR.sub <- K562_ArchR[cellsPass, ]
p3 <- plotFragmentSizes(ArchRProj = K562_ArchR.sub)
p3
p4 <- plotTSSEnrichment(ArchRProj = K562_ArchR.sub)
p4
#######
#
#
setwd("/home/sunkeyong/SPATAC/Figure1/K562/bedtools")
# for correlation analysis
# first calculated the coverage in the peakset using the following bash script 
#nohup bedtools coverage -a K562_peak.bed -b Encode-rep1.bed > Encode-rep1.coverage.bed &
# then load all coverage file
DNase_rep1 <- read.table("DNase-rep1.coverage.bed")
DNase_rep2 <- read.table("DNase-rep2.coverage.bed")
DNase_rep3 <- read.table("DNase-rep3.coverage.bed")
Encode_rep1 <- read.table("Encode-rep1.coverage.bed")
Encode_rep2 <- read.table("Encode-rep2.coverage.bed")
K562_rep1 <- read.table("K562-rep1.coverage.bed")
K562_rep2 <- read.table("K562-rep2.coverage.bed")
SL1 <- read.table("SL1.coverage.bed")
SL2 <- read.table("SL2.coverage.bed")
SL3 <- read.table("SL3.coverage.bed")
SL4 <- read.table("SL4.coverage.bed")
SL5 <- read.table("SL5.coverage.bed")
SL6 <- read.table("SL6.coverage.bed")
SL7 <- read.table("SL7.coverage.bed")
SL8 <- read.table("SL8.coverage.bed")
Jason <- read.table("Jason.coverage.bed")
TC1 <- read.table("T_C1.coverage.bed")
Tplate <- read.table("T_plate.coverage.bed")
dim(T_plate)
dim(T_C1)
#
# other.seq for negative control
other.seq <- cbind(sample(c(rep(0,dim(T_plate)[1]),rep(1,dim(T_plate)[1])),dim(T_plate)[1]),
                   sample(c(rep(0,dim(T_plate)[1]),rep(1,dim(T_plate)[1])),dim(T_plate)[1]))
colnames(other.seq) <- c("other1","other2")
head(other.seq)
dim(other.seq)
#
Coverage.merge <- cbind(DNase_rep1$V4,DNase_rep2$V4,DNase_rep3$V4,Encode_rep1$V4,Encode_rep2$V4,K562_rep1$V4,K562_rep2$V4,
                 SL1$V4,SL3$V4,SL3$V4,SL4$V4,SL5$V4,SL6$V4,SL7$V4,SL8$V4,
                 Jason$V4,T_C1$V4,T_plate$V4,other.seq)
#
head(Coverage.merge)
dim(Coverage.merge)
dim(Coverage.merge)
colnames(Coverage.merge) <- c("DNase_rep1","DNase_rep2","DNase_rep3",
                              "Encode_rep1","Encode_rep2",
                              "K562_rep1","K562_rep2",
                              "SL1","SL3","SL3","SL4","SL5","SL6","SL7","SL8",
                              "Jason","T_C1","T_plate","other1","other2")
head(Coverage.merge)
#
# For Figure 1E
library(pheatmap)
pheatmap(cor(Coverage.merge))
pheatmap(cor(log2(Coverage.merge+1)))
pheatmap(cor(log2(Coverage.merge+1)),cluster_rows = T,cluster_cols = T,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(300))
pheatmap(cor(log2(Coverage.merge+1)),cluster_rows = F,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(300))
#
SPATAC.SL8.cov.cpm <- SL8$V4/sum(SL8$V4)*1000000
K562_rep1.cov.cpm <- K562_rep1$V4/sum(K562_rep1$V4)*1000000
cor(log2(SPATAC.SL8.cov.cpm+1), log2(K562_rep1.cov.cpm+1))
#
# Figure 1D
p <- ggPoint(
  x = log2(SPATAC.SL8.cov.cpm+1), 
  y = log2(K562_rep1.cov.cpm+1), 
  colorDensity = TRUE,
  continuousSet = "sambaNight",
  xlabel = "SPATAC-seq SL8 log2(Unique Fragments)",
  ylabel = "K562_rep1 Bulk ATAC-seq log2(Unique Fragments)",
  rastr =T)#
p
#
# save Source Data
K562_ArchR.meta <- as.data.frame(K562_ArchR@cellColData)
table(K562_ArchR.meta$Sample)
#
cellsrename <- c(paste(rep("SPATACseq_SL1_cells_",798),1:798,sep = ""),
                 paste(rep("SPATACseq_SL2_cells_",826),1:826,sep = ""),
                 paste(rep("SPATACseq_SL3_cells_",1534),1:1534,sep = ""),
                 paste(rep("SPATACseq_SL4_cells_",1595),1:1595,sep = ""),
                 paste(rep("SPATACseq_SL5_cells_",562),1:562,sep = ""),
                 paste(rep("SPATACseq_SL6_cells_",755),1:755,sep = ""),
                 paste(rep("SPATACseq_SL7_cells_",1560),1:1560,sep = ""),
                 paste(rep("SPATACseq_SL8_cells_",1474),1:1474,sep = ""),
                 paste(rep("Jason_Fludigm_cells_",472),1:472,sep = ""),
                 paste(rep("Chen_Fludigm_cells_",181),1:181,sep = ""),
                 paste(rep("Chen_Plate_cells_",188),1:188,sep = ""))
rownames(K562_ArchR.meta) <- cellsrename
#
head(K562_ArchR.meta)
K562_ArchR.meta2 <- K562_ArchR.meta[,c(2,14,15,16)]
head(K562_ArchR.meta2)
write.csv(K562_ArchR.meta2,file = "SourceData/Fig1F.csv")
#
write.csv(cor(log2(Coverage.merge+1))[1:18,1:18],file = "SourceData/FigS1F.csv")
write.csv(cor(log2(Coverage.merge+1))[1:11,1:11],file = "SourceData/Fig1E.csv")
#
FigS1H <- cbind(as.numeric(ax.f.filter.log.filter[,"c"]), as.numeric(ax.f.filter.log.filter[,"d"]))
dim(FigS1H)
head(FigS1H)
colnames(FigS1H) <- c("SPATACseq_SL3","SPATACseq_SL4")
write.csv(FigS1H,file = "SourceData/FigS1H.csv")
#
Fig1D <- cbind(as.numeric(log2(SPATAC.SL8.cov.cpm+1)), as.numeric(log2(K562_rep1.cov.cpm+1)))
colnames(Fig1D) <- c("SPATACseq_SL8","K562_Bulk_ATACseq")
write.csv(Fig1D,file = "SourceData/Fig1D.csv")
#
# For Figure 2B
# chr19
# 35639001-35780000
# (35780000-3563900)/100 = 322161
#
rep("chr19",322161)
seq(35639001,(35780000-99),100)
seq(35639001+99,35780000,100)
#
head(seq(35639001,(35780000-99),100))
head(head(seq(35639001+99,35780000,100)))
#
tail(seq(35639001,(35780000-99),100))
tail(seq(35639001+99,35780000,100))
#
length(seq(35639001,(35780000-99),100))
length(seq(35639001+99,35780000,100))
#
target_bin <- cbind(rep("chr19",322161),
                    seq(35639001,(35780000-99),100),
                    seq(35639001+99,35780000,100))
head(target_bin)
write.table(target_bin,file = "target_bin.csv",quote = F, sep = "\t",col.names = F,row.names = F)
#
#
library(ChIPseeker)
DotLine_bin <- readPeakFile("target_bin.csv")
DotLine_bin
#
K562_ArchR <- addPeakSet(K562_ArchR,peakSet = DotLine_bin,force = T)
K562_ArchR <- addPeakMatrix(K562_ArchR,force = T)
plotGroups(ArchRProj = K562_ArchR, groupBy = "Sample", colorBy = "cellColData", name = "FRIP",plotAs = "violin",alpha = 0.4,addBoxPlot = TRUE)
#
Get.Peak.mat <- getMatrixFromProject(
  ArchRProj = K562_ArchR,
  useMatrix = "PeakMatrix")
Get.Peak.mat.matrix <- Get.Peak.mat@assays@data$PeakMatrix
####
peak_chr_name <- rep(as.character(Get.Peak.mat@rowRanges@seqnames@values),Get.Peak.mat@rowRanges@seqnames@lengths)
peak_start_site <- Get.Peak.mat@rowRanges@ranges@start
peak_end_site <- Get.Peak.mat@rowRanges@ranges@start + Get.Peak.mat@rowRanges@ranges@width
peak.names <- paste(peak_chr_name,peak_start_site,peak_end_site,sep = "-")
###
rownames(Get.Peak.mat.matrix) <- peak.names
head(Get.Peak.mat.matrix)
dim(Get.Peak.mat.matrix)
max(Get.Peak.mat.matrix)
head(Get.Peak.mat.matrix)
#
#
library(pheatmap)
#
Get.Peak.mat.matrix2 <- t(Get.Peak.mat.matrix)
Get.Peak.mat.matrix2 <- Get.Peak.mat.matrix2[1:100,]
dim(Get.Peak.mat.matrix2)
#
DotLine_bin_V <- c(0,0)
for (i in 1:100) {
  for (j in 1:322160) {
    if (Get.Peak.mat.matrix2[i,j] > 0) {
      DotLine_bin_V <- rbind(DotLine_bin_V,c(i,j))
    }
  }
}
plot(DotLine_bin_V[,1],DotLine_bin_V[,2])
#
dim(DotLine_bin_V)
head(DotLine_bin_V)
plot(DotLine_bin_V[,2],DotLine_bin_V[,1])
plot(DotLine_bin_V[,2],DotLine_bin_V[,1],pch=20,col ="red",cex=0.01)#######牛
#
pheatmap::pheatmap(Get.Peak.mat.matrix[1:30000,],cluster_rows = F,cluster_cols = F,show_rownames=F,
                   show_colnames = F)
pheatmap::pheatmap(Get.Peak.mat.matrix,cluster_rows = F,cluster_cols = F,show_rownames=F,
                   show_colnames = F)
##
Get.Peak.mat.matrix2 <- t(Get.Peak.mat.matrix)
Get.Peak.mat.matrix2[which(Get.Peak.mat.matrix2>1)] <- 1
pheatmap::pheatmap(Get.Peak.mat.matrix2[,1:30000],cluster_rows = F,cluster_cols = F,show_rownames=F,
                   show_colnames = F,color = c("white","red"))
pheatmap::pheatmap(Get.Peak.mat.matrix2,cluster_rows = F,cluster_cols = F,show_rownames=F,
                   show_colnames = F,color = c("white","red"))
#
#
##################################################################
##################################################################
##################################################################
#
setwd("/home/sunkeyong/data_storage/SPATAC/Figure1/hepa/")
library(Seurat)
library(ArchR)
addArchRGenome("mm10")
library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
addArchRThreads(threads = 13)
###
Bed_files= c("./mouse/SPATACseqSL1.sort.bam.fragments.filter.sort.bed.gz",
             "./mouse/SPATACseqSL2.sort.bam.fragments.filter.sort.bed.gz",
             "./mouse/SPATACseqSL3.sort.bam.fragments.filter.sort.bed.gz",
             "./mouse/SPATACseqSL4.sort.bam.fragments.filter.sort.bed.gz",
             "./mouse/SPATACseqSL5.sort.bam.fragments.filter.sort.bed.gz",
             "./mouse/SPATACseqSL6.sort.bam.fragments.filter.sort.bed.gz",
             "./mouse/SPATACseqSL7.sort.bam.fragments.filter.sort.bed.gz",
             "./mouse/SPATACseqSL8.sort.bam.fragments.filter.sort.bed.gz",
             "./mouse/Jason.sort.bed.gz",
             "./mouse/T_C1.sort.bed.gz",
             "./mouse/T_plate.sort.bed.gz")
names(Bed_files) = c("SPATACseqSL1","SPATACseqSL2","SPATACseqSL3","SPATACseqSL4",
                     "SPATACseqSL5","SPATACseqSL6","SPATACseqSL7","SPATACseqSL8",
                     "Jason","T_C1","T_plate")
#
ArrowFiles <- createArrowFiles(inputFiles = Bed_files,
                               sampleNames = names(Bed_files),
                               minTSS = 4,
                               minFrags = 500,
                               addTileMat = TRUE,
                               addGeneScoreMat = TRUE,
                               force = T)
#
Hepa_ArchR <- ArchRProject(ArrowFiles = ArrowFiles,
                           copyArrows = TRUE)
Hepa_ArchR
#
table(Hepa_ArchR$Sample)
#
Hepa_ArchR$ReadsInTSS.ratio <- Hepa_ArchR$ReadsInTSS/Hepa_ArchR$nFrags
Hepa_ArchR$ReadsInPromoter.ratio <- Hepa_ArchR$ReadsInPromoter/Hepa_ArchR$nFrags
Hepa_ArchR$nFrags.log2 <- log2(Hepa_ArchR$nFrags)
####### Quality plot 
plotGroups(ArchRProj = Hepa_ArchR, groupBy = "Sample", colorBy = "cellColData", name = "TSSEnrichment",plotAs = "violin",alpha = 0.4,addBoxPlot = TRUE)
plotGroups(ArchRProj = Hepa_ArchR, groupBy = "Sample", colorBy = "cellColData", name = "nFrags.log2",plotAs = "violin",alpha = 0.4,addBoxPlot = TRUE)
plotGroups(ArchRProj = Hepa_ArchR, groupBy = "Sample", colorBy = "cellColData", name = "ReadsInTSS.ratio",plotAs = "violin",alpha = 0.4,addBoxPlot = TRUE)
plotGroups(ArchRProj = Hepa_ArchR, groupBy = "Sample", colorBy = "cellColData", name = "ReadsInPromoter.ratio",plotAs = "violin",alpha = 0.4,addBoxPlot = TRUE)
#
#### plot DNA fragment distribution and TSS enrichment
#
p1 <- plotFragmentSizes(ArchRProj = Hepa_ArchR,returnDF =F,maxSize = 1000)
p1
p2 <- plotTSSEnrichment(ArchRProj = Hepa_ArchR)
p2
idxPass <- which(Hepa_ArchR$Sample=="SPATACseqSL8" | Hepa_ArchR$Sample=="Jason" | Hepa_ArchR$Sample=="T_C1" | Hepa_ArchR$Sample=="T_plate")
cellsPass <- Hepa_ArchR$cellNames[idxPass]
Hepa_ArchR.sub <- Hepa_ArchR[cellsPass, ]
p3 <- plotFragmentSizes(ArchRProj = Hepa_ArchR.sub)
p3
p4 <- plotTSSEnrichment(ArchRProj = Hepa_ArchR.sub)
p4
#
#
setwd("/home/sunkeyong/SPATAC/Figure1/Hepa/bedtools")
# for correlation analysis
# first calculated the coverage in the peakset using the following bash script 
#nohup bedtools coverage -a Hepa_peak.bed -b SL1.bed > SL1.coverage.bed &
# then load all coverage file
SL1 <- read.table("SL1.coverage.bed")
SL2 <- read.table("SL2.coverage.bed")
SL3 <- read.table("SL3.coverage.bed")
SL4 <- read.table("SL4.coverage.bed")
SL5 <- read.table("SL5.coverage.bed")
SL6 <- read.table("SL6.coverage.bed")
SL7 <- read.table("SL7.coverage.bed")
SL8 <- read.table("SL8.coverage.bed")
Hepa_rep1 <- read.table("Hepa_rep1.coverage.bed")
Hepa_rep2 <- read.table("Hepa_rep2.coverage.bed")
#
# other.seq for negative control
other.seq <- cbind(sample(c(rep(0,dim(SL8)[1]),rep(1,dim(SL8)[1])),dim(SL8)[1]),
                   sample(c(rep(0,dim(SL8)[1]),rep(1,dim(SL8)[1])),dim(SL8)[1]))
colnames(other.seq) <- c("other1","other2")
head(other.seq)
dim(other.seq)
#
Coverage.merge <- cbind(SL1$V4,SL3$V4,SL3$V4,SL4$V4,SL5$V4,SL6$V4,SL7$V4,SL8$V4,
                        Hepa_rep1$V4,Hepa_rep2$V4,other.seq)
#
head(Coverage.merge)
dim(Coverage.merge)
dim(Coverage.merge)
colnames(Coverage.merge) <- c("SL1","SL3","SL3","SL4","SL5","SL6","SL7","SL8",
                              "Hepa_rep1","Hepa_rep2","other1","other2")
head(Coverage.merge)
#
# For Figure S1I
library(pheatmap)
pheatmap(cor(Coverage.merge))
pheatmap(cor(log2(Coverage.merge+1)))
pheatmap(cor(log2(Coverage.merge+1)),cluster_rows = T,cluster_cols = T,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(300))
pheatmap(cor(log2(Coverage.merge+1)),cluster_rows = F,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(300))
#

##########################
##########################
# analyze doublets
SL1.meta <- read.csv("/home/sunkeyong/data_storage/MOPA_Mouse_Embryo_project/SPATAC-seq_Doublets/processed_data_file/SL1_Cell_Ranger_metadata.csv")
SL2.meta <- read.csv("/home/sunkeyong/data_storage/MOPA_Mouse_Embryo_project/SPATAC-seq_Doublets/processed_data_file/SL2_Cell_Ranger_metadata.csv")
SL3.meta <- read.csv("/home/sunkeyong/data_storage/MOPA_Mouse_Embryo_project/SPATAC-seq_Doublets/processed_data_file/SL3_Cell_Ranger_metadata.csv")
SL4.meta <- read.csv("/home/sunkeyong/data_storage/MOPA_Mouse_Embryo_project/SPATAC-seq_Doublets/processed_data_file/SL4_Cell_Ranger_metadata.csv")
SL5.meta <- read.csv("/home/sunkeyong/data_storage/MOPA_Mouse_Embryo_project/SPATAC-seq_Doublets/processed_data_file/SL5_Cell_Ranger_metadata.csv")
SL6.meta <- read.csv("/home/sunkeyong/data_storage/MOPA_Mouse_Embryo_project/SPATAC-seq_Doublets/processed_data_file/SL6_Cell_Ranger_metadata.csv")
SL7.meta <- read.csv("/home/sunkeyong/data_storage/MOPA_Mouse_Embryo_project/SPATAC-seq_Doublets/processed_data_file/SL7_Cell_Ranger_metadata.csv")
SL8.meta <- read.csv("/home/sunkeyong/data_storage/MOPA_Mouse_Embryo_project/SPATAC-seq_Doublets/processed_data_file/SL8_Cell_Ranger_metadata.csv")
#
SL1.meta.r <- SL1.meta
SL1.meta.r$Debris <- "Debris"
dim(SL1.meta.r)
SL1.meta.r[which(SL1.meta.r$passed_filters_GRCh38 >= 500 | SL1.meta.r$passed_filters_mm10 >= 500),24] <- "Cell"
table(SL1.meta.r$Debris)
SL1.meta.r$Fra <- SL1.meta.r$passed_filters_GRCh38+SL1.meta.r$passed_filters_mm10
SL1.meta.r$h.m.ratio <- SL1.meta.r$passed_filters_GRCh38/SL1.meta.r$Fra
SL1.meta.r$species <- SL1.meta.r$h.m.ratio
dim(SL1.meta.r)
SL1.meta.r[which(SL1.meta.r$h.m.ratio >= 0.8 & SL1.meta.r$Debris == "Cell"),27] <- "Human"
SL1.meta.r[which(SL1.meta.r$h.m.ratio <= 0.2 & SL1.meta.r$Debris == "Cell"),27] <- "Mouse"
SL1.meta.r[which(SL1.meta.r$h.m.ratio < 0.8 & SL1.meta.r$h.m.ratio > 0.2 & SL1.meta.r$Debris == "Cell"),27] <- "Doublets"
SL1.meta.r[which(SL1.meta.r$Debris == "Debris"),27] <- "Debris"
table(SL1.meta.r$species)
1-sum(subset(SL1.meta.r,species=="Debris")$total)/sum(SL1.meta.r$total)
2*table(SL1.meta.r$species)[2]/sum(table(SL1.meta.r$species)[2],table(SL1.meta.r$species)[3],table(SL1.meta.r$species)[4])*100
#
SL1.meta.r.f <- SL1.meta.r
SL1.meta.r.f[which(SL1.meta.r$passed_filters_GRCh38 > 30000),12] <- 30000
SL1.meta.r.f[which(SL1.meta.r$passed_filters_mm10 > 30000),14] <- 30000
ggplot(data=SL1.meta.r.f,aes(x=passed_filters_GRCh38/10000,y=passed_filters_mm10/10000,color=species))+
  geom_point()+theme_test()+labs(x = "hg38 unique fragments", y = "mm10 unique fragments")+
  scale_color_manual(values = c("gray","#7CAF33","#C77CF8","#58C0C4"))
####
ggplot(data=SL1.meta.r,aes(x=log2(passed_filters_GRCh38+1),y=log2(passed_filters_mm10+1),color=species))+
  geom_point()+theme_test()+labs(x = "hg38 unique fragments", y = "mm10 unique fragments")+
  scale_color_manual(values = c("gray","#7CAF33","#C77CF8","#58C0C4"))
####
##
# for Figure 1G
# read each metadata and loading the ratio of fragments in cells
# and loading the ratio of fragment in debris.
boxplot(c(75.17016,74.29343,74.27785,74.06365,72.87156,73.82872,72.75904,72.24117),# the ratio of fragments in cells
        c(1.831552,1.813454,2.090286,1.948361,2.818903,2.196699,2.523711,2.341398)) # the ratio of fragment in debris.
#
#
