# Figure 4
# analyse the peak
setwd("/home/sunkeyong/ZEPA/Peakset")
#
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
# for the Figure 4A and 4B
# perform peak calling in each stage using macs2
# Taking 10 hpf as an example
# in the environment of the Figure 2
# Making Pseudo-bulk Replicates
addArchRThreads(threads = 1) 
zhpf10 <- addGroupCoverages(ArchRProj = zhpf10, groupBy = "ID")
# Peak calling
zhpf10 <- addReproduciblePeakSet(
  ArchRProj = zhpf10, 
  groupBy = "ID",
  genomeSize = 1.5e9,
  excludeChr = c("MT",paste("chr",chr_rm$V1,sep = "")), 
  genomeAnnotation = genomeAnnotation,
  geneAnnotation = GRCz11.geneAnnotation,
  pathToMacs2 = "/home/sunkeyong/.local/bin/macs2", force=T
)
# get peakset at 10 hpf and save it
z10hpf.peakset <- getPeakSet(zhpf10)
save(zhpf10.peakset,file="/home/sunkeyong/ZEPA/OBO/hpf10/zhpf10.peakset.RData")
#
# Here, we load all the peaksets of 20 stages, directly
load("/home/sunkeyong/ZEPA/OBO/hpf4/z4hpf.peakset.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf5/z5hpf.peakset.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf6/z6hpf.peakset.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf7/z7hpf.peakset.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf8/zhpf8.peakset.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf9/z9hpf.peakset.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf10/z10hpf.peakset.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf11/z11hpf.peakset.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf12/z12hpf.peakset.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf14/z14hpf.peakset.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf18/z18hpf.peakset.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf20/z20hpf.peakset.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf22/z22hpf.peakset.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf24/z24hpf.peakset.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf30/z30hpf.peakset.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf34/z34hpf.peakset.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf38/z38hpf.peakset.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf42/z42hpf.peakset.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf48/z48hpf.peakset.RData")
load("/home/sunkeyong/ZEPA/OBO/hpf72/z72hpf.peakset.RData")
load("/home/sunkeyong/ZEPA/All_cells/ZEPA.peakset.RData")
#
### For the Figure B
table(ZEPA.peakset$peakType)
info = c(370357, 444816, 85503, 58364)
# name
names = c("4.Distal","3.Intronic","2.Exonic","1.Promoter")
#
# color (options)
cols = c("#ED1C24","#22B14C","#FFC90E","#3f48CC")
# calculating the ratio
piepercent = paste(round(100*info/sum(info)), "%")
# map
pie(info, labels=piepercent, main = "zebrafish-all-peaktype", col=cols, family='GB1')
# add legend
legend("topright", names, cex=0.8, fill=cols)
###
library(ggplot2)
library(ggpubr)
peak.info <- as.data.frame(info)
peak.info$name <- names
peak.info
ggdonutchart(peak.info,"info", label = "name",fill = "name",
             color = "white",palette = c("#ED1C24","#22B14C","#FFC90E","#3f48CC"))
###
hpf04.peakset <- z4hpf.peakset
hpf05.peakset <- z5hpf.peakset
hpf06.peakset <- z6hpf.peakset
hpf07.peakset <- z7hpf.peakset
hpf08.peakset <- zhpf8.peakset
hpf09.peakset <- z9hpf.peakset
hpf10.peakset <- z10hpf.peakset
hpf11.peakset <- z11hpf.peakset
hpf12.peakset <- z12hpf.peakset
hpf14.peakset <- z14hpf.peakset
hpf18.peakset <- z18hpf.peakset
hpf20.peakset <- z20hpf.peakset
hpf22.peakset <- z22hpf.peakset
hpf24.peakset <- z24hpf.peakset
hpf30.peakset <- z30hpf.peakset
hpf34.peakset <- z34hpf.peakset
hpf38.peakset <- z38hpf.peakset
hpf42.peakset <- z42hpf.peakset
hpf48.peakset <- z48hpf.peakset
hpf72.peakset <- z72hpf.peakset
#
ZEPA.stage.peakset <- list(hpf04.peakset,hpf05.peakset,hpf06.peakset,hpf07.peakset,hpf08.peakset,
                           hpf09.peakset,hpf10.peakset,hpf11.peakset,hpf12.peakset,hpf14.peakset,
                           hpf18.peakset,hpf20.peakset,hpf22.peakset,hpf24.peakset,hpf30.peakset,
                           hpf34.peakset,hpf38.peakset,hpf42.peakset,hpf48.peakset,hpf72.peakset)
#
sample.list <- c("hpf04","hpf05","hpf06","hpf07","hpf08",
                 "hpf09","hpf10","hpf11","hpf12","hpf14",
                 "hpf18","hpf20","hpf22","hpf24","hpf30",
                 "hpf34","hpf38","hpf42","hpf48","hpf72")
#
zTime_peaktype.tmp <- c()
zTime_peaktype <- c()
for (i in 1:21) {
  sample.tmp <- ZEPA.stage.peakset[[i]]
  zTime_peaktype.tmp <- c(zTime_peaktype.tmp,as.data.frame(table(sample.tmp$peakType))$Freq)
}
#
zTime_peaktype <- cbind(rep(sample.list,each=4),rep(c("Distal","Exonic","Intronic","Promoter"),21),zTime_peaktype.tmp)
zTime_peaktype <- as.data.frame(zTime_peaktype)
colnames(zTime_peaktype) <- c("hpf","type","number")
zTime_peaktype$number <- as.numeric(zTime_peaktype$number)
zTime_peaktype$number.log2 <- log2(zTime_peaktype$number)
zTime_peaktype$number.log10 <- log10(zTime_peaktype$number)
zTime_peaktype
# plot the Figure 4A
ggplot(data=zTime_peaktype, aes(x=hpf, y=number, fill=type)) +
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("#459ACD",
                             "#A0CC59",
                             "#EB772F",
                             "#854C9E"))
#

# plot Figure 4G
# using the deeptools to calculate the conservation 
# here loadding the results of deeptools
ZEPA.peakset.Distal.phastCons <- read.table("/home/sunkeyong/ZEPA/Peakset/conservation/Distal.bed.phastCons.center.2")
ZEPA.peakset.Exon.phastCons <- read.table("/home/sunkeyong/ZEPA/Peakset/conservation/Exonic.bed.phastCons.center.2")
ZEPA.peakset.Inton.phastCons <- read.table("/home/sunkeyong/ZEPA/Peakset/conservation/Intronic.bed.phastCons.center.2")
ZEPA.peakset.Promoter.phastCons <- read.table("/home/sunkeyong/ZEPA/Peakset/conservation/Promoter.bed.phastCons.center.2")
Tile.Random.phastCons <- read.table("/home/sunkeyong/ZEPA/Peakset/conservation/Tile.Random.bed")
#
boxplot(ZEPA.peakset.Exon.phastCons$V7,
        ZEPA.peakset.Distal.phastCons$V7,
        ZEPA.peakset.Inton.phastCons$V7,
        ZEPA.peakset.Promoter.phastCons$V7,
        Tile.Random.phastCons$V7,
        outline =F,ylim = c(0, 1),
        main="ZEPA.peakset.conservation",names = c("Exon","Distal","Inton","Promoter","Shuffle"),ylab = "phastCons score")
#
wilcox.test(ZEPA.peakset.Promoter.phastCons$V7,Tile.Random.phastCons$V7)
#
# For Figure 4C
# using the bedtools to calculate the overlapped base pairs or peak number
# here loadding the results of bedtools
sample.list <- c("Peaks_sc","Peaks_bulk","Peaks_adult","Enhancer","TSS")
row.values <- c()
for (i in 1:length(sample.list)) {
  row.values.tmp <- c()
  for (j in 1:length(sample.list)) {
    bed.x <- read.table(file = paste("/home/sunkeyong/ZEPA/Peakset/overlap/",sample.list[i],".bed",sep = ""))
    row.values.tmp <- c(row.values.tmp,length(unique(paste(bed.x$V1,bed.x$V2,bed.x$V3))))
  }
  row.values <- cbind(row.values,row.values.tmp)
}
rownames(row.values) <- sample.list
colnames(row.values) <- sample.list
#
pheatmap(t(row.values)/sample.bp.length,cluster_rows = F,cluster_cols = F)
#
pheatmap(t(as.matrix(t(row.values)/sample.bp.length)),
         cluster_rows = F,cluster_cols = F,
         color = colorRampPalette(c(c("#FFFFFF","#4699C8")))(10))
#
# Figure 4D
# plot the cell ratio 
# plot the track views of the marker genes
# in the environment of hpf8 (Figure)
cluster.levels <- c("hpf8:Prechordal plate","hpf8:EVL","hpf8:Epidermal",
                    "hpf8:Neural plate posterior","hpf8:Neural plate anterior",
                    "hpf8:Tailbud spinal cord","hpf8:Tailbud mesoderm",
                    "hpf8:Mesoderm adaxial cells","hpf8:Mesoderm lateral plate",
                    "hpf8:Endoderm","hpf8:Notochord","hpf8:YSL")
color.hpf08h.track <- c("#6F94E6","#F5CFB0","#71CECB",
                        "#5C509D","#BDE1CE",
                        "#489657","#EDCCE2",
                        "#BDA6CC","#323690",
                        "#EF8632","#DAC386","gray")
#
hpf08.ratio <- as.data.frame(table(zhpf8$ID))
rownames(hpf08.ratio) <- hpf08.ratio$Var1
hpf08.ratio <- hpf08.ratio[cluster.levels,]
barplot(hpf08.ratio$Freq/sum(hpf08.ratio$Freq),col =color.hpf08h.track )
barplot(hpf08.ratio$Freq/sum(hpf08.ratio$Freq),col =color.hpf08h.track,space =0 )
#
Figure4D <- cbind(cluster.levels,hpf08.ratio$Freq/sum(hpf08.ratio$Freq))
Figure4D
write.csv(Figure4D,file = "/home/sunkeyong/ZEPA/SourceData/Fig4D.csv")
#
p <- plotBrowserTrack(ArchRProj = zhpf8, groupBy = "ID", geneSymbol = "klf2b", plotSummary = c("bulkTrack", "featureTrack", "geneTrack"),
                      sizes = c(10, 1, 1),minCells =1,ylim=c(0.000001,0.99999),useGroups=cluster.f,
                      pal= color.hpf08h.track,
                      upstream = 4000,downstream = 10000)
grid::grid.newpage()
grid::grid.draw(p$klf2b)
plotPDF(plotList = p, name = "klf2b.pdf", ArchRProj = zhpf8, addDOC = FALSE, width = 5, height = 8)
#
p <- plotBrowserTrack(ArchRProj = zhpf8, groupBy = "ID", geneSymbol = "sox32", plotSummary = c("bulkTrack", "featureTrack","geneTrack"),
                      sizes = c(10,1, 1),minCells =1,ylim=c(0,1),useGroups=cluster.f,
                      pal= color.hpf08h.track,
                      upstream = 5000,downstream = 2000)
grid::grid.newpage()
grid::grid.draw(p$sox32)
plotPDF(plotList = p, name = "sox32.pdf", ArchRProj = zhpf8, addDOC = FALSE, width = 5, height = 8)
#
#
p <- plotBrowserTrack(ArchRProj = zhpf8, groupBy = "ID", geneSymbol = "krt92", plotSummary = c("bulkTrack", "featureTrack", "geneTrack"),
                      sizes = c(10, 1, 1),minCells =1,ylim=c(0,1),useGroups=cluster.f,
                      pal= color.hpf08h.track,
                      upstream = 5000,downstream = 19000)
grid::grid.newpage()
grid::grid.draw(p$krt92)
plotPDF(plotList = p, name = "krt92.pdf", ArchRProj = zhpf8, addDOC = FALSE, width = 5, height = 8)
#
p <- plotBrowserTrack(ArchRProj = zhpf8, groupBy = "ID", geneSymbol = "hesx1", plotSummary = c("bulkTrack", "featureTrack", "geneTrack"),
                      sizes = c(10, 1, 1),minCells =1,ylim=c(0,1),useGroups=cluster.f,
                      pal= color.hpf08h.track,
                      upstream = 2000,downstream = 3500)
grid::grid.newpage()
grid::grid.draw(p$hesx1)
plotPDF(plotList = p, name = "hesx1.pdf", ArchRProj = zhpf8, addDOC = FALSE, width = 5, height = 8)
#
#
p <- plotBrowserTrack(ArchRProj = zhpf8, groupBy = "ID", geneSymbol = "twist2", plotSummary = c("bulkTrack", "featureTrack", "geneTrack"),
                      sizes = c(10, 1, 1),minCells =1,ylim=c(0,1),useGroups=cluster.f,
                      pal= color.hpf08h.track,
                      upstream = 2500,downstream = 10000)
grid::grid.newpage()
grid::grid.draw(p$twist2)
plotPDF(plotList = p, name = "twist2.pdf", ArchRProj = zhpf8, addDOC = FALSE, width = 5, height = 8)
#
p <- plotBrowserTrack(ArchRProj = zhpf8, groupBy = "ID", geneSymbol = "frzb", plotSummary = c("bulkTrack", "featureTrack", "geneTrack"),
                      sizes = c(10, 1, 1),minCells =1,ylim=c(0,1),useGroups=cluster.f,
                      pal= color.hpf08h.track,
                      upstream = 2000,downstream = 10000)
grid::grid.newpage()
grid::grid.draw(p$frzb)
plotPDF(plotList = p, name = "frzb.pdf", ArchRProj = zhpf8, addDOC = FALSE, width = 5, height = 8)
#
#
p <- plotBrowserTrack(ArchRProj = zhpf8, groupBy = "ID", geneSymbol = "mir206-1", plotSummary = c("bulkTrack", "featureTrack", "geneTrack"),
                      sizes = c(10, 1, 1),minCells =1,ylim=c(0,1),useGroups=cluster.f,
                      pal= color.hpf08h.track,
                      upstream = 4000,downstream = 3000)
grid::grid.newpage()
grid::grid.draw(p$`mir206-1`)
plotPDF(plotList = p, name = "mir206.1.pdf", ArchRProj = zhpf8, addDOC = FALSE, width = 5, height = 8)
#
p <- plotBrowserTrack(ArchRProj = zhpf8, groupBy = "ID", geneSymbol = "actb1", plotSummary = c("bulkTrack", "featureTrack", "geneTrack"),
                      sizes = c(10, 1, 1),minCells =1,ylim=c(0,1),useGroups=cluster.f,
                      pal= color.hpf08h.track,
                      upstream = 3800,downstream = 2000)
grid::grid.newpage()
grid::grid.draw(p$actb1)
plotPDF(plotList = p, name = "actb1.pdf", ArchRProj = zhpf8, addDOC = FALSE, width = 5, height = 8)
#
#
#########
# for the Figure 4E and 4F
# calculating the odds ratio in different chromHMM regions
#  taking 24 hpf as an example
library(rtracklayer)
options(scipen = 200)
load("/home/sunkeyong/ZEPA/OBO/hpf24/z24hpf.peakset.RData")
hpf24.chromHMMraw <- import("/home/sunkeyong/ZEPA/DCC/raw_file_from_UCSC/Prim5_ChromHMM.bb")
#
hpf24.chromHMM.metadata <- as.data.frame(hpf24.chromHMMraw)
dim(hpf24.chromHMM.metadata)
hpf24.chromHMM.metadata <- hpf24.chromHMM.metadata[1:(which(hpf24.chromHMM.metadata$seqnames=="chrM")[1]-1),]
dim(hpf24.chromHMM.metadata)
hpf24.chromHMM.metadata$seqnames <- as.character(hpf24.chromHMM.metadata$seqnames)
chromHMM.type <- as.data.frame(table(hpf24.chromHMM.metadata$name))
chromHMM.type[11,] <- chromHMM.type[2,]
chromHMM.type <- chromHMM.type[-2,] 
chromHMM.type$Var1 <- as.character(chromHMM.type$Var1)
rownames(chromHMM.type) <- 1:10
chromHMM.type
chromHMM.type$region <- 1
chromHMM.type$bp <- 1
#
hpf24.chromHMM <- hpf24.chromHMM.metadata
for (i in 1:10) {
  hpf24.chromHMM.x <- hpf24.chromHMM[which(hpf24.chromHMM$name==chromHMM.type$Var1[i]),]
  hpf24.chromHMM.x <- hpf24.chromHMM.x[,1:3]
  chromHMM.type$region[i] <- nrow(hpf24.chromHMM.x)
  chromHMM.type$bp[i] <- sum(abs(hpf24.chromHMM.x$end-hpf24.chromHMM.x$start))+nrow(hpf24.chromHMM.x)
  write.table(hpf24.chromHMM.x,file = paste("/home/sunkeyong/ZEPA/DCC/cal_chromHMM/hpf24/region/E",i,".csv",sep = ""),quote=F, sep = "\t",row.names = F,col.names = F)
}
#
head(chromHMM.type)
#
z24hpf.peakset.bed <- cbind(as.character(z24hpf.peakset@seqnames),
                          z24hpf.peakset@ranges@start,
                          z24hpf.peakset@ranges@start+500)
head(z24hpf.peakset.bed)
write.table(z24hpf.peakset.bed,file = "/home/sunkeyong/ZEPA/DCC/cal_chromHMM/hpf24/region/z24hpf.peakset.bed",quote = F, sep = "\t",col.names = F,row.names = F)
#
#using the bedtools to calculated the overlapped base pair and peak number
#
chromHMM.type$O.region <- 1
chromHMM.type$O.bp <- 1
for (i in 1:10) {
  hpf24.Self.Ex.bed.tmp <- read.table(paste("/home/sunkeyong/ZEPA/DCC/cal_chromHMM/hpf24/region/Self.E",i,".bed",sep = ""))
  hpf24.Self.Ex.wa.bed.tmp <- read.table(paste("/home/sunkeyong/ZEPA/DCC/cal_chromHMM/hpf24/region/Self.E",i,".wa.bed",sep = ""))
  chromHMM.type$O.region[i] <- length(unique(paste(hpf24.Self.Ex.wa.bed.tmp$V1,hpf24.Self.Ex.wa.bed.tmp$V2,hpf24.Self.Ex.wa.bed.tmp$V3)))
  chromHMM.type$O.bp[i] <- sum(abs(hpf24.Self.Ex.bed.tmp$V2-hpf24.Self.Ex.bed.tmp$V3))
}
chromHMM.type$O.region.test <- 1
chromHMM.type$O.region.test.p <- 1
chromHMM.type$O.bp.test <- 1
chromHMM.type$O.bp.test.p <- 1
#
for (i in 1:10) {
  chromHMM.type_O.region.test.tmp <- fisher.test(matrix(c(chromHMM.type$O.region[i],
                                                          sum(chromHMM.type$O.region),
                                                          chromHMM.type$region[i],
                                                          sum(chromHMM.type$region)),
                                                        ncol=2,dimnames = list(c('Male','Female'),c('Placebo','Treated'))))
  chromHMM.type$O.region.test.p[i] <- chromHMM.type_O.region.test.tmp$p.value
  chromHMM.type$O.region.test[i] <- as.numeric(chromHMM.type_O.region.test.tmp$estimate)
  #
  chromHMM.type_O.bp.test.tmp <- fisher.test(matrix(c(chromHMM.type$O.bp[i],
                                                      sum(chromHMM.type$O.bp),
                                                      (chromHMM.type$bp[i]-chromHMM.type$O.bp[i]),
                                                      (sum(chromHMM.type$bp)-chromHMM.type$bp[i]+chromHMM.type$O.bp[i])),
                                                    ncol=2,dimnames = list(c('Male','Female'),c('Placebo','Treated'))))
  chromHMM.type$O.bp.test.p[i] <- chromHMM.type_O.bp.test.tmp$p.value
  chromHMM.type$O.bp.test[i] <- as.numeric(chromHMM.type_O.bp.test.tmp$estimate)
}
#
chromHMM.type.hpf24 <- chromHMM.type
save(chromHMM.type.hpf24,file = "chromHMM.type.hpf24.RData")
head(chromHMM.type.hpf24)
chromHMM.type
#
rownames(chromHMM.type.hpf24) <- chromHMM.type.hpf24$Var1
#
chromHMM.type.hpf24$O.bp.test <- as.numeric(chromHMM.type.hpf24$O.bp.test)
chromHMM.type.hpf24$O.bp.test.p <- as.numeric(chromHMM.type.hpf24$O.bp.test.p)
chromHMM.type.hpf24$O.region.test <- as.numeric(chromHMM.type.hpf24$O.region.test)
chromHMM.type.hpf24$O.region.test.p <- as.numeric(chromHMM.type.hpf24$O.region.test.p)
barplot(chromHMM.type.hpf24$O.bp.test)
barplot(chromHMM.type.hpf24$O.region.test)
save(chromHMM.type.hpf24,file = "/home/sunkeyong/ZEPA/DCC/cal_chromHMM/chromHMM.type.hpf24.RData")
# repeat above procedures for the data at 4 hpf, 8 hpf, 12 hpf, 24 hpf, 48 hpf
# 
# then loading all results
load("/home/sunkeyong/ZEPA/DCC/cal_chromHMM/chromHMM.type.hpf4.RData")
load("/home/sunkeyong/ZEPA/DCC/cal_chromHMM/chromHMM.type.hpf8.RData")
load("/home/sunkeyong/ZEPA/DCC/cal_chromHMM/chromHMM.type.hpf12.RData")
load("/home/sunkeyong/ZEPA/DCC/cal_chromHMM/chromHMM.type.hpf24.RData")
load("/home/sunkeyong/ZEPA/DCC/cal_chromHMM/chromHMM.type.hpf48.RData")
#
hpf4.PADREs <- import("/home/sunkeyong/ZEPA/DCC/raw_file_from_UCSC/Dome_PADREs.bb")
hpf8.PADREs <- import("/home/sunkeyong/ZEPA/DCC/raw_file_from_UCSC/Epi75_PADREs.bb")
hpf12.PADREs <- import("/home/sunkeyong/ZEPA/DCC/raw_file_from_UCSC/Hpf12_PADREs.bb")
hpf24.PADREs <- import("/home/sunkeyong/ZEPA/DCC/raw_file_from_UCSC/Prim5_PADREs.bb")
hpf48.PADREs <- import("/home/sunkeyong/ZEPA/DCC/raw_file_from_UCSC/LongPec_PADREs.bb")
#
zPADREs <- as.data.frame(cbind(table(hpf4.PADREs$name),
                               table(hpf8.PADREs$name),
                               table(hpf12.PADREs$name),
                               table(hpf24.PADREs$name),
                               table(hpf48.PADREs$name)))
zPADREs[11,] <- zPADREs[2,]
zPADREs <- zPADREs[-2,]
rownames(zPADREs) <- chromHMM.type.hpf4$Var1
colnames(zPADREs) <- c("Dome","Epi75","Hpf12","Prim5","LongPec")
head(zPADREs)
#
library(reshape2)
zPADREs <- melt(zPADREs)
zPADREs$type <- paste("C",1:10,sep = "")
zPADREs$type
head(zPADREs)
ggplot(zPADREs, aes(x=type,y=value, fill = variable))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_x_discrete(limits=paste("C",1:10,sep = ""))+
  scale_fill_manual(values=New.ZtimeColor[c(1,5,9,14,19)])+
  theme_classic()
#
barplot(c(zPADREs[which(zPADREs$variable=="Dome"),]$value[10]/(sum(zPADREs[which(zPADREs$variable=="Dome"),]$value)),
          zPADREs[which(zPADREs$variable=="Epi75"),]$value[10]/(sum(zPADREs[which(zPADREs$variable=="Epi75"),]$value)),
          zPADREs[which(zPADREs$variable=="Hpf12"),]$value[10]/(sum(zPADREs[which(zPADREs$variable=="Hpf12"),]$value)),
          zPADREs[which(zPADREs$variable=="Prim5"),]$value[10]/(sum(zPADREs[which(zPADREs$variable=="Prim5"),]$value)),
          zPADREs[which(zPADREs$variable=="LongPec"),]$value[10]/(sum(zPADREs[which(zPADREs$variable=="LongPec"),]$value))))
#
barplot(chromHMM.type.hpf4$O.region)
#
chromHMM.type.hpf4$hpf <- "hpf04"
chromHMM.type.hpf8$hpf <- "hpf08"
chromHMM.type.hpf12$hpf <- "hpf12"
chromHMM.type.hpf24$hpf <- "hpf24"
chromHMM.type.hpf48$hpf <- "hpf48"
#
chromHMM.type.merge <- rbind(chromHMM.type.hpf4,
                             chromHMM.type.hpf8,
                             chromHMM.type.hpf12,
                             chromHMM.type.hpf24,
                             chromHMM.type.hpf48)
chromHMM.type.merge <- as.data.frame(chromHMM.type.merge)
#
chromHMM.type.merge$O.bp.test.p
-log10(chromHMM.type.merge$O.bp.test)
#
darkblue <- "#303B7F"
darkred <- "#D51113"
yellow <- "#EECA1F"
#
chromHMM.type.merge$Var1 <- as.character(chromHMM.type.merge$Var1)
levels(chromHMM.type.merge$Var1)
levels(chromHMM.type.merge$Var1) <- rev(chromHMM.type.hpf4$Var1)
#
ggplot(chromHMM.type.merge, aes(x=Var1,y=O.region, fill = hpf))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_x_discrete(limits=as.character(chromHMM.type.hpf4$Var1))+
  scale_fill_manual(values=New.ZtimeColor[c(1,5,9,14,19)])+
  theme_classic()
#
ggplot(chromHMM.type.merge, aes(x=Var1,y=O.bp, fill = hpf))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_x_discrete(limits=as.character(chromHMM.type.hpf4$Var1))+
  scale_fill_manual(values=New.ZtimeColor[c(1,5,9,14,19)])+
  theme_classic()
#
#
ggplot(chromHMM.type.merge, aes(x=Var1,y=O.bp.test, fill = hpf))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_x_discrete(limits=as.character(chromHMM.type.hpf4$Var1))+
  scale_fill_manual(values=New.ZtimeColor[c(1,5,9,14,19)])+
  theme_classic()
#
ggplot(chromHMM.type.merge, aes(x=Var1,y=O.region.test, fill = hpf))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_x_discrete(limits=as.character(chromHMM.type.hpf4$Var1))+
  scale_fill_manual(values=New.ZtimeColor[c(1,5,9,14,19)])+
  theme_classic()
#
barplot(c(chromHMM.type.merge[which(chromHMM.type.merge$hpf=="hpf04"),]$O.region[10]/(sum(chromHMM.type.merge[which(chromHMM.type.merge$hpf=="hpf04"),]$O.region)),
          chromHMM.type.merge[which(chromHMM.type.merge$hpf=="hpf08"),]$O.region[10]/(sum(chromHMM.type.merge[which(chromHMM.type.merge$hpf=="hpf08"),]$O.region)),
          chromHMM.type.merge[which(chromHMM.type.merge$hpf=="hpf12"),]$O.region[10]/(sum(chromHMM.type.merge[which(chromHMM.type.merge$hpf=="hpf12"),]$O.region)),
          chromHMM.type.merge[which(chromHMM.type.merge$hpf=="hpf24"),]$O.region[10]/(sum(chromHMM.type.merge[which(chromHMM.type.merge$hpf=="hpf24"),]$O.region)),
          chromHMM.type.merge[which(chromHMM.type.merge$hpf=="hpf48"),]$O.region[10]/(sum(chromHMM.type.merge[which(chromHMM.type.merge$hpf=="hpf48"),]$O.region))))
#
my_palette <- colorRampPalette(c(darkblue,yellow,darkred), alpha=TRUE)(n=128)
ggplot(chromHMM.type.merge, aes(x=hpf,y=Var1)) +
  geom_point(aes(size=O.bp.test,color=O.bp.test
                 #,color=-log10(O.bp.test)
  )) +
  scale_color_gradientn('NES',
                        colors=my_palette) +
  scale_size_continuous(range = c(1,4)) + #scale of point sizes
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.x = element_text(angle = 45, size = 12, hjust = 0.3, vjust = 0.5, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title = element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"),
        legend.position = "bottom",
        plot.margin = unit(c(1,1,1,1), "lines"))+
  scale_y_discrete(limits=as.character(rev(chromHMM.type.hpf4$Var1)))
#
### analyse the bulk ATAC-seq and scATAC-seq data
# for Figure 4H
# loading the coverage results calculated by bedtools
cor.bulk.list <- c("4hpf","5hpf","6hpf","7hpf","8hpf","9hpf","10hpf",
                   "11hpf","12hpf","14hpf","18hpf","20hpf","22hpf",
                   "24hpf","30hpf","34hpf","38hpf",
                   "42hpf","48hpf","72hpf")
cor.bulk.value.tmp
cor.bulk.value <- c()
for (i in 1:21) {
  cor.bulk.value.tmp <- read.table(paste("/home/sunkeyong/ZEPA/Bulk/bulk_peak_based/",cor.bulk.list[i],".bed",sep = ""))
  cor.bulk.value <- cbind(cor.bulk.value,as.numeric(cor.bulk.value.tmp$V7))
}
#
cor.bulk.list2 <- paste(c("4hpf","5hpf","6hpf","7hpf","8hpf","9hpf","10hpf",
                          "11hpf","12hpf","14hpf","18hpf","20hpf","22hpf",
                          "24hpf","30hpf","34hpf","38hpf",
                          "42hpf","48hpf","72hpf"),"_sc",sep ="" )
cor.bulk.value <- as.data.frame(cor.bulk.value)
head(cor.bulk.value)
dim(cor.bulk.value)
colnames(cor.bulk.value) <- cor.bulk.list2
head(cor.bulk.value)
rownames(cor.bulk.value) <- paste(cor.bulk.value.tmp$V1,cor.bulk.value.tmp$V2,cor.bulk.value.tmp$V3,sep="-")
#
cor.bulk.list2.rank <- paste(c("4hpf","5hpf","6hpf","7hpf","8hpf","9hpf","10hpf",
                               "11hpf","12hpf","14hpf","18hpf","20hpf","22hpf",
                               "24hpf","30hpf","34hpf","38hpf",
                               "42hpf","48hpf","72hpf"),"_sc",sep ="" )
#
cor.bulk.value <- cor.bulk.value[,cor.bulk.list2.rank]
head(cor.bulk.value)
#
cor(cor.bulk.value)
#
pheatmap(cor(cor.bulk.value),cluster_cols = F,cluster_rows = F)
#
bulkATACseqvalues <- read.table("/home/sunkeyong/ZEPA/Bulk/bulk_peak_based/Bulkpeak.read.count.txt")
dim(bulkATACseqvalues)
bulkATACseqvalues2 <- bulkATACseqvalues[,7:46]
dim(bulkATACseqvalues2)
head(bulkATACseqvalues2)
#
cor.bulk.list3.rank <- paste(rep(c("4hpf","5hpf","6hpf","7hpf","8hpf","9hpf","10hpf",
                                   "11hpf","12hpf","14hpf","18hpf","20hpf","22hpf","24hpf","30hpf","34hpf","38hpf",
                                   "42hpf","48hpf","72hpf"),each=2),rep(c("_bulk1","_bulk2"),time=20),sep ="" )
#
colnames(bulkATACseqvalues2) <- cor.bulk.list3.rank
rownames(bulkATACseqvalues2) <- paste(cor.bulk.value.tmp$V1,cor.bulk.value.tmp$V2,cor.bulk.value.tmp$V3,sep="-")
head(bulkATACseqvalues2)
#
pheatmap(cor(bulkATACseqvalues2),cluster_cols = F,cluster_rows = F)
pheatmap(cor(bulkATACseqvalues2),cluster_cols = F,cluster_rows = F)
#
BPM <- cbind(bulkATACseqvalues2,cor.bulk.value)
pheatmap(cor(BPM),cluster_cols = F,cluster_rows = F)
pheatmap(cor(BPM),cluster_cols = T,cluster_rows = T)
#
zBulk_sc <- CreateSeuratObject(counts = BPM, project = "BPM", min.cells = 0, min.features = 0)
zBulk_sc
zBulk_sc <- NormalizeData(zBulk_sc, normalization.method = "RC", scale.factor = 1000000)
zBulk_sc <- FindVariableFeatures(zBulk_sc, selection.method = "vst", nfeatures = 30000)
all.genes <- rownames(zBulk_sc)
#zBulk_sc <- ScaleData(zBulk_sc, features = all.genes)
zBulk_sc <- RunPCA(zBulk_sc, features = VariableFeatures(object = zBulk_sc))
DimPlot(zBulk_sc, reduction = "pca",label = F)
zBulk_sc$sample <- Cells(zBulk_sc)
zBulk_sc$time <- unlist(strsplit(Cells(zBulk_sc),"_"))[seq(1,2*60,2)]
zBulk_sc$batch <- unlist(strsplit(Cells(zBulk_sc),"_"))[seq(2,2*60,2)]
DimPlot(zBulk_sc, reduction = "pca",label = T,repel = T,group.by = "sample",dims = c(1,2))+NoLegend()
DimPlot(zBulk_sc, reduction = "pca",label = T,repel = T,group.by = "sample",dims = c(1,3))+NoLegend()
#
zBulk_sc@active.ident <- as.factor(zBulk_sc$time)
levels(zBulk_sc) <- c("4hpf","5hpf","6hpf","7hpf","8hpf","9hpf","10hpf",
                  "11hpf","12hpf","14hpf","18hpf","20hpf","22hpf","24hpf",
                  "30hpf","34hpf","38hpf","42hpf","48hpf","72hpf")
ZEPA.timeColor <- c("#FAD6A4","#F7B264","#E1520D","#AD3417",
                    "#F7C7D2","#EA66A1","#E50C84",
                    "#D2E2F2","#227DBD","#004695",
                    "#CFE6C0","#A8D392","#0C9E4A","#006D32",
                    "#C6C4E0","#857BB8","#694799",
                    "#EFBF00","#B78510","#935913")
p <- DimPlot(zBulk_sc, reduction = "pca",group.by = "time",shape.by = "batch",pt.size =2,cols = ZEPA.timeColor,
             label=F,repel = T,dims = c(1,3))+NoLegend()+
  theme(plot.title = element_text(color="black", size=14))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 1))+
  theme(panel.border = element_blank())
ggsave("/home/sunkeyong/ZEPA/Bulk/Z.sc.bulk.pca.png",plot=p,device="png",dpi=400,units = "cm",width = 10,height = 10)
#
# for the figure S10D
BPM.x <- BPM[,c(seq(1,2*20,2),41:60)]
BPM.x.cor <- cor(BPM.x)[1:20,21:40]
pheatmap(BPM.x.cor,cluster_cols = F,cluster_rows = F)
col2 = colorRampPalette(rev(c('#67001F', '#B2182B', '#D6604D', '#F4A582',
                              '#FDDBC7', '#FFFFFF', '#D1E5F0', '#92C5DE',
                              '#4393C3', '#2166AC', '#053061')))
Heatmap(BPM.x.cor,
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = T,
        show_row_names = T,
        gap = unit(0.5, "mm"), 
        border = "black",
        col =col2(200))
# for the figure S2L
BPM.y <- BPM[,1:40]
BPM.y.cor <- cor(BPM.y)[seq(1,2*20,2),seq(2,2*20,2)]
pheatmap(BPM.y.cor,cluster_cols = F,cluster_rows = F)
Heatmap(BPM.y.cor,
        cluster_rows = F,
        cluster_columns = F,
        show_column_names = T,
        show_row_names = T,
        gap = unit(0.5, "mm"), 
        border = "black",
        col =col2(200))
#for the figure S10E
boxplot(diag(BPM.y.cor[1:20,1:20]),
        diag(BPM.y.cor[2:20,1:19]),
        diag(BPM.y.cor[3:20,1:18]),
        diag(BPM.y.cor[4:20,1:17]),
        diag(BPM.y.cor[5:20,1:16]),
        diag(BPM.y.cor[6:20,1:15]),
        diag(BPM.y.cor[7:20,1:14]),
        diag(BPM.y.cor[8:20,1:13]),
        diag(BPM.y.cor[9:20,1:12]),
        diag(BPM.y.cor[10:20,1:11]),
        diag(BPM.y.cor[11:20,1:10]),
        diag(BPM.y.cor[12:20,1:9]),
        diag(BPM.y.cor[13:20,1:8]),
        diag(BPM.y.cor[14:20,1:7]),
        diag(BPM.y.cor[15:20,1:6]),
        diag(BPM.y.cor[16:20,1:5]),
        diag(BPM.y.cor[17:20,1:4]),
        diag(BPM.y.cor[18:20,1:3]),
        diag(BPM.y.cor[19:20,1:2]),
        BPM.y.cor[20:20,1:1],outline=F,col=New.ZtimeColor 
)
#
# for the Figure 4I
# plot the trackview near sp5l using scATAC-seq in archr
# plot the trackview near sp5l using bulk ATAC-seq by IGV
# in the R environment of ZEPA all cells 
ZEPA.timeColor <- c("#FAD6A4","#F7B264","#E1520D","#AD3417",
                    "#F7C7D2","#EA66A1","#E50C84",
                    "#D2E2F2","#227DBD","#004695",
                    "#CFE6C0","#A8D392","#0C9E4A","#006D32",
                    "#C6C4E0","#857BB8","#694799",
                    "#EFBF00","#B78510","#935913")
p <- plotBrowserTrack(ArchRProj = ZEPA, groupBy = "time",pal=ZEPA.timeColor,
                      geneSymbol = "sp5l", upstream = 10000,downstream = 10000,
                      ylim = c(0.001,0.999999))
grid::grid.newpage()
grid::grid.draw(p$sp5l)
plotPDF(plotList = p, name = "sp5l.pdf", ArchRProj = ZEPA, addDOC = FALSE, width = 5, height = 8)
#
################
# analyse the temporal peaks in the Figure 4J and 4K
# loading the coverage results calculated by bedtools
cor.bulk.list <- c("4hpf","5hpf","6hpf","7hpf","8hpf","9hpf","10hpf",
                   "11hpf","12hpf","14hpf","18hpf","20hpf","22hpf",
                   "24hpf","30hpf","34hpf","38hpf",
                   "42hpf","48hpf","72hpf")
cor.bulk.value.tmp
cor.bulk.value <- c()
for (i in 1:21) {
  cor.bulk.value.tmp <- read.table(paste("/home/sunkeyong/ZEPA/Bulk/bulk_peak_based/",cor.bulk.list[i],".bed",sep = ""))
  cor.bulk.value <- cbind(cor.bulk.value,as.numeric(cor.bulk.value.tmp$V7))
}
#
cor.bulk.list2 <- paste(c("4hpf","5hpf","6hpf","7hpf","8hpf","9hpf","10hpf",
                          "11hpf","12hpf","14hpf","18hpf","20hpf","22hpf",
                          "24hpf","30hpf","34hpf","38hpf",
                          "42hpf","48hpf","72hpf"),"_sc",sep ="" )
cor.bulk.value <- as.data.frame(cor.bulk.value)
head(cor.bulk.value)
dim(cor.bulk.value)
colnames(cor.bulk.value) <- cor.bulk.list2
head(cor.bulk.value)
rownames(cor.bulk.value) <- paste(cor.bulk.value.tmp$V1,cor.bulk.value.tmp$V2,cor.bulk.value.tmp$V3,sep="-")
#
cor.bulk.list2.rank <- paste(c("4hpf","5hpf","6hpf","7hpf","8hpf","9hpf","10hpf",
                               "11hpf","12hpf","14hpf","18hpf","20hpf","22hpf",
                               "24hpf","30hpf","34hpf","38hpf",
                               "42hpf","48hpf","72hpf"),"_sc",sep ="" )
#
cor.bulk.value <- cor.bulk.value[,cor.bulk.list2.rank]
head(cor.bulk.value)
#
pheatmap(cor.bulk.value[1:10,],cluster_rows = F,cluster_cols = F)
#
axsum <- cor.bulk.value 
dim(axsum)
hist(log2(rowSums(axsum.filter+1)))
axsum.filter <- axsum
dim(axsum.filter)
axsum.filter <- axsum.filter[(which((rowSums(axsum.filter) > 1))),]
dim(axsum.filter)
hist(rowSums(axsum.filter))
axsum.filter <- axsum.filter[(which((rowSums(axsum.filter) < 400))),]
dim(axsum.filter)
hist(rowSums(axsum.filter))
#
axsum.filter.row.max <- apply(axsum.filter , 1, max)
df <- axsum.filter/axsum.filter.row.max
dim(df)
hist(apply(df , 1, sd))
#
set.seed(123)
km_result <- kmeans(df,30,nstart = 24)
#
dim(df)
#
library(pheatmap)
axsum.filter <- axsum.filter[rownames(df),]
dd <- cbind(axsum.filter,cluster=km_result$cluster)
head(dd)
dd <- as.data.frame(dd)
min(table(dd$cluster))
table(dd$cluster)
barplot(log2(table(dd$cluster)))
#
# pick one cluster to show
df.1 <- df[rownames(subset(dd,dd$cluster== "4")),]
df.1.sub <- df.1[sample(dim(df.1)[1],200),]
head(df.1.sub)
pheatmap(df.1.sub,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=F,show_colnames=F)
pheatmap(df.1.sub,cluster_rows = T,cluster_cols = F,annotation_names_row = T,show_rownames=F,show_colnames=T)
#
# pick 200 peaks in each module to show
df.x.sub <- c()  
for(i in 1:length(unique(dd$cluster))){
  df.x <- df[rownames(subset(dd,dd$cluster== i)),]
  df.x.sub <- c(df.x.sub,rownames(df.x)[sample(dim(df.x)[1],200)])
}
df.x.sub.m <- df[df.x.sub,]
dim(df.x.sub.m)
pheatmap(df.x.sub.m,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=F,show_colnames=F)
pheatmap(df.x.sub.m,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=F,show_colnames=F)
#
# calculate the average chromatin accessibility 
df.x.m <- c()
df.m <- c()
for(i in 1:length(unique(dd$cluster))){
  df.x.m <- colMeans(df[rownames(subset(dd,dd$cluster== i)),])
  df.m <- rbind(df.m,df.x.m)
}
dim(df.m)
head(df.m)
#
hist(apply(df.m , 1, sd))
#
pheatmap(df.m,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=F,show_colnames=F)
pheatmap(df.m,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=T,show_colnames=T)
#
df.m.1 <- df.m
rownames(df.m.1) <- paste("M",seq(1:30),sep = "")
pheatmap(df.m.1,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=F,show_colnames=T)
#
pheatmap(df.m.1,cluster_rows = T,cluster_cols = T,annotation_names_row = F,show_rownames=T,show_colnames=T)
pheatmap(df.m.1,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=T,show_colnames=T)
#
library(ComplexHeatmap)
head(df.m.1)
Heatmap(df.m.1,row_dend_reorder = T)
pheatmap(cor(t(df.m.1)))
pheatmap(cor(df.m.1))
#
##
# merge similar module
head(df.m.1) 
head(dd)
module.merge1 <- read.csv("Time.module.total.csv")
module.merge1 <- as.data.frame(module.merge1)
head(module.merge1)
dd.mm1 <- dd
dd.mm1$peak <- rownames(dd.mm1)
dd.mm1 <- merge(dd.mm1,module.merge1,by="cluster")
head(dd.mm1)
rownames(dd.mm1) <- dd.mm1$peak
head(dd.mm1)
#
table(dd.mm1$cluster1)
##
# calculate the average chromatin accessibility  after merging similar module
dd.mm1.x.m <- c()
dd.mm1.m <- c()
for(i in 1:length(unique(dd.mm1$cluster1))){
  dd.mm1.x.m  <- colMeans(df[rownames(subset(dd.mm1,dd.mm1$cluster1== i)),])
  dd.mm1.m <- rbind(dd.mm1.m,dd.mm1.x.m )
}
dim(dd.mm1.m) 
head(dd.mm1.m)
pheatmap(dd.mm1.m,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=T,show_colnames=T)
pheatmap(dd.mm1.m,cluster_rows = T,cluster_cols = F,annotation_names_row = F,show_rownames=F,show_colnames=T)
pheatmap(dd.mm1.m,cluster_rows = T,cluster_cols = T,annotation_names_row = F,show_rownames=F,show_colnames=T)
#
dd.mm1.mx <- dd.mm1.m
rownames(dd.mm1.mx) <- paste("M",seq(1:17),sep = "")
pheatmap(dd.mm1.mx,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=T,show_colnames=T)
#
df.rank <- df
df.rank.x.sub <- c()  
df.rank.x <- c()
table(dd.mm1$cluster1)
for(i in 1:length(table(dd.mm1$cluster1))){
  df.rank.x <- df.rank[rownames(subset(dd.mm1,dd.mm1$cluster1== i)),]
  df.rank.x.sub <- c(df.rank.x.sub,rownames(df.rank.x)[sample(dim(df.rank.x)[1],40)])
}
df.rank.x.sub.m <- df.rank[df.rank.x.sub,]
dim(df.rank.x.sub.m)
pheatmap(df.rank.x.sub.m,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=F,show_colnames=F)
#
#
####  filter noise in each final module after merging module
df.rank <- df
head(df.rank)
colnames(df.rank)
#
i=4
df.rank.x <- df.rank[rownames(subset(dd.mm1,dd.mm1$cluster1== i)),]
dim(df.rank.x)
pheatmap(df.rank.x ,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=F,show_colnames=T)
#pheatmap(df.rank.x ,cluster_rows = T,cluster_cols = F,annotation_names_row = F,show_rownames=F,show_colnames=T)
df.rank.x.m <- colMeans(df.rank[rownames(subset(dd.mm1,dd.mm1$cluster1== i)),])
df.rank.x.m
cor.tmp <- c()
cor.1 <- c()
dim(df.rank.x)[1]
for(i in 1:dim(df.rank.x)[1]){
  cor.tmp <- cor(as.numeric(df.rank.x[i,]),df.rank.x.m)
  cor.1 <- c(cor.1,cor.tmp)
}
hist(cor.1)
mean(cor.1)
median(cor.1)
length(which(cor.1>0.8))
df.rank.x.filter <- df.rank.x[which(cor.1>0.8),] 
dim(df.rank.x.filter)
pheatmap(df.rank.x.filter[1:1000,] ,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=F,show_colnames=F)
#
#
# filter
module.rank.s <- 1:17
#
df.rank.x.filter.merge <- c()
df.rank.x.filter.merge.length <- c()
for (j in 1:17) {
  df.rank.x <- df.rank[rownames(subset(dd.mm1,dd.mm1$cluster1== module.rank.s[j])),]
  df.rank.x.m <- colMeans(df.rank.x)
  cor.tmp <- c()
  cor.1 <- c()
  for(i in 1:dim(df.rank.x)[1]){
    cor.tmp <- cor(as.numeric(df.rank.x[i,]),df.rank.x.m)
    cor.1 <- c(cor.1,cor.tmp)
  }
  df.rank.x.filter.merge.length <- c(df.rank.x.filter.merge.length,length(which(cor.1>0.8))) 
  df.rank.x.filter.merge  <- rbind(df.rank.x.filter.merge,df.rank.x[which(cor.1>0.8),]) 
}
##
#
dim(df.rank.x.filter.merge)
barplot(log2(df.rank.x.filter.merge.length))
min(df.rank.x.filter.merge.length)
head(df.rank.x.filter.merge)
#
#
df.rank.x.filter.merge.length
length(df.rank.x.filter.merge.length)
table(df.rank.x.filter.merge.length)
barplot(log2(df.rank.x.filter.merge.length))
df.rank.x.filter.merge.cluster <- df.rank.x.filter.merge
head(df.rank.x.filter.merge)
#df.rank.x.filter.merge.cluster$cluster <- rep(seq(1:55),df.rank.x.filter.merge.length)
df.rank.x.filter.merge.cluster <- as.data.frame(df.rank.x.filter.merge.cluster)
df.rank.x.filter.merge.cluster$cluster <- rep(module.rank.s,df.rank.x.filter.merge.length)
table(df.rank.x.filter.merge.cluster$cluster )
#
head(df.rank.x.filter.merge.cluster)
dim(df.rank.x.filter.merge.cluster)
#
df.rank.x.filter.merge.cluster$cluster <- rep(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17),
                                              as.numeric(table(df.rank.x.filter.merge.cluster$cluster)))
df.rank.x.filter.merge.cluster <- df.rank.x.filter.merge.cluster[60616:163012,]
df.rank.x.filter.merge.cluster$cluster <- rep(c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16),
                                              as.numeric(table(df.rank.x.filter.merge.cluster$cluster)))
df.rank.x.filter.merge <- df.rank.x.filter.merge[rownames(df.rank.x.filter.merge.cluster),]
##
write.table(df.rank.x.filter.merge.cluster,"/home/sunkeyong/ZEPA/Peakset/Bulk_kmeans/zebrafish.dynamic.csv",quote=F, sep = " ",row.names = T,col.names = T)
#
##
# calculating the average chromatin accessibility after filtering 
df.rank.x.filter.merge.cluster.x.m <- c()
df.rank.x.filter.merge.cluster.m <- c()
for(i in 1:length(unique(df.rank.x.filter.merge.cluster$cluster))){
  df.rank.x.filter.merge.cluster.x.m <- colMeans(df.rank.x.filter.merge[rownames(subset(df.rank.x.filter.merge.cluster,df.rank.x.filter.merge.cluster$cluster== module.rank.s[i])),])
  df.rank.x.filter.merge.cluster.m <- rbind(df.rank.x.filter.merge.cluster.m,df.rank.x.filter.merge.cluster.x.m)
}
dim(df.rank.x.filter.merge.cluster.m)
head(df.rank.x.filter.merge.cluster.m)
rownames(df.rank.x.filter.merge.cluster.m) <- paste("M",1:16,sep = "")
pheatmap(df.rank.x.filter.merge.cluster.m,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=T,show_colnames=T)
pheatmap(df.rank.x.filter.merge.cluster.m,cluster_rows = T,cluster_cols = F,annotation_names_row = F,show_rownames=F,show_colnames=T)
pheatmap(df.rank.x.filter.merge.cluster.m,cluster_rows = T,cluster_cols = T,annotation_names_row = F,show_rownames=F,show_colnames=T)
#
df.rank.x.filter.merge.cluster.m.1 <- df.rank.x.filter.merge.cluster.m
rownames(df.rank.x.filter.merge.cluster.m.1) <- 1:16
head(df.rank.x.filter.merge.cluster.m.1)
pheatmap(df.rank.x.filter.merge.cluster.m.1,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=T,show_colnames=T)
pheatmap(df.rank.x.filter.merge.cluster.m.1,cluster_rows = T,cluster_cols = F,annotation_names_row = F,show_rownames=T,show_colnames=T)
#
#
# peak heatmap 
# pick 500 picks in each final peak module to show
head(df.rank.x.filter.merge)
dim(df.rank.x.filter.merge)
head(df.rank.x.filter.merge.cluster)
table(df.rank.x.filter.merge.cluster$cluster)
min(table(df.rank.x.filter.merge.cluster$cluster))
#

#
module.rank.s <- 1:16
df.rank.x.sub <- c()  
df.rank.x <- c()
for(i in 1:length(module.rank.s)){
  df.rank.x <- df.rank.x.filter.merge[rownames(subset(df.rank.x.filter.merge.cluster,df.rank.x.filter.merge.cluster$cluster== module.rank.s[i])),]
  df.rank.x.sub <- c(df.rank.x.sub,rownames(df.rank.x)[sample(dim(df.rank.x)[1],500)])  ### barplot(log2(df.rank.x.filter.merge.length))
}
length(df.rank.x.sub)
df.rank.x.sub.m <- df.rank.x.filter.merge[df.rank.x.sub,]
dim(df.rank.x.sub.m)
pheatmap(df.rank.x.sub.m,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=F,show_colnames=T)
#
#
# pick one module to shown
df.rank.1 <- df.rank.x.filter.merge[rownames(subset(df.rank.x.filter.merge.cluster,df.rank.x.filter.merge.cluster$cluster== 9)),]
dim(df.rank.1)
df.rank.1.sub <- df.rank.1[sample(dim(df.rank.1)[1],10),]
#pheatmap(df.rank.1.sub ,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=F,show_colnames=F)
pheatmap(df.rank.1.sub ,cluster_rows = F,cluster_cols = F,annotation_names_row = F,show_rownames=F,show_colnames=T)
#
df.rank.x.sub.m[1:3,1:3]
df.rank.x.sub.m.for.complexheatmap <- df.rank.x.sub.m
dim(df.rank.x.sub.m.for.complexheatmap)
dim(df.rank.x.sub.m.for.complexheatmap) 
split.module <- paste("M",rep(module.rank.s,as.integer(rep(500,16))),sep = "")
split.module  <- factor(split.module,levels =paste("M",module.rank.s,sep = "") )
#
#
cluster.rank <- c("4hpf","5hpf","6hpf","7hpf","8hpf","9hpf","10hpf",
                  "11hpf","12hpf","14hpf","18hpf","20hpf","22hpf","24hpf","30hpf","34hpf","38hpf",
                  "42hpf","48hpf","72hpf")
ZEPA.timeColor <- c("#FAD6A4","#F7B264","#E1520D","#AD3417",
                    "#F7C7D2","#EA66A1","#E50C84",
                    "#D2E2F2","#227DBD","#004695",
                    "#CFE6C0","#A8D392","#0C9E4A","#006D32",
                    "#C6C4E0","#857BB8","#694799",
                    "#EFBF00","#B78510","#935913")
#####
color.use.1 <-  ZEPA.timeColor
color.use.1 
sample_info.z = data.frame(CellType = cluster.rank )
col.use.1.tmp <- color.use.1
names(col.use.1.tmp) <- cluster.rank
col.use.1.tmp
class(col.use.1.tmp)
col.use.1.tmp2 <- list("CellType"=col.use.1.tmp)
col.use.1.tmp2
top_color <- HeatmapAnnotation(df=sample_info.z,col=col.use.1.tmp2, annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"))
##
peak.plot1 <- Heatmap(df.rank.x.sub.m.for.complexheatmap,
                      cluster_rows = F,
                      cluster_columns = F,
                      row_split = split.module,
                      show_column_names = T,
                      show_row_names = F,
                      gap = unit(0.5, "mm"), 
                      border = "black",
                      top_annotation=top_color, 
                      col = colorRampPalette(c("#439796","white", "#E93323"))(200))
draw(peak.plot1, heatmap_legend_side = "right", show_annotation_legend = T)
# 
df.rank.x.filter.merge.cluster$cluster <- as.numeric(df.rank.x.filter.merge.cluster$cluster)
table(df.rank.x.filter.merge.cluster$cluster )
# save the peaks information for homer analysis
for(i in 1:length(unique(df.rank.x.filter.merge.cluster$cluster))){
  x.p1 <- unlist(strsplit(rownames(subset(df.rank.x.filter.merge.cluster,df.rank.x.filter.merge.cluster$cluster== i)),"-"))
  x.p2 <- as.data.frame(cbind(x.p1[seq(1,length(x.p1),3)],x.p1[seq(2,length(x.p1),3)],x.p1[seq(3,length(x.p1),3)]))
  write.table(x.p2,file = paste("/home/sunkeyong/ZEPA/Peakset/Bulk_kmeans/cCREs/","M",i,".csv",sep = ""),quote=F, sep = "\t",row.names = F,col.names = F)
}
#
# analyse the bulk RNA-seq data from elife
#
TSSname <- as.data.frame(GRCz11.geneAnnotation$TSS)
TSSnamey <- as.data.frame(GRCz11.geneAnnotation$genes)
head(TSSname)
head(TSSnamey)
TSSnamey1 <- TSSnamey[which(TSSnamey$strand=="-"),]
TSSnamey2 <- TSSnamey[which(TSSnamey$strand=="+"),]
TSSnamey2$start <- TSSnamey2$end
#
TSSnamey3 <- rbind(TSSnamey1,TSSnamey2)
TSSnamey3$start <- as.numeric(TSSnamey3$start)
#
df.rank.x.filter.merge.cluster$cluster <- as.numeric(df.rank.x.filter.merge.cluster$cluster)
table(df.rank.x.filter.merge.cluster$cluster )
head(df.rank.x.filter.merge.cluster)
df.rank.x.filter.merge.cluster$chr <- unlist(strsplit(rownames(df.rank.x.filter.merge.cluster),"-"))[seq(1,3*nrow(df.rank.x.filter.merge.cluster),3)]
df.rank.x.filter.merge.cluster$start <- unlist(strsplit(rownames(df.rank.x.filter.merge.cluster),"-"))[seq(2,3*nrow(df.rank.x.filter.merge.cluster),3)]
df.rank.x.filter.merge.cluster$end <- unlist(strsplit(rownames(df.rank.x.filter.merge.cluster),"-"))[seq(3,3*nrow(df.rank.x.filter.merge.cluster),3)]
df.rank.x.filter.merge.cluster$center <- (as.numeric(df.rank.x.filter.merge.cluster$start)+as.numeric(df.rank.x.filter.merge.cluster$end))/2
#
table(df.rank.x.filter.merge.cluster$cluster)
# find the nearest genes of each peak in each module
x.genes <- as.list(c())
for (j in 1:16) {
  x1 <- df.rank.x.filter.merge.cluster[which(df.rank.x.filter.merge.cluster$cluster==j),]
  x1.genes <- c()
  x1.genes.tmp <- c()
  for (i in 1:nrow(x1)) {
    TSSnamey3.tmp <- TSSnamey3[which(TSSnamey3$seqnames==x1[i,]$chr),]
    x1.genes.tmp <- TSSnamey3.tmp[which.min(abs(TSSnamey3.tmp$start-x1[i,]$center)),]$symbol
    x1.genes <- c(x1.genes,x1.genes.tmp)
  }
  x.genes[[j]] <- as.data.frame(table(x1.genes))
}
#
x1.genes <- as.data.frame(table(x1.genes))
head(x1.genes)
hist(x1.genes$Freq)
x1.genes[order(x1.genes$Freq,decreasing = T),]
# filter genes with frequency lower than 2
x1.genes[which(x1.genes$Freq> 2),]$x1.genes 
#
load("/home/sunkeyong/ZEPA/E-ERAD-475-atlasExperimentSummary.Rdata")
experimentSummary$rnaseq@assays$data$counts
dim(experimentSummary$rnaseq@assays$data$counts)
xxx.anno <- read.csv("/home/sunkeyong/ZEPA/GRCz11_anno/Lawson/v4_3_2geneinfo3.csv")
dim(xxx.anno)
xxx.anno$Ens99geneIDversion2 <- unlist(strsplit(xxx.anno$Ens99geneIDversion,".",fixed=T))[seq(1,2*nrow(xxx.anno),2)]
length(unique(xxx.anno$Ens99geneIDversion2))
xxx.anno <- xxx.anno[!duplicated(xxx.anno[,"Ens99geneIDversion2"]),]
rownames(xxx.anno) <- xxx.anno$Ens99geneIDversion2
#
RNAseqcount <- experimentSummary$rnaseq@assays$data$counts
dim(RNAseqcount)
head(RNAseqcount)
RNAseqcount <- RNAseqcount[xxx.anno$Ens99geneIDversion2,]
rownames(RNAseqcount) <- xxx.anno$LLgeneSymbol
RNAseqcount.meta <- read.table("/home/sunkeyong/ZEPA/E-ERAD-475-experiment-design.tsv",header = T)
length(unique(RNAseqcount.meta$Run))
length(RNAseqcount.meta$Run)
rownames(RNAseqcount.meta) <- RNAseqcount.meta$Run
RNAseqcount.meta <- RNAseqcount.meta[colnames(RNAseqcount),]
#
ZBulkRNAseq <- CreateSeuratObject(counts = RNAseqcount, project = "zBulkRNAseq", min.cells = 0, min.features = 0)
ZBulkRNAseq
ZBulkRNAseq$stage <- RNAseqcount.meta$stage
ZBulkRNAseq@active.ident <- as.factor(ZBulkRNAseq$stage)
levels(ZBulkRNAseq) <- c("zygote","cleavage2-cell","blastula128-cell","blastula1k-cell","blastuladome","gastrula50%-epiboly",
                  "gastrulashield","gastrula75%-epiboly","segmentation1-4somites","segmentation14-19somites","segmentation20-25somites",
                  "pharyngulaprim-5","pharyngulaprim-15","pharyngulaprim-25","hatchinglong-pec",
                  "larvalprotrudingmouth","larvalday4","larvalday5")
ZBulkRNAseq <- subset(ZBulkRNAseq,idents=c("blastuladome","gastrula50%-epiboly",
                             "gastrulashield","gastrula75%-epiboly","segmentation1-4somites","segmentation14-19somites","segmentation20-25somites",
                             "pharyngulaprim-5","pharyngulaprim-15","pharyngulaprim-25","hatchinglong-pec",
                             "larvalprotrudingmouth"))
#
ZBulkRNAseq <- NormalizeData(ZBulkRNAseq, normalization.method = "LogNormalize", scale.factor = 10000)
VlnPlot(ZBulkRNAseq,features = "myod1",group.by = "stage")+NoLegend()
all.genes <- rownames(ZBulkRNAseq)
ZBulkRNAseq <- ScaleData(ZBulkRNAseq, features = all.genes)
#
# Create a list of a vector of gene names for score computing
# genes.for.scoring <- list(intersect(as.character(x1.genes[which(x1.genes$Freq> 3),]$x1.genes),rownames(ZBulkRNAseq)))
#
# ZBulkRNAseq <- AddModuleScore(object = ZBulkRNAseq, features = genes.for.scoring, name = "test", random.seed = 1)
# VlnPlot(ZBulkRNAseq,features = "test1")+NoLegend()
#
for (i in 1:16) {
  genes.for.scoring <- list(intersect(as.character(x.genes[[i]][which(x.genes[[i]]$Freq> 0),]$x1.genes),rownames(ZBulkRNAseq)))
  ZBulkRNAseq <- AddModuleScore(object = ZBulkRNAseq, features = genes.for.scoring, name = paste(i,"C",sep = ""), random.seed = 1)
}
VlnPlot(ZBulkRNAseq,features = "X1C1")+NoLegend()
VlnPlot(ZBulkRNAseq,features = "X16C1")+NoLegend()
#
ZBulkRNAseq.meta2 <- ZBulkRNAseq@meta.data
dim(ZBulkRNAseq.meta2)
#
ZBulkRNAseq2 <- CreateSeuratObject(counts = t(ZBulkRNAseq.meta2[,5:20]), project = "zBulkRNAseq", min.cells = 0, min.features = 0)
ZBulkRNAseq2
ZBulkRNAseq2$stage <- ZBulkRNAseq$stage
ZBulkRNAseq2@active.ident <- as.factor(ZBulkRNAseq2$stage)
levels(ZBulkRNAseq2) <- c("blastuladome","gastrula50%-epiboly",
                   "gastrulashield","gastrula75%-epiboly","segmentation1-4somites","segmentation14-19somites","segmentation20-25somites",
                   "pharyngulaprim-5","pharyngulaprim-15","pharyngulaprim-25","hatchinglong-pec",
                   "larvalprotrudingmouth")
#
ZBulkRNAseq2.avg <- as.data.frame(AverageExpression(ZBulkRNAseq2,slot="counts"))
pheatmap(ZBulkRNAseq2.avg,cluster_rows = F,cluster_cols = F)
pheatmap(ZBulkRNAseq2.avg,cluster_rows = F,cluster_cols = F,scale = "row",color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100))
#
normalizex <- function(x){
  return ((x-min(x))/(max(x)-min(x)))
}
pheatmap(t(apply(ZBulkRNAseq2.avg, 1, normalizex)),cluster_rows = F,cluster_cols = F,color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100))
#
#
