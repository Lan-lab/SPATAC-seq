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
library(purrr)
library(chromVARmotifs)
library('ggseqlogo')
#
addArchRThreads(threads = 30) 
#
addArchRGenome("mm10")
#
# Based on the results of Figure 1 and Figure 2
load("MOPA.RData")
MOPA <- addMotifAnnotations(ArchRProj = MOPA, motifSet = "cisbp", name = "Motif",force=T)
MOPA <- addBgdPeaks(MOPA,force=T)
MOPA <- addDeviationsMatrix(
  ArchRProj = MOPA, 
  peakAnnotation = "Motif",
  force = TRUE
)
#
# plot some key TF motif
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "GeneScoreMatrix", name = c("Eomes"), quantCut = c(0.25, 0.99),
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("solarExtra"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/TFfeatureplot/Eomes.genescore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "MotifMatrix", name = c("z:Eomes_768"), quantCut = c(0.1, 0.99),
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("blueYellow"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/TFfeatureplot/Eomes.deviationscore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "GeneScoreMatrix", name = c("Mef2c"), 
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("solarExtra"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/TFfeatureplot/Mef2c.genescore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "MotifMatrix", name = c("z:Mef2c_638"), 
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("blueYellow"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/TFfeatureplot/Mef2c.deviationscore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "MotifMatrix", name = c("z:Sox2_750"), 
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("blueYellow"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/TFfeatureplot/Sox2.deviationscore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "MotifMatrix", name = c("z:Pou5f1_618"), 
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("blueYellow"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/TFfeatureplot/Pou5f1.deviation.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "MotifMatrix", name = c("z:Gata1_387"), 
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("blueYellow"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/TFfeatureplot/gata1.deviationscore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "MotifMatrix", name = c("z:Trp63_853"), 
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("blueYellow"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/TFfeatureplot/trp63.deviationscore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "MotifMatrix", name = c("z:Pax6_614"), 
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("blueYellow"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/TFfeatureplot/Pax6.deviationscore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "MotifMatrix", name = c("z:Pax8_854"), 
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("blueYellow"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/TFfeatureplot/Pax8.deviationscore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "MotifMatrix", name = c("z:Emx1_491"), 
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("blueYellow"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/TFfeatureplot/Emx1.deviationscore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "MotifMatrix", name = c("z:Erg_287"), 
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("blueYellow"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/TFfeatureplot/Erg.deviationscore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "MotifMatrix", name = c("z:Hnf4g_664"), 
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("blueYellow"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/TFfeatureplot/Hnf4g.deviationscore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "MotifMatrix", name = c("z:Hnf4g_664"), 
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("blueYellow"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/TFfeatureplot/Hnf4g.deviationscore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "MotifMatrix", name = c("z:Neurod2_69"), 
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("blueYellow"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/TFfeatureplot/Neurod2.deviationscore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "MotifMatrix", name = c("z:Tbr1_769"), 
                   quantCut = c(0.25, 0.99), 
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("blueYellow"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/TFfeatureplot/Tbr1.deviationscore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "MotifMatrix", name = c("z:Hnf4a_665"), 
                   quantCut = c(0.35, 0.99), 
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("blueYellow"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/TFfeatureplot/HNf4a.deviationscore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "MotifMatrix", name = c("z:Klf1_214"), 
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("blueYellow"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/motifseq/Klf1.deviationscore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "MotifMatrix", name = c("z:Klf2_819"), 
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("blueYellow"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/motifseq/Klf2.deviationscore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "MotifMatrix", name = c("z:Klf4_143"), 
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("blueYellow"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/motifseq/Klf4.deviationscore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "MotifMatrix", name = c("z:Klf3_804"), 
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("blueYellow"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/motifseq/Klf3.deviationscore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "MotifMatrix", name = c("z:Klf6_794"), 
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("blueYellow"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/motifseq/Klf6.deviationscore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "MotifMatrix", name = c("z:Klf10_810"), 
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("blueYellow"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/motifseq/Klf10.deviationscore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
p <- plotEmbedding(ArchRProj = MOPA, colorBy = "MotifMatrix", name = c("z:Klf11_798"), 
                   embedding = "PeakMatrix_IterativeLSI_tSNE",imputeWeights = getImputeWeights(MOPA),
                   pal = paletteContinuous("blueYellow"))+
  NoLegend()+ theme(plot.title = element_blank())+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p 
ggsave("/motifseq/Klf11.deviationscore.png",plot=p,device="png",dpi=300,units = "cm",width = 7,height = 7)
#
##
m <- read.table("motifseq/klf1.csv",sep = ",")
pwm <- makePWM(t(m)/colSums(as.matrix(t(m))))
seqLogo(pwm)
#
m <- read.table("motifseq/csv/klf1.csv",sep = ",")
klf1.m <- t(m)/colSums(as.matrix(t(m)))
rownames(klf1.m) <- c("A","C","G","T")
klf1.motif <-new("pcm", mat=as.matrix(klf1.m), name="klf1")
plot(klf1.motif)
klf <- list(klf1.motif,klf1.motif)
klf
motifStack(klf, layout="stack", ncex=1.0)
#
m1 <- read.table("motifseq/csv/klf1.csv",sep = ",");klf1.m <- t(m1);rownames(klf1.m) <- c("A","C","G","T")
m2 <- read.table("motifseq/csv/klf2.csv",sep = ",");klf2.m <- t(m2);rownames(klf2.m) <- c("A","C","G","T")
m3 <- read.table("motifseq/csv/klf3.csv",sep = ",");klf3.m <- t(m3);rownames(klf3.m) <- c("A","C","G","T")
m4 <- read.table("motifseq/csv/klf4.csv",sep = ",");klf4.m <- t(m4);rownames(klf4.m) <- c("A","C","G","T")
m5 <- read.table("motifseq/csv/klf5.csv",sep = ",");klf5.m <- t(m5);rownames(klf5.m) <- c("A","C","G","T")
m6 <- read.table("motifseq/csv/klf6.csv",sep = ",");klf6.m <- t(m6);rownames(klf6.m) <- c("A","C","G","T")
m9 <- read.table("motifseq/csv/klf9.csv",sep = ",");klf9.m <- t(m9);rownames(klf9.m) <- c("A","C","G","T")
m10 <- read.table("motifseq/csv/klf10.csv",sep = ",");klf10.m <- t(m10);rownames(klf10.m) <- c("A","C","G","T")
m11 <- read.table("motifseq/csv/klf11.csv",sep = ",");klf11.m <- t(m11);rownames(klf11.m) <- c("A","C","G","T")
m12 <- read.table("motifseq/csv/klf12.csv",sep = ",");klf12.m <- t(m12);rownames(klf12.m) <- c("A","C","G","T")
m13 <- read.table("motifseq/csv/klf13.csv",sep = ",");klf13.m <- t(m13);rownames(klf13.m) <- c("A","C","G","T")
m14 <- read.table("motifseq/csv/klf14.csv",sep = ",");klf14.m <- t(m14);rownames(klf14.m) <- c("A","C","G","T")
m15 <- read.table("motifseq/csv/klf15.csv",sep = ",");klf15.m <- t(m15);rownames(klf15.m) <- c("A","C","G","T")
m16 <- read.table("motifseq/csv/klf16.csv",sep = ",");klf16.m <- t(m16);rownames(klf16.m) <- c("A","C","G","T")
m17 <- read.table("motifseq/csv/klf17.csv",sep = ",");klf17.m <- t(m17);rownames(klf17.m) <- c("A","C","G","T")
#
klf1.motif <-new("pcm", mat=as.matrix(klf1.m), name="klf1");plot(klf1.motif)
klf2.motif <-new("pcm", mat=as.matrix(klf2.m), name="klf2");plot(klf2.motif)
klf3.motif <-new("pcm", mat=as.matrix(klf3.m), name="klf3");plot(klf3.motif)
klf4.motif <-new("pcm", mat=as.matrix(klf4.m), name="klf4");plot(klf4.motif)
klf5.motif <-new("pcm", mat=as.matrix(klf5.m), name="klf5");plot(klf5.motif)
klf6.motif <-new("pcm", mat=as.matrix(klf6.m), name="klf6");plot(klf6.motif)
klf9.motif <-new("pcm", mat=as.matrix(klf9.m), name="klf9");plot(klf9.motif)
klf10.motif <-new("pcm", mat=as.matrix(klf10.m), name="klf10");plot(klf10.motif)
klf11.motif <-new("pcm", mat=as.matrix(klf11.m), name="klf11");plot(klf11.motif)
klf12.motif <-new("pcm", mat=as.matrix(klf12.m), name="klf12");plot(klf12.motif)
klf13.motif <-new("pcm", mat=as.matrix(klf13.m), name="klf13");plot(klf13.motif)
klf14.motif <-new("pcm", mat=as.matrix(klf14.m), name="klf14");plot(klf14.motif)
klf15.motif <-new("pcm", mat=as.matrix(klf15.m), name="klf15");plot(klf15.motif)
klf16.motif <-new("pcm", mat=as.matrix(klf16.m), name="klf16");plot(klf16.motif)
klf17.motif <-new("pcm", mat=as.matrix(klf17.m), name="klf17");plot(klf17.motif)
#
#
#
klf <- list(klf1.motif,klf2.motif,klf3.motif,klf4.motif,
            klf5.motif,klf6.motif,klf9.motif,klf10.motif,
            klf11.motif,klf12.motif,klf13.motif,klf14.motif,
            klf15.motif,klf16.motif,klf17.motif)
klf
motifStack(klf, layout="stack", ncex=1.0)
#
VlnPlot(MOCA,features = c("Emx2") ,pt.size = 0,ncol=1)+NoLegend()
StackedVlnPlot(MOCA,features = c("Emx2"),pt.size = 0,ncol=1)+NoLegend()
StackedVlnPlot(MOCA,features = c("Klf12","Klf13","Klf14","Klf15","Klf16"),pt.size = 0,cols =color.use.em)+NoLegend()
StackedVlnPlot(MOCA,features = c("Klf1","Klf2","Klf3","Klf4","Klf5","Klf6","Klf9","Klf10","Klf11"),pt.size = 0,cols =color.use.em)+NoLegend()
#
#
# import the motif deviation matrix
Get.motif.mat <- getMatrixFromProject(ArchRProj = MOPA,useMatrix = "MotifMatrix")
rownames(Get.motif.mat@assays@data$deviations) <- Get.motif.mat@elementMetadata$name
Get.deviations.matrix <- Get.motif.mat@assays@data$deviations
head(Get.deviations.matrix)
dim(Get.deviations.matrix)
rownames(Get.motif.mat@assays@data$z) <- Get.motif.mat@elementMetadata$name
Get.z.matrix <- Get.motif.mat@assays@data$z
head(Get.z.matrix)
dim(Get.z.matrix)
#
save(Get.z.matrix,file = "matrix/Get.z.matrix.RData")
#
plotVarDev <- getVarDeviations(MOPA.filter, name = "MotifMatrix", plot = TRUE)
plotVarDev
plotVarDev$data$name[1:20]
#
Get.deviations <- CreateSeuratObject(counts = Get.z.matrix, assay = "deviations", project = "z.deviations", min.cells = 0, min.features = 0)
Get.deviations
##
MOPA.metadata <- as.data.frame(MOPA@cellColData)
MOPA.metadata <- MOPA.metadata[rownames(Get.deviations@meta.data),]
Get.deviations$ID <- MOPA.metadata$ID
#
Get.deviations@active.ident <- as.factor(Get.deviations$ID)
downsample.num <- 300
Get.deviations.sub <- subset(Get.deviations, downsample = downsample.num)
Get.deviations.sub
Get.deviations.sub <- ScaleData(Get.deviations.sub, features = rownames(Get.deviations.sub))
#
# load the PseudoCell function
#
#
PseudoCell.num <- 60
MOPA.motif <- PseudoCell(Get.deviations.sub, "deviations","data","ID",PseudoCell.num)
MOPA.motif
rownames(MOPA.motif@meta.data)
MOPA.motif$ID <- unlist(strsplit(rownames(MOPA.motif@meta.data),"_Cell"))[seq(1,2*length(rownames(MOPA.motif@meta.data)),2)]
MOPA.motif@active.ident <- as.factor(MOPA.motif$ID)
levels(MOPA.motif) <- paste("C",1:33,sep="")
MOPA.motif <- ScaleData(MOPA.motif, features = rownames(MOPA.motif))
#
DoHeatmap(MOPA.motif, features = "Zfp263-877") + NoLegend()
#
#
x.dev <- MOPA.motif@assays$RNA@scale.data
rownames(x.dev) <- paste(unlist(strsplit(rownames(x.dev),"-"))[seq(1,2*dim(x.dev)[1],2)],
                         unlist(strsplit(rownames(x.dev),"-"))[seq(2,2*dim(x.dev)[1],2)],sep = "_")

dim(x.dev)
x.dev[1:3,1:3]

rep(c("Cell0","Cell1","Cell2","Cell3","Cell4"),33)
rep(paste("C",1:33,sep=""),each=5)
cellrank <- paste(rep(paste("C",1:33,sep=""),each=5),
                  rep(c("Cell0","Cell1","Cell2","Cell3","Cell4"),33),
                  sep = "_")
x.dev <- x.dev[,cellrank]
x.dev[1:3,1:3]
#
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
#
col_fun = colorRamp2(seq(-127,128,1)/64, paletteContinuous("blueYellow"))
#
pheatmap(x.dev)
pheatmap(x.dev2)
##
Heatmap(x.dev2,
        cluster_rows = T,
        col = col_fun,
        cluster_columns = T,
        show_column_names = F,
        show_row_names = F)
#
sample_info = data.frame(CellType = rep(as.character(as.data.frame(table(MOPA.motif@active.ident))$Var1),each=downsample.num/PseudoCell.num))
#sample_info = rep(paste("C",1:33,sep=""),each=5)
#
col.use.1.tmp <- color.use.em
names(col.use.1.tmp) <- as.data.frame(table(MOPA.motif@active.ident))$Var1
col.use.1.tmp
class(col.use.1.tmp)
col.use.1.tmp2 <- list("CellType"=col.use.1.tmp)
col.use.1.tmp2
top_color <- HeatmapAnnotation(df=sample_info,col=col.use.1.tmp2, annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"))
#
motif.plot1 <- Heatmap(x.dev,
                       row_km=0,
                       col = col_fun,
                       cluster_rows = T,
                       cluster_columns = T,
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
                       top_annotation=top_color,
                       show_column_dend=F,
                       show_row_dend=F, row_dend_reorder = TRUE)
#
draw(motif.plot1, heatmap_legend_side = "right", show_annotation_legend = T)
x.dev.cor <- cor(x.dev)
Heatmap(x.dev.cor,
        cluster_rows = T,
        cluster_columns = T,
        top_annotation=top_color,
        col = colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100),
        show_column_dend=F,
        show_row_dend=F,
        show_column_names = F,
        use_raster=T,
        show_row_names = F)
#
dim(x.dev)
x.dev[1:3,1:3]
x.dev[x.dev < -2]<- -2
x.dev[x.dev > 2]<- 2
#
length(which (plotVarDev$data$combinedVars > 1.6)) ## 411
length(which (plotVarDev$data$combinedVars > 1.5)) ## 450
length(which (plotVarDev$data$combinedVars > 2)) ## 303
# 
x.dev.var <- x.dev[plotVarDev$data$idx[1:450],] 
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
                       use_raster=T,
                       show_row_dend=F, row_dend_reorder = TRUE)
#
draw(motif.plot2, heatmap_legend_side = "right", show_annotation_legend = T)
##
x.dev.var2 <- x.dev.var[rownames(x.dev.var)[row_order(motif.plot2)],]
rownames(x.dev.var2)
marker_TFs <- c("Lef1-734")
marker_TFs <- as.character(plotVarDev$data$name[1:40])
#marker_TFs <- subset(cor.metadata.merge,maxDelta > maxDelta_cutoff & cor > cor_cutoff)$MotifMatrix_name
marker_TFs <- c("Rfx2_707","Rfx5_704","Sox2_750","Sox12_741","Hoxb4_511","Hoxa5_510","Neurog1_85","Neurod1_784","NP_602","NP_596",
                "Hoxa10_395","Hoxb13_547","Foxl1_361","Foxl1_360","Mef2c_638","Mef2b_882","Gata2_383","Gata3_384",
                "Nr2f1_690","Nr1h2_688","Klf16_874","Klf1_214")
marker_TFs <- c("Zic5_195","Rfx2_707","Pou5f1_618","Sox12_741","Dlx4_427","Hoxa2_420","Esx1_444","Emx2_535","Nfib_859","Twist2_21",
                "Neurog1_85","Eomes_768","Myog_45","Ebf1_90","Foxd4_833","Mef2c_638","Gata6_382","Hnf4g_664","Otx2_437","Pbx3_513","Klf4_143")
marker_TFs <- intersect(rownames(x.dev.var2),unique(marker_TFs))
marker_TFs
TF_pos <- which(rownames(x.dev.var2) %in% marker_TFs)
row_anno <-  rowAnnotation(marker_TFs = anno_mark(at = TF_pos, 
                                                  labels = marker_TFs))
motif.plot3 <- Heatmap(x.dev.var2,
                       row_km=0,
                       #col = colorRamp2(seq(-127,128,1)/64, paletteContinuous("horizonExtra")),
                       #col = colorRamp2(seq(-127,128,1)/64, paletteContinuous("solarExtra")),
                       col = colorRamp2(seq(-127,128,1)/64, paletteContinuous("blueYellow")),
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
                       use_raster=T,
                       right_annotation = row_anno)
#
draw(motif.plot3, heatmap_legend_side = "right", show_annotation_legend = F)
##
TF.names <- unlist(strsplit(rownames(x.dev.var2),"_"))[seq(1,2*450,2)]
TF.allnames <- unlist(strsplit(rownames(mypbmc),"-"))[seq(1,2*884,2)]
length(unique(TF.allnames))
length(TF.allnames)
#
TF.names.inter <- intersect(rownames(pbmc.Main_cell_type.Stage.pseudo.seurat),TF.names)
TF.allnames.inter <- intersect(rownames(pbmc.Main_cell_type.Stage.pseudo.seurat),TF.allnames)
length(unique(TF.allnames.inter))
length(TF.allnames.inter)
length(unique(TF.names.inter))
length(TF.names.inter)
#
x.dev.var.x <- x.dev.var2
rownames(x.dev.var.x) <- unlist(strsplit(rownames(x.dev.var.x),"_"))[seq(1,2*450,2)]
x.dev.var.x <- x.dev.var.x[TF.names.inter,]
dim(x.dev.var.x)
motif.plot4 <- Heatmap(x.dev.var.x,
                       row_km=0,
                       #col = colorRamp2(seq(-127,128,1)/64, paletteContinuous("horizonExtra")),
                       #col = colorRamp2(seq(-127,128,1)/64, paletteContinuous("solarExtra")),
                       col = colorRamp2(seq(-127,128,1)/64, paletteContinuous("blueYellow")),
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
                       use_raster=T,
                       right_annotation = row_anno)
#
draw(motif.plot4, heatmap_legend_side = "right", show_annotation_legend = F)
##
dim(x.dev)
x.dev[1:3,1:3]
x.dev.cor <- cor(x.dev)
Heatmap(x.dev.cor,
        cluster_rows = T,
        cluster_columns = T,
        top_annotation=top_color,
        col = colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100),
        show_column_dend=F,
        show_row_dend=F,
        show_column_names = F,
        use_raster=T,
        show_row_names = F)
write.csv(x.dev.cor,file = "Source_Data/Figure3A.csv")
#
#
x.dev.cor.t <- cor(t(x.dev))
Heatmap(x.dev.cor.t,
        show_column_dend=T, ## 
        show_row_dend=T, ## 
        show_column_names = F,
        use_raster=T,
        show_row_names = F,
        #row_km=15,
        #column_km=15,
        col = colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100))
write.csv(x.dev.cor.t,file = "Source_Data/FigureS3A.csv")
write.csv(x.dev.var2,file = "Source_Data/Figure3B.csv")
#
##
#0620 #
motifPositions <- getPositions(MOPA)
motifs <- c("Myo")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs
#
seFoot <- getFootprints(
  ArchRProj = MOPA, 
  positions = motifPositions[plotVarDev$data$name[751:884]], 
  groupBy = "cnum"
)
plotFootprints(
  seFoot = seFoot,
  ArchRProj = MOPA, 
  normMethod = "Subtract",
  plotName = "plotVarDev-footprinting751-884",
  addDOC = FALSE,
  smoothWindow = 5
)
#
seFoot <- getFootprints(
  ArchRProj = MOPA, 
  positions = motifPositions[c("Mef2a_640","Sfpi1_265","Pou4f1_841","Pou4f2_840","Myf6_63","Grhl1_390","Klf1_214")], 
  groupBy = "cnum"
)
color.use.em <- c("#FDCDAC","#FFFF99","#B3CDE3","#6495ED","#00FFFF", ## 1--5 
                  "#FBB4AE","#F4CAE4","#ABDDA4","#B15928","#D6604D", ## 6--10
                  "#999999","#00BFFF","#CC99FF","#48D1CC","#FFED6F", ## 11--15 
                  "#9E0142","#0099FF","#B3E2CD","#276419","#66C2A5", ## 16--20
                  "#B3EE3A","#6A3D9A","#87CEFA","#ABD9E9","#B2ABD2", ## 21--25 
                  "#0000FF","#FF33CC","#66BD63","#E08214","#DE77AE", ## 26--30
                  "#9970AB","#EE82EE","#99CCCC") 
##
plotFootprints(
  seFoot = seFoot,
  pal = color.use.em,
  ArchRProj = MOPA, 
  normMethod = "Subtract",
  plotName = "plotVarDev-footprinting-0705.all",
  addDOC = FALSE,
  smoothWindow = 5
)
#
#
data("mouse_pwms_v1") 
data("mouse_pwms_v2") 
MouseNr5a2.V2 <- mouse_pwms_v2$ENSMUSG00000026398_LINE1611_Nr5a2_D_N1
MouseNr5a2.V1_1 <- mouse_pwms_v1$ENSMUSG00000026398_LINE1611_Nr5a2_D_N1
MouseNr5a2.V1_2 <- mouse_pwms_v1$ENSMUSG00000026398_LINE1612_Nr5a2_D_N1
#
#
ggseqlogo(MouseNr5a2.V2@profileMatrix)
ggseqlogo(danRer11.Motif.raw2$`nr5a1a+nr5a2+nr5a5`@profileMatrix)
#
#
#
#
zA <- c(10,2,1,193,194,1,0,25,1,137)
zC <- c(54,74,185,1,0,0,1,64,181,8)
zG <- c(56,1,8,0,1,193,193,12,1,21)
zT <- c(75,118,1,1,0,1,1,94,12,29)
#
zNr5a2 <- as.matrix(rbind(zA,zT,zG,zC))
rownames(zNr5a2) <- c("A","T","G","C")
#ggseqlogo(as.matrix(pfm2ppm(zNr5a2)))
#
motif <- new("pcm", mat=as.matrix(zNr5a2), name="zNr5a2")
plot(motif)
#
### jaspar
mA <- c(645,985,276,122,71,8,1685,1679,39,0,169,0,1437,271,249)
mC <- c(224,170,198,378,351,1616,14,0,0,0,584,1568,20,345,693)
mG <- c(634,336,990,366,23,74,0,23,1663,1702,61,12,90,740,454)
mT <- c(199,211,238,836,1257,4,3,0,0,0,888,122,155,346,306)
#
mNr5a2 <- as.matrix(rbind(mA,mT,mG,mC))
rownames(mNr5a2) <- c("A","T","G","C")
#ggseqlogo(as.matrix(pfm2ppm(mNr5a2)))
#
motif.m <- new("pcm", mat=as.matrix(mNr5a2), name="mNr5a2")
plot(motif.m)
#
#
mA <- c(0.184210526315789,0.823529411764706,0.949152542372881,0.0169491525423729,
        0.0169491525423729,0.0847457627118644,0.0178571428571429,0.7)
mC <- c(0.552631578947368,0.0196078431372549,0.0169491525423729,0.0169491525423729,
        0.0169491525423729,0.135593220338983,0.910714285714286,0.025)
mG <- c(0.184210526315789,0.137254901960784,0.0169491525423729,0.949152542372881,
        0.949152542372881,0.0169491525423729,0.0178571428571429,0.25)
mT <- c(0.0789473684210526,0.0196078431372549,0.0169491525423729,0.0169491525423729,
        0.0169491525423729,0.76271186440678,0.0535714285714286,0.025)
#
mNr5a2 <- as.matrix(rbind(mA,mT,mG,mC))
rownames(mNr5a2) <- c("A","T","G","C")
#ggseqlogo(as.matrix(pfm2ppm(mNr5a2)))
#
motif.m <- new("pcm", mat=as.matrix(mNr5a2), name="mNr5a2")
plot(motif.m)
#
save.image("MOPA.Figure3.for_motif_analysis.RData")
#
