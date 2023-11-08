setwd("/home/sunkeyong/MOPA_project/Label_trasfer")
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
addArchRThreads(threads = 13) 
#
# processing scRNA-seq data
### Download TOME MOCA scRNA-seq and load
NE10.5 <- readRDS("/home/sunkeyong/MOPA_project/tome/TOME_raw_RDS/seurat_object_E10.5.rds") 
NE11.5 <- readRDS("/home/sunkeyong/MOPA_project/tome/TOME_raw_RDS/seurat_object_E11.5.rds") 
NE12.5 <- readRDS("/home/sunkeyong/MOPA_project/tome/TOME_raw_RDS/seurat_object_E12.5.rds") 
NE13.5 <- readRDS("/home/sunkeyong/MOPA_project/tome/TOME_raw_RDS/seurat_object_E13.5.rds") 
#
TOME_M <- merge(NE10.5,y=c(NE11.5,NE12.5,NE13.5))
rm(NE10.5)
rm(NE11.5)
rm(NE12.5)
rm(NE13.5)
#
gc()
#
TOME_M.metadata <- TOME_M@meta.data
head(TOME_M.metadata)
#
save(TOME_M.metadata,file = "/home/sunkeyong/MOPA_project/tome/TOME_process_data/TOME_M.metadata")
#
#
Nature_annotation <- read.csv("/home/sunkeyong/MOPA_project/tome/sci_RNA_seq3/cell_annotate.csv",header = T,sep = ",")
head(Nature_annotation)
Nature_annotation <- as.data.frame(Nature_annotation)
rownames(Nature_annotation) <- as.character(Nature_annotation$sample)
Nature_annotation$Barcode <- as.character(Nature_annotation$sample)
#
Nature_gene_annotation <- read.csv("/home/sunkeyong/MOPA_project/tome/sci_RNA_seq3/gene_annotate.csv",header = T,sep = ",")
head(Nature_gene_annotation)
#
table(Nature_gene_annotation$gene_short_name)
length(Nature_gene_annotation$gene_short_name)
length(unique(Nature_gene_annotation$gene_short_name))
#
head(Nature_gene_annotation)
Nature_gene_annotation$gene_id <- as.character(Nature_gene_annotation$gene_id)
Nature_gene_annotation$gene_id2 <- unlist(strsplit(Nature_gene_annotation$gene_id,".",fixed = T))[seq(1,2*dim(Nature_gene_annotation)[1],2)]
length(unique(Nature_gene_annotation$gene_id2))
rownames(Nature_gene_annotation) <- Nature_gene_annotation$gene_id2
head(Nature_gene_annotation)
#
#
#
TOME_M.GeneName <- rownames(TOME_M@assays$RNA@data)
head(TOME_M.GeneName)
length(TOME_M.GeneName)
#
length(intersect(Nature_gene_annotation$gene_id2,TOME_M.GeneName))
length(unique(intersect(Nature_gene_annotation$gene_id2,TOME_M.GeneName)))
length(setdiff(Nature_gene_annotation$gene_id2,TOME_M.GeneName))
length(setdiff(TOME_M.GeneName,Nature_gene_annotation$gene_id2))
#
intersect(Mouse_feature,Nature_gene_annotation[intersect(Nature_gene_annotation$gene_id2,TOME_M.GeneName),]$gene_short_name)
intersect(Mouse_feature,Nature_gene_annotation$gene_short_name)
#
TOME_Mmatrix <- TOME_M@assays$RNA@counts
TOME_Mmatrix
dim(TOME_Mmatrix)
#
TOME_Mmatrix2 <- TOME_Mmatrix[intersect(Nature_gene_annotation$gene_id2,TOME_M.GeneName),]
TOME_Mmatrix2
dim(TOME_Mmatrix2)
#
dim(Nature_gene_annotation)
Nature_gene_annotation.tmp <- Nature_gene_annotation[intersect(Nature_gene_annotation$gene_id2,TOME_M.GeneName),]
dim(Nature_gene_annotation.tmp)
#
rownames(TOME_Mmatrix2) <- Nature_gene_annotation.tmp$gene_short_name
TOME_Mmatrix2[1:3,1:3]
#
dim(TOME_Mmatrix2)
#
TOME_Mmatrix3 <- TOME_Mmatrix2[intersect(Mouse_feature,rownames(TOME_Mmatrix2)),]
TOME_Mmatrix3[1:3,1:3]
#
save(TOME_Mmatrix3,file = "/home/sunkeyong/MOPA_project/tome/TOME_process_data/TOME_Mmatrix3.Final.RData")
#
co_cells <- intersect(TOME_M.metadata$sample,Nature_annotation$Barcode)
#
TOME_Mmatrix4 <- TOME_Mmatrix3[,co_cells]
dim(TOME_Mmatrix4)
rm(TOME_Mmatrix3)
gc()
save(TOME_Mmatrix4,file = "/home/sunkeyong/MOPA_project/tome/TOME_process_data/TOME_Mmatrix4.Final.RData")
#
#
TOME_M <- CreateSeuratObject(counts = TOME_Mmatrix4, project = "TOME_M", min.cells = 0, min.features = 0)
TOME_M
#
head(TOME_M.metadata)
TOME_M.metadata2 <- TOME_M.metadata[co_cells,]
dim(TOME_M.metadata2)
head(Nature_annotation)
Nature_annotation2 <- Nature_annotation[co_cells,]
dim(Nature_annotation2)
#
#
#
TOME_M$NG.cell_state <- TOME_M.metadata2$cell_state
TOME_M$NG.cell_type <- TOME_M.metadata2$cell_type
TOME_M$Nature.all_read_count <- Nature_annotation2$all_read_count
TOME_M$Nature.Total_mRNAs <- Nature_annotation2$Total_mRNAs
TOME_M$development_stage <- Nature_annotation2$development_stage
TOME_M$Nature.num_genes_expressed <- Nature_annotation2$num_genes_expressed
#
TOME_M$Size_Factor <- Nature_annotation2$Size_Factor
TOME_M$Main_Cluster <- Nature_annotation2$Main_Cluster
TOME_M$Main_cluster_tsne_1 <- Nature_annotation2$Main_cluster_tsne_1
TOME_M$Main_cluster_tsne_2 <- Nature_annotation2$Main_cluster_tsne_2
TOME_M$Sub_cluster <- Nature_annotation2$Sub_cluster
TOME_M$Sub_cluster_tsne_1 <- Nature_annotation2$Sub_cluster_tsne_1
TOME_M$Sub_cluster_tsne_2 <- Nature_annotation2$Sub_cluster_tsne_2
TOME_M$sub_cluster_id <- Nature_annotation2$sub_cluster_id
TOME_M$Main_cell_type <- Nature_annotation2$Main_cell_type
TOME_M$Main_trajectory <- Nature_annotation2$Main_trajectory
TOME_M$Main_trajectory_umap_1 <- Nature_annotation2$Main_trajectory_umap_1
TOME_M$Main_trajectory_umap_2 <- Nature_annotation2$Main_trajectory_umap_2
TOME_M$Main_trajectory_umap_3 <- Nature_annotation2$Main_trajectory_umap_3
TOME_M$Main_trajectory_refined_by_cluster <- Nature_annotation2$Main_trajectory_refined_by_cluster
TOME_M$Main_trajectory_refined_umap_1 <- Nature_annotation2$Main_trajectory_refined_umap_1
TOME_M$Main_trajectory_refined_umap_2 <- Nature_annotation2$Main_trajectory_refined_umap_2
TOME_M$Main_trajectory_refined_umap_3 <- Nature_annotation2$Main_trajectory_refined_umap_3
TOME_M$Sub_trajectory_name <- Nature_annotation2$Sub_trajectory_name
TOME_M$Sub_trajectory_umap_1 <- Nature_annotation2$Sub_trajectory_umap_1
TOME_M$Sub_trajectory_umap_2 <- Nature_annotation2$Sub_trajectory_umap_2
TOME_M$Sub_trajectory_louvain_component <- Nature_annotation2$Sub_trajectory_louvain_component
TOME_M$Sub_trajectory_Pseudotime <- Nature_annotation2$Sub_trajectory_Pseudotime
TOME_M$Barcode <- Nature_annotation2$Barcode
#
TOME_M@active.ident <- as.factor(TOME_M$Main_cell_type)
TOME_M.V2 <- subset(TOME_M,idents=as.character(as.data.frame(table(TOME_M$Main_cell_type))$Var1))
TOME_M.V2
#
table(TOME_M.V2@active.ident)
TOME_M.V2@active.ident <- as.factor(TOME_M.V2$Main_cell_type)
#
#
TOME_M.V2$development_stage2 <- as.factor(paste("E",TOME_M.V2$development_stage,sep = "_"))
TOME_M.V2@active.ident <- as.factor(TOME_M.V2$development_stage2)
table(TOME_M.V2@active.ident)
TOME_M.V2.E10.5 <- subset(TOME_M.V2,idents="E_10.5")
TOME_M.V2.E11.5 <- subset(TOME_M.V2,idents="E_11.5")
TOME_M.V2.E12.5 <- subset(TOME_M.V2,idents="E_12.5")
TOME_M.V2.E13.5 <- subset(TOME_M.V2,idents="E_13.5")
#
TOME_M.V2.E10.5@active.ident <- as.factor(TOME_M.V2.E10.5$Main_cell_type)
TOME_M.V2.E11.5@active.ident <- as.factor(TOME_M.V2.E11.5$Main_cell_type)
TOME_M.V2.E12.5@active.ident <- as.factor(TOME_M.V2.E12.5$Main_cell_type)
TOME_M.V2.E13.5@active.ident <- as.factor(TOME_M.V2.E13.5$Main_cell_type)
#
TOME_M.V2.E10.5.downsample <- subset(TOME_M.V2.E10.5,downsample= 75*5)
TOME_M.V2.E11.5.downsample <- subset(TOME_M.V2.E11.5,downsample= 75*5)
TOME_M.V2.E12.5.downsample <- subset(TOME_M.V2.E12.5,downsample= 75*5)
TOME_M.V2.E13.5.downsample <- subset(TOME_M.V2.E13.5,downsample= 75*5)
#
TOME_M.V2.downsample <- merge(TOME_M.V2.E10.5.downsample,y=c(TOME_M.V2.E11.5.downsample,TOME_M.V2.E12.5.downsample,TOME_M.V2.E13.5.downsample))
TOME_M.V2.downsample
TOME_M.V2.downsample@active.ident <- as.factor(TOME_M.V2.downsample$Main_cell_type)
table(TOME_M.V2.downsample@active.ident)
##
TOME_M.V2.downsample.33 <- subset(TOME_M.V2.downsample,idents = c("Lens","Stromal cells","Epithelial","Megakaryocytes","Neutrophils"),invert=T)
TOME_M.V2.downsample.33
TOME_M.V2.downsample.33@active.ident <- as.factor(paste(TOME_M.V2.downsample.33$Main_cell_type,TOME_M.V2.downsample.33$development_stage2,sep = "#"))
table(TOME_M.V2.downsample.33@active.ident)
#
rm(TOME_M.V2.E10.5)
rm(TOME_M.V2.E11.5)
rm(TOME_M.V2.E12.5)
rm(TOME_M.V2.E13.5)
gc()
#
#
TOME_M.V2.downsample.33
TOME_M.V2.downsample.33$pseudoID <- TOME_M.V2.downsample.33@active.ident
TOME_M.V2.downsample.33.pseudo <- PseudoCell(TOME_M.V2.downsample.33, "RNA","data","pseudoID",5)
#
save(TOME_M.V2.downsample,file = "/home/sunkeyong/MOPA_project/tome/TOME_process_data/33.pseudo/TOME_M.V2.downsample.RData")
save(TOME_M.V2.downsample.33,file = "/home/sunkeyong/MOPA_project/tome/TOME_process_data/33.pseudo/TOME_M.V2.downsample.33.RData")
save(TOME_M.V2.downsample.33.pseudo,file = "/home/sunkeyong/MOPA_project/tome/TOME_process_data/33.pseudo/TOME_M.V2.downsample.33.pseudo.RData")
#
#
setwd("/home/sunkeyong/MOPA_project/Figure3")
ArrowFiles <- c("/home/sunkeyong/MOPA_project/fragment_total/E1A.arrow",
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
#
MOPA <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  copyArrows = TRUE)
MOPA
load("../Figure1/MOPA.seurat_type")
MOPA.seurat_type@active.ident <- as.factor(MOPA.seurat_type$ID)
cluster.rank.new <- c("Neural Tube","Notochord cells","Neural progenitor cells", "Radial glia","Postmitotic premature neurons",
                      "Sensory neurons","Oligodendrocyte Progenitors","Premature oligodendrocyte","Schwann cell precursor",
                      "Inhibitory neuron progenitors","Inhibitory interneurons","Inhibitory neurons","Excitatory neurons",
                      "Cholinergic neurons","Granule neurons","Ependymal cell","Isthmic organizer cells","Connective tissue progenitors",
                      "Chondrocytes & osteoblasts","Chondroctye progenitors","Intermediate Mesoderm","Limb mesenchyme",
                      "Jaw and tooth progenitors","Early mesenchyme","Myocytes","Cardiac muscle lineages","Epithelial cells",
                      "Melanocytes","Hepatocytes","Endothelial cells","White blood cells",
                      "Definitive erythroid lineage","Primitive erythroid lineage")
levels(MOPA.seurat_type) <-  cluster.rank.new
MOPA.seurat_type <- subset(MOPA.seurat_type,downsample=300)
MOPA.seurat_type
#
MOPA <- MOPA[Cells(MOPA.seurat_type),]
MOPA
#
load("../Figure2/MOPA_peakset.RData")
MOPA <- addPeakSet(MOPA,peakSet =MOPA_peakset,force=T)
MOPA <- addPeakMatrix(MOPA,force = T)
getAvailableMatrices(MOPA)
#
MOPA <- addIterativeLSI(
  ArchRProj = MOPA,
  useMatrix = "PeakMatrix", 
  name = "PeakMatrix_LSI", 
  iterations = 3, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(3,4), 
    sampleCells = 10000, 
    n.start = 10,
    maxClusters = 100
  ), 
  varFeatures = 150000, 
  dimsToUse = 1:100,
  force=T
)
#
MOPA <- addImputeWeights(MOPA, reducedDims = "PeakMatrix_LSI")
MOPA <- addUMAP(
  ArchRProj = MOPA, 
  reducedDims = "PeakMatrix_LSI", 
  name = "PeakMatrix_LSI_UMAP1", 
  nNeighbors = 30, 
  dimsToUse = 1:100, 
  minDist = 0.3, 
  metric = "cosine",force = T
)
#
addArchRThreads(threads = 1)
#
#
MOPA <- addGeneIntegrationMatrix(
  ArchRProj = MOPA, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GIMsc",
  reducedDims = "PeakMatrix_LSI",
  seRNA = TOME_M.V2.downsample.33,
  addToArrow = F,
  force = TRUE,
  dimsToUse = 1:100,
  UMAPParams = list(n_neighbors = 40, min_dist = 0.2, metric = "cosine", verbose = FALSE),
  nGenes = 3000,
  sampleCellsATAC = 50000,
  sampleCellsRNA = 50000,
  groupRNA = "Main_cell_type",
  nameCell = "predictedCell_sc",
  nameGroup = "predictedGroup_sc",
  nameScore = "predictedScore_sc"
)
#
MOPA <- addGeneIntegrationMatrix(
  ArchRProj = MOPA, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GIMpseu",
  reducedDims = "PeakMatrix_LSI",
  seRNA = TOME_M.V2.downsample.33.pseudo,
  addToArrow = F,
  force = TRUE,
  dimsToUse = 1:100,
  UMAPParams = list(n_neighbors = 40, min_dist = 0.2, metric = "cosine", verbose = FALSE),
  nGenes = 3000,
  sampleCellsATAC = 50000,
  sampleCellsRNA = 50000,
  groupRNA = "Main_cell_types",
  nameCell = "predictedCell_pseu",
  nameGroup = "predictedGroup_pseu",
  nameScore = "predictedScore_pseu"
)
save.image("pro.filter.RData")
#
pseu1 <- readRDS("ArchROutput/RNAIntegration/GIMpseu1/Save-Block1-JointCCA.rds")
sc1 <- readRDS("ArchROutput/RNAIntegration/GIMsc1/Save-Block1-JointCCA.rds")
#
#
pseu1.ATAC <- subset(pseu1,Assay=="ATAC")
sc1.ATAC <- subset(sc1,Assay=="ATAC")
#
hist(pseu1$Score)
hist(sc1.ATAC$Score)
#
midu1 <-as.data.frame(pseu1.ATAC$Score)
midu1$method  <- "pseudocell"
colnames(midu1) <- c("predicted_score","method")
midu3 <- as.data.frame(sc1.ATAC$Score)
midu3$method  <- "singlecell"
colnames(midu3) <- c("predicted_score","method")
head(midu1)
#
midu <- rbind(midu1,midu3)
head(midu)
#
p<-ggplot(midu, aes(x = predicted_score))
p + geom_density(color = "black", fill = "gray")
p + geom_density(aes(color = method))
p + geom_density(aes(fill = method), alpha=0.4)
#
median(pseu1.ATAC$Score) # 0.8946582
median(sc1.ATAC$Score) # 0.7433351
mean(pseu1.ATAC$Score)# 0.8157211
mean(sc1.ATAC$Score) # 0.7235369
#
median.values <- as.data.frame(c(0.7433351,0.8946582))
colnames(median.values) <- "median"
median.values
median.values$method <- c("singlecell","pseudobulk")
median.values
p + geom_density(aes(fill = method), alpha=0.6)+
  geom_vline(data = median.values, aes(xintercept = median, color=method),linetype="dashed")+
  #theme_classic()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
#
#
unique(projEu$emid)
unique(TOME_M.V2.downsample.300.33$Main_cell_type)
#
intersect(unique(projEu$emid),unique(TOME_M.V2.downsample.300.F$Main_cell_type))
setdiff(unique(projEu$emid),unique(projEu$predictedGroup_pseu1))
#
library(ggplot2)
library(cowplot)
#
table(projEu$emid==projEu$predictedGroup_pseu1) # 1-2553/9900 = 0.7421212
predictions  <- table(projEu$emid,projEu$predictedGroup_pseu1)
dim(predictions)
head(predictions)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- predictions[rev(cluster.rank.new),rev(cluster.rank.new)]
head(predictions)
head(predictions)
predictions <- as.data.frame(predictions)
head(predictions)
dim(predictions)
p1.pseu1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
                                                                                                  low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p1.pseu1
#
write.csv(predictions,file = "Source_Data/Fig4C.heatmap.csv")# 20231015
#
table(projEu$emid==projEu$predictedGroup_sc1) ## 1-2908/9900 = 0.7062626
predictions  <- table(projEu$emid,projEu$predictedGroup_sc1)
dim(predictions)
head(predictions)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- predictions[rev(cluster.rank.new),rev(cluster.rank.new)]
predictions <- as.data.frame(predictions)
head(predictions)
p1.sc1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
                                                                                                low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
  theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
p1.sc1
write.csv(predictions,file = "Source_Data/FigS4F.heatmap.csv")# 20231015
#
#
#
# plot 2D umap 
pseu1.ATAC <- subset(pseu1,Assay=="ATAC")
sc1.ATAC <- subset(sc1,Assay=="ATAC")
#
pseu1.RNA <- subset(pseu1,Assay=="RNA")
sc1.RNA <- subset(sc1,Assay=="RNA")
#
rownames(pseu1.ATAC) <- unlist(strsplit(rownames(pseu1.ATAC),"_query"))
rownames(sc1.ATAC) <- unlist(strsplit(rownames(sc1.ATAC),"_query"))
#
meta.em2 <- meta.em[rownames(pseu1.ATAC),]
dim(meta.em2)
meta.em3 <- meta.em[rownames(sc1.ATAC),]
#
pseu1.ATAC$Group <- meta.em2$emID.F2
sc1.ATAC$Group <- meta.em3$emID.F2
#
pseu1 <- rbind(pseu1.RNA,pseu1.ATAC)
dim(pseu1)
sc1 <- rbind(sc1.RNA,sc1.ATAC)
dim(sc1)
#
#
#
x.num <- c(rep(1,800000),rep(0,800000),rep(2,800000))
head(x.num)
x.num.sample <- cbind(sample(x.num,nrow(pseu1)),
                      sample(x.num,nrow(pseu1)),
                      sample(x.num,nrow(pseu1)),
                      sample(x.num,nrow(pseu1)),
                      sample(x.num,nrow(pseu1)))
head(x.num.sample)
colnames(x.num.sample) <- c("F1","F2","F3","F4","F5")
rownames(x.num.sample) <- rownames(pseu1)
#
pseu1.seurat <- CreateSeuratObject(counts = t(x.num.sample), assay = "ATAC", project = "CardiacMuscle", min.cells = 0, min.features = 0)
pseu1.seurat
pseu1.seurat <- NormalizeData(pseu1.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
pseu1.seurat <- FindVariableFeatures(pseu1.seurat, selection.method = "vst", nfeatures = 4)
pseu1.seurat <- ScaleData(pseu1.seurat, features = rownames(pseu1.seurat))
variablegene <- VariableFeatures(object = pseu1.seurat)
pseu1.seurat <- RunPCA(pseu1.seurat, features = variablegene,npcs =2)
DimPlot(pseu1.seurat)
#
s.UMAP <- cbind(pseu1$UMAP1,pseu1$UMAP2)
#
head(s.UMAP )
rownames(s.UMAP ) <- rownames(pseu1)
head(s.UMAP )
colnames(s.UMAP ) <- c("PC_1","PC_2")
head(s.UMAP )
class(s.UMAP )
#s.UMAP <- as.data.frame(s.UMAP)
pseu1.seurat@reductions$pca@cell.embeddings <- s.UMAP 
#DimPlot(seurat.projEu,label = T,repel = T)+NoLegend()
DimPlot(pseu1.seurat,label = T,repel = T,pt.size = 0.2)+NoLegend()
#
##
x.num <- c(rep(1,800000),rep(0,800000),rep(2,800000))
head(x.num)
x.num.sample <- cbind(sample(x.num,nrow(sc1)),
                      sample(x.num,nrow(sc1)),
                      sample(x.num,nrow(sc1)),
                      sample(x.num,nrow(sc1)),
                      sample(x.num,nrow(sc1)))
head(x.num.sample)
colnames(x.num.sample) <- c("F1","F2","F3","F4","F5")
rownames(x.num.sample) <- rownames(sc1)
#
sc1.seurat <- CreateSeuratObject(counts = t(x.num.sample), assay = "ATAC", project = "CardiacMuscle", min.cells = 0, min.features = 0)
sc1.seurat
sc1.seurat <- NormalizeData(sc1.seurat, normalization.method = "LogNormalize", scale.factor = 10000)
sc1.seurat <- FindVariableFeatures(sc1.seurat, selection.method = "vst", nfeatures = 4)
sc1.seurat <- ScaleData(sc1.seurat, features = rownames(sc1.seurat))
variablegene <- VariableFeatures(object = sc1.seurat)
sc1.seurat <- RunPCA(sc1.seurat, features = variablegene,npcs =2)
DimPlot(sc1.seurat)
#
s.UMAP <- cbind(sc1$UMAP1,sc1$UMAP2)
#
head(s.UMAP )
rownames(s.UMAP ) <- rownames(sc1)
head(s.UMAP )
colnames(s.UMAP ) <- c("PC_1","PC_2")
head(s.UMAP )
class(s.UMAP )
#s.UMAP <- as.data.frame(s.UMAP)
sc1.seurat@reductions$pca@cell.embeddings <- s.UMAP 
#DimPlot(seurat.projEu,label = T,repel = T)+NoLegend()
DimPlot(sc1.seurat,label = T,repel = T,pt.size = 0.2)+NoLegend()
#
pseu1.seurat$celltype <- pseu1$Group
sc1.seurat$celltype <- sc1$Group
#
pseu1.seurat$tech <- pseu1$Assay
sc1.seurat$tech <- sc1$Assay
#
pseu1.seurat@active.ident <- as.factor(pseu1.seurat$celltype)
levels(pseu1.seurat) <- cluster.rank.new
DimPlot(pseu1.seurat,cols = color.use.em,pt.size = 0.000001,repel = T,label = T)+NoLegend()
p <- DimPlot(pseu1.seurat,cols = color.use.em,pt.size = 0.000001,repel = T,label = F)+NoLegend()+
  theme(plot.title = element_text(color="black", size=14))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
#
ggsave("pseu1.coembedding.20230410.png",plot=p,device="png",dpi=400,units = "cm",width = 13,height = 13)
#
sc1.seurat@active.ident <- as.factor(sc1.seurat$celltype)
levels(sc1.seurat) <- cluster.rank.new
DimPlot(sc1.seurat,cols = color.use.em,pt.size = 0.000001,repel = T,label = T)+NoLegend()
p <- DimPlot(sc1.seurat,cols = color.use.em,pt.size = 0.000001,repel = T,label = F)+NoLegend()+
  theme(plot.title = element_text(color="black", size=14))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
#
ggsave("sc1.coembedding.20230410.png",plot=p,device="png",dpi=400,units = "cm",width = 13,height = 13)
#
#
#
pseu1.seurat@active.ident <- as.factor(pseu1.seurat$tech)
levels(pseu1.seurat) <- c("RNA","ATAC")
DimPlot(pseu1.seurat,pt.size = 0.000001,cols =c("#FF6347","#40E0D0"),shuffle=T,repel = T,label = T)
p <- DimPlot(pseu1.seurat,cols =c("#40E0D0","#FF6347"),shuffle=T,pt.size = 0.000001,repel = T,label = F)+NoLegend()+
  theme(plot.title = element_text(color="black", size=14))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
#
ggsave("pseu1.coembedding.tech.20230410.png",plot=p,device="png",dpi=400,units = "cm",width = 13,height = 13)
#
sc1.seurat@active.ident <- as.factor(sc1.seurat$tech)
levels(pseu1.seurat) <- c("RNA","ATAC")
DimPlot(sc1.seurat,cols =c("#40E0D0","#FF6347"),shuffle=T,pt.size = 0.000001,repel = T,label = T)+NoLegend()
p <- DimPlot(sc1.seurat,cols =c("#40E0D0","#FF6347"),shuffle=T,pt.size = 0.000001,repel = T,label = F)+NoLegend()+
  theme(plot.title = element_text(color="black", size=14))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
#
ggsave("sc1.coembedding.tech.20230410.png",plot=p,device="png",dpi=400,units = "cm",width = 13,height = 13)
#
TOME_M.V2.downsample.300.33@active.ident <- as.factor(TOME_M.V2.downsample.300.33$Main_cell_type)
TOME_M.V2.downsample.33.pseudo.F@active.ident <- as.factor(TOME_M.V2.downsample.33.pseudo.F$Main_cell_types)
levels(TOME_M.V2.downsample.300.33) <- cluster.rank.new
levels(TOME_M.V2.downsample.33.pseudo.F) <- cluster.rank.new
#
VlnPlot(TOME_M.V2.downsample.300.33,features = c("nCount_RNA"),cols = color.use.em ,pt.size = 0,y.max = 4000)+NoLegend()
VlnPlot(TOME_M.V2.downsample.300.33,features = c("nFeature_RNA"),cols = color.use.em ,pt.size = 0,y.max = 2500)+NoLegend()
VlnPlot(TOME_M.V2.downsample.33.pseudo.F,features = c("nCount_RNA"),cols = color.use.em ,pt.size = 0,y.max = 4000)+NoLegend()
VlnPlot(TOME_M.V2.downsample.33.pseudo.F,features = c("nFeature_RNA"),cols = color.use.em ,pt.size = 0,y.max = 4500)+NoLegend()
#
#
DimPlot(TOME_M.V2.downsample.300.33,cols = color.use.em,pt.size = 0.000001,repel = T,label = T)+NoLegend()
p <- DimPlot(TOME_M.V2.downsample.300.33,cols = color.use.em,pt.size = 0.000001,repel = T,label = F)+NoLegend()+
  theme(plot.title = element_text(color="black", size=14))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
#
ggsave("TOME_M.V2.downsample.300.33.20230410.png",plot=p,device="png",dpi=400,units = "cm",width = 13,height = 13)
#
DimPlot(TOME_M.V2.downsample.33.pseudo.F,cols = color.use.em,pt.size = 0.000001,repel = T,label = T)+NoLegend()
p <- DimPlot(TOME_M.V2.downsample.33.pseudo.F,cols = color.use.em,pt.size = 0.000001,repel = T,label = F)+NoLegend()+
  theme(plot.title = element_text(color="black", size=14))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.x = element_blank(),axis.title.y = element_blank())+
  theme(axis.text.x=element_blank(),axis.text.y=element_blank())+
  theme(axis.ticks.y=element_blank(),axis.ticks.x=element_blank())+
  theme(axis.line.x=element_blank(),axis.line.y=element_blank())+
  theme(panel.background = element_rect( colour = "black", size = 0))+
  theme(panel.border = element_blank())
p
#
ggsave("TOME_M.V2.downsample.33.pseudo.F.20230410.png",plot=p,device="png",dpi=400,units = "cm",width = 13,height = 13)
#
save.image("/home/sunkeyong/MOPA_revision1/Fig1_33MainCellTypes/33Main_cell_types.RData")
#
#
#
### define the subtypes 
## taking White blood cells as an examples
# processing scRNA-seq data
setwd("/home/sunkeyong/MOPA_project/subcluster/onebyone/Whitebloodcell")
### Download TOME MOCA scRNA-seq and load
load("MOCA.RData.RData")
MOCA
MOCA@active.ident <- as.factor(MOCA$Main_cell_type)
table(MOCA@active.ident)
MOCA.blood <- subset(MOCA,idents = "White blood cells")
rm(MOCA)
#
MOCA.blood  <- NormalizeData(MOCA.blood , normalization.method = "LogNormalize", scale.factor = 10000)
MOCA.blood  <- FindVariableFeatures(MOCA.blood , selection.method = "vst", nfeatures = 3000)
MOCA.blood  <- ScaleData(MOCA.blood , features = rownames(MOCA.blood ))
variablegene <- VariableFeatures(object = MOCA.blood )
MOCA.blood  <- RunPCA(MOCA.blood , features = variablegene)
ElbowPlot(MOCA.blood,ndims = 50)
MOCA.blood  <- FindNeighbors(MOCA.blood , dims = 1:30)
MOCA.blood  <- FindClusters(MOCA.blood , resolution = 1)
MOCA.blood  <- RunUMAP(MOCA.blood , dims = 1:30)
MOCA.blood  <- RunTSNE(MOCA.blood , dims = 1:30)
DimPlot(MOCA.blood ,reduction = "tsne",label=T,repel=T)
DimPlot(MOCA.blood ,reduction = "umap",label=T,repel=T)
#
table(MOCA.blood$seurat_clusters,MOCA.blood$development_stage)
#
MOCA.blood$development_stage <- MOCA.blood.raw$development_stage
MOCA.blood$sub_cluster_id <- MOCA.blood.raw$sub_cluster_id
#
DimPlot(MOCA.blood,reduction = "umap",label=T,repel=T,group.by = "development_stage")
DimPlot(MOCA.blood,reduction = "umap",label=T,repel=T,group.by = "sub_cluster_id")
#
library(dplyr)
MOCA.blood.markers <- FindAllMarkers(MOCA.blood, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MOCA.blood.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top10 <- MOCA.blood.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(MOCA.blood, features = top10$gene) + NoLegend()
#
VlnPlot(MOCA.blood,features = c("Mki67"),group.by = "development_stage")
VlnPlot(MOCA.blood,features = c("Pcna"),group.by = "development_stage")
VlnPlot(MOCA.blood,features = c("Top2a"),group.by = "development_stage")
VlnPlot(MOCA.blood,features = c("Mcm6"),group.by = "development_stage")
#
#
MOCA.blood <- RenameIdents(MOCA.blood, '0' = 'Promonocyte')
MOCA.blood <- RenameIdents(MOCA.blood, '1' = 'Macrophage')
MOCA.blood <- RenameIdents(MOCA.blood, '2' = 'Macrophage')
MOCA.blood <- RenameIdents(MOCA.blood, '3' = 'monocyte')
MOCA.blood <- RenameIdents(MOCA.blood, '4' = 'Kupffer cell')
MOCA.blood <- RenameIdents(MOCA.blood, '5' = 'HSPC_2')
MOCA.blood <- RenameIdents(MOCA.blood, '6' = 'Microglia')
MOCA.blood <- RenameIdents(MOCA.blood, '7' = 'Basophil/mast')
MOCA.blood <- RenameIdents(MOCA.blood, '8' = 'HSPC_1 (Liver)')
MOCA.blood <- RenameIdents(MOCA.blood, '9' = 'HSPC_3')
MOCA.blood <- RenameIdents(MOCA.blood, '10' = 'Low-quality')
MOCA.blood <- RenameIdents(MOCA.blood, '11' = 'Proerythroblast')
MOCA.blood <- RenameIdents(MOCA.blood, '12' = 'Lymphoid lineage')
MOCA.blood <- RenameIdents(MOCA.blood, '13' = 'Lymphoid lineage')
MOCA.blood <- RenameIdents(MOCA.blood, '14' = 'granulocyte')
#
DimPlot(MOCA.blood,label = T,repel = T)+NoLegend()
#
MOCA.blood$ID <- MOCA.blood@active.ident
#
#
ArrowFiles <- c("/home/sunkeyong/MOPA_project/Major_cell_types_arrowfile/White.blood.cells.arrow")
#
MOPA.blood.ArchR <- ArchRProject(ArrowFiles = ArrowFiles, copyArrows = TRUE)
#
MOPA.blood.ArchR <- addIterativeLSI(
  ArchRProj = MOPA.blood.ArchR,
  useMatrix = "PeakMatrix", 
  name = "PeakMatrix_LSI", 
  iterations = 3, 
  clusterParams = list(
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
MOPA.blood.ArchR <- addClusters(
  input = MOPA.blood.ArchR,
  reducedDims = "PeakMatrix_LSI",
  method = "Seurat",
  name = "PeakMatrix_LSI_R0.8",
  resolution = 3, 
  dimsToUse = 1:100,
  maxClusters= 300,
  force = T,
  sampleCells = NULL
)
#
MOPA.blood.ArchR <- addUMAP(
  ArchRProj = MOPA.blood.ArchR, 
  reducedDims = "PeakMatrix_LSI", 
  name = "PeakMatrix_LSI_UMAP_1", 
  nNeighbors = 100, 
  dimsToUse = 1:100, 
  minDist = 0.4, 
  metric = "cosine",force = T
)
#
plotEmbedding(ArchRProj = MOPA.blood.ArchR, colorBy = "cellColData", name = "Sample", embedding = "PeakMatrix_LSI_UMAP_1")
plotEmbedding(ArchRProj = MOPA.blood.ArchR, colorBy = "cellColData", name = "PeakMatrix_LSI_R0.8", embedding = "PeakMatrix_LSI_UMAP_1")
#
###### create seurat file
x.num <- c(rep(1,800000),rep(0,800000),rep(2,800000))
head(x.num)
x.num.sample <- cbind(sample(x.num,length(MOPA.blood.ArchR$cellNames)),
                      sample(x.num,length(MOPA.blood.ArchR$cellNames)),
                      sample(x.num,length(MOPA.blood.ArchR$cellNames)),
                      sample(x.num,length(MOPA.blood.ArchR$cellNames)),
                      sample(x.num,length(MOPA.blood.ArchR$cellNames)))
head(x.num.sample)
colnames(x.num.sample) <- c("F1","F2","F3","F4","F5")
rownames(x.num.sample) <- MOPA.blood.ArchR$cellNames
#
MOPA.blood.Seurat <- CreateSeuratObject(counts = t(x.num.sample), assay = "ATAC", project = "CardiacMuscle", min.cells = 0, min.features = 0)
MOPA.blood.Seurat
MOPA.blood.Seurat <- NormalizeData(MOPA.blood.Seurat, normalization.method = "LogNormalize", scale.factor = 10000)
MOPA.blood.Seurat <- FindVariableFeatures(MOPA.blood.Seurat, selection.method = "vst", nfeatures = 4)
MOPA.blood.Seurat <- ScaleData(MOPA.blood.Seurat, features = rownames(MOPA.blood.Seurat))
variablegene <- VariableFeatures(object = MOPA.blood.Seurat)
MOPA.blood.Seurat <- RunPCA(MOPA.blood.Seurat, features = variablegene,npcs =2)
DimPlot(MOPA.blood.Seurat)
#
x <- MOPA.blood.ArchR@embeddings$PeakMatrix_LSI_UMAP_1
y <- x$df
s.UMAP <- cbind(y$`PeakMatrix_LSI#UMAP_Dimension_1`,y$`PeakMatrix_LSI#UMAP_Dimension_2`)
#
head(s.UMAP )
rownames(s.UMAP ) <- rownames(y)
head(s.UMAP )
colnames(s.UMAP ) <- c("PC_1","PC_2")
head(s.UMAP )
class(s.UMAP )
#s.UMAP <- as.data.frame(s.UMAP)
MOPA.blood.Seurat@reductions$pca@cell.embeddings <- s.UMAP 
#DimPlot(seurat.MOPA.blood.ArchR,label = T,repel = T)+NoLegend()
DimPlot(MOPA.blood.Seurat,label = T,repel = T,pt.size = 0.2)+NoLegend()
#
#
###
MOPA.blood.Seurat$ID1 <- MOPA.blood.ArchR$PeakMatrix_LSI_R1
MOPA.blood.Seurat$ID2 <- MOPA.blood.ArchR$PeakMatrix_LSI_R1.3
levels(MOPA.blood.Seurat)
DimPlot(MOPA.blood.Seurat,label = T,repel = T,group.by = "ID1",pt.size = 0.2)
DimPlot(MOPA.blood.Seurat,label = T,repel = T,group.by = "ID2",pt.size = 0.2)
#
MOPA.blood.Seurat@active.ident <- as.factor(MOPA.blood.Seurat$ID2)
MOPA.blood.Seurat <- RenameIdents(MOPA.blood.Seurat, "C6"="C2")
MOPA.blood.Seurat <- RenameIdents(MOPA.blood.Seurat, "C7"="C6")
#
levels(MOPA.blood.Seurat) <- c("C1","C2","C3","C4","C5","C6")
MOPA.blood.Seurat$ID <- MOPA.blood.Seurat@active.ident
DimPlot(MOPA.blood.Seurat,label = T,repel = T,pt.size = 0.2)
#
MOPA.blood.Seurat.Whitebloodcell <- MOPA.blood.Seurat
#
#
MOCA.blood.pse1 <- PseudoCell(MOCA.blood, "RNA","data","ID",5)
MOCA.blood.pse2 <- PseudoCell(MOCA.blood, "RNA","data","ID",5)
MOCA.blood.pseudo <- merge(MOCA.blood.pse1,MOCA.blood.pse2)
MOCA.blood.pseudo 
MOCA.blood.pseudo$pseudo.ID <- unlist(strsplit(Cells(mypbmc),"_"))[seq(1,3*2744,3)]
#
MOCA.blood.pseudo  <- NormalizeData(MOCA.blood.pseudo , normalization.method = "LogNormalize", scale.factor = 10000)
MOCA.blood.pseudo  <- FindVariableFeatures(MOCA.blood.pseudo , selection.method = "vst", nfeatures = 3000)
MOCA.blood.pseudo  <- ScaleData(MOCA.blood.pseudo , features = rownames(MOCA.blood.pseudo ))
variablegene <- VariableFeatures(object = MOCA.blood.pseudo )
MOCA.blood.pseudo  <- RunPCA(MOCA.blood.pseudo , features = variablegene)
ElbowPlot(MOCA.blood.pseudo,ndims = 50)
MOCA.blood.pseudo  <- FindNeighbors(MOCA.blood.pseudo , dims = 1:30)
MOCA.blood.pseudo  <- FindClusters(MOCA.blood.pseudo , resolution = 1)
MOCA.blood.pseudo  <- RunUMAP(MOCA.blood.pseudo , dims = 1:30)
MOCA.blood.pseudo  <- RunTSNE(MOCA.blood.pseudo , dims = 1:30)
DimPlot(MOCA.blood.pseudo ,reduction = "tsne",label=T,repel=T)
DimPlot(MOCA.blood.pseudo ,reduction = "umap",label=T,repel=T)
#
MOPA.blood.ArchR <- addGeneIntegrationMatrix(
  ArchRProj = MOPA.blood.ArchR, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GIM.pseudo",
  reducedDims = "PeakMatrix_LSI",
  seRNA = MOCA.blood.pseudo,
  addToArrow = T,
  force = TRUE,
  dimsToUse = 1:50,
  UMAPParams = list(n_neighbors = 40, min_dist = 0.4, metric = "cosine", verbose = FALSE),
  nGenes = 3000,
  sampleCellsATAC = 10000,
  sampleCellsRNA = 10000,
  groupRNA = "define",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore")
#
plotEmbedding(ArchRProj = MOPA.blood.ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "PeakMatrix_LSI_UMAP_1")
#
getAvailableMatrices(MOPA.blood.ArchR)
Get.GIMpseudo.mat <- getMatrixFromProject(
  ArchRProj = MOPA.blood.ArchR,
  useMatrix = "GIM.pseudo")
rownames(Get.GIMpseudo.mat@assays@data$GIM.pseudo) <- Get.GIMpseudo.mat@elementMetadata$name
Get.GIMpseudo.mat.matrix <- Get.GIMpseudo.mat@assays@data$GIM.pseudo
head(Get.GIMpseudo.mat.matrix)
dim(Get.GIMpseudo.mat.matrix)
#
#
cogenes <- intersect(rownames(Get.GIMpseudo.mat.matrix),rownames(MOCA.blood.pseudo))
Get.GIMpseudo.mat.matrix.seurat <- CreateSeuratObject(Get.GIMpseudo.mat.matrix,assay = "RNA")
Get.GIMpseudo.mat.matrix.seurat
Get.GIMpseudo.mat.matrix.seurat$tech <- "ATAC"
#
MOPA.blood.ArchR.metadata <- as.data.frame(MOPA.blood.ArchR@cellColData)
MOPA.blood.ArchR.metadata <- MOPA.blood.ArchR.metadata[rownames(Get.GIMpseudo.mat.matrix.seurat@meta.data),]
Get.GIMpseudo.mat.matrix.seurat$ID <- MOPA.blood.ArchR.metadata$predictedGroup
Get.GIMpseudo.mat.matrix.seurat
#
coembed <- merge(x = MOCA.blood.pseudo, y = Get.GIMpseudo.mat.matrix.seurat)
coembed 
#
coembed.list <- SplitObject(coembed, split.by = "tech")
coembed.list
#
coembed.list <- lapply(X = coembed.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
#
features <- SelectIntegrationFeatures(object.list = coembed.list,nfeatures = 3000)
coembed.anchors <- FindIntegrationAnchors(object.list = coembed.list, anchor.features = features)
coembed.combined <- IntegrateData(anchorset = coembed.anchors)
coembed.combined <- ScaleData(coembed.combined, features = features, do.scale = T)
coembed.combined <- RunPCA(coembed.combined, features = features, verbose = FALSE)
coembedcombined <- RunUMAP(coembed.combined, dims = 1:50,min.dist = 0.1,n.neighbors=50)
DimPlot(coembed.combined, group.by = "tech")
DimPlot(coembed.combined, group.by = "cluster",label = T)+NoLegend()
#
DimPlot(coembed.combined, group.by = "cluster",label = T,repel = T,split.by = "tech")+NoLegend()
#
levels(coembed.combined)
coembed.combined@active.ident <- as.factor(immune.combined$ID)
levels(coembed.combined) 
levels(coembed.combined) <- c("Basophilmast","HPSC1Liver","HPSC2","HPSC3","Lymphoidlineage","Kupffer cell",
                             "Macrophage","Microglia","monocyte","Proerythroblast")
#
color.for.immune <- c("#4EAEF0","#228B22","#74DDD0","#FFA500","#008B8B","#8B4513","#8A2BE2","#7CFC00","#EE82EE","#6A5ACD")
DimPlot(immune.combined, label = T,cols =color.for.immune )
DimPlot(immune.combined, label = T,cols = c("#69E4D6","#ED6A51"),group.by = "tech" )
#
MOCA.blood.pseudo.markers <- FindAllMarkers(MOCA.blood.pseudo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
MOCA.blood.pseudo.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top10 <- MOCA.blood.pseudo.markers %>% group_by(cluster) %>% top_n(n = 66, wt = avg_log2FC)
length(unique(top10$gene))
DoHeatmap(MOCA.blood.pseudo, features = top10$gene) + NoLegend()
#
#
MOCA.blood.pseudo.averages <- AverageExpression(MOCA.blood.pseudo, slot = 'data', return.seurat = TRUE)
length(unique(top10$gene))
DoHeatmap(MOCA.blood.pseudo.averages, features = unique(top10$gene),draw.lines = F) + NoLegend()+
  scale_fill_gradientn(colors = paletteContinuous("horizonExtra"))
#
Get.GS.mat <- getMatrixFromProject(
  ArchRProj = MOPA.blood.ArchR,
  useMatrix = "GeneScoreMatrix")
rownames(Get.GS.mat@assays@data$GeneScoreMatrix) <- Get.GS.mat@elementMetadata$name
Get.GS.mat.matrix <- Get.GS.mat@assays@data$GeneScoreMatrix
head(Get.GS.mat.matrix)
dim(Get.GS.mat.matrix)
#
Get.GS.mat.matrix.seurat <-CreateSeuratObject(Get.GS.mat.matrix[cogenes,],assay = "GS")
Get.GS.mat.matrix.seurat$ID <- MOPA$predictedGroup
table(Get.GS.mat.matrix.seurat$ID )
Get.GS.mat.matrix.seurat
#
Get.GS.mat.matrix.seurat
Get.GS.mat.matrix.seurat  <- NormalizeData(Get.GS.mat.matrix.seurat , normalization.method = "LogNormalize", scale.factor = 10000)
Get.GS.mat.matrix.seurat  <- FindVariableFeatures(Get.GS.mat.matrix.seurat , selection.method = "vst", nfeatures = 3000)
Get.GS.mat.matrix.seurat  <- ScaleData(Get.GS.mat.matrix.seurat , features = rownames(Get.GS.mat.matrix.seurat))
#
Get.GS.mat.matrix.seurat@active.ident <- as.factor(Get.GS.mat.matrix.seurat$ID)
levels(Get.GS.mat.matrix.seurat) <- levels(coembed.combined)
DoHeatmap(Get.GS.mat.matrix.seurat, features = top10$gene) + NoLegend()
#
#
Get.GS.mat.matrix.seurat.averages <- AverageExpression(Get.GS.mat.matrix.seurat, slot = 'data', return.seurat = TRUE)
length(unique(top10$gene))
DoHeatmap(Get.GS.mat.matrix.seurat.averages, features = unique(top10$gene),draw.lines = F,group.colors=color.for.immune) + NoLegend()+
  scale_fill_gradientn(colors = paletteContinuous("blueYellow"))
###
###
MOCA.blood.pseudo.averages.sub <- GetAssayData(MOCA.blood.pseudo.averages, slot = "scale.data")[unique(top10$gene),]
pheatmap::pheatmap(cor(MOCA.blood.pseudo.averages.sub))
Get.GS.mat.matrix.seurat.averages.sub <- GetAssayData(Get.GS.mat.matrix.seurat.averages, slot = "scale.data")[unique(top10$gene),]
pheatmap::pheatmap(cor(Get.GS.mat.matrix.seurat.averages.sub))
##
colnames(MOCA.blood.pseudo.averages.sub) <- paste("RNA",colnames(MOCA.blood.pseudo.averages.sub),sep= "_")
colnames(Get.GS.mat.matrix.seurat.averages.sub) <- paste("ATAC",colnames(Get.GS.mat.matrix.seurat.averages.sub),sep= "_")
cc1 <- cbind(MOCA.blood.pseudo.averages.sub,Get.GS.mat.matrix.seurat.averages.sub)
cc1.t <- cor(cc1,method = c("spearman"))
dim(cc1.t)
head(cc1.t)
library(pheatmap)
pheatmap::pheatmap(cc1.t[10:18,1:9],border_color = NA,
                   cluster_rows = F,
                   cluster_cols = F,
                   color = colorRampPalette(rev(brewer.pal(n = 8, name = "RdBu")))(100))
#
save.image("MOPA.blood.RData")
#
#
