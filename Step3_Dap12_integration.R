
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)

#load in data from Cell Ranger or other counts data ====

#for loading Cell Ranger counts:
setwd("/athena/ganlab/scratch/lif4001/Dap12/DF_2ndRound")
WT_1 <- readRDS(file = "WT_1_singlets_PCA.rds")
WT_2 <- readRDS(file = "WT_2_singlets_PCA.rds")
WT_3 <- readRDS(file = "WT_3_singlets_PCA.rds")

DAP12KO_1 <- readRDS(file = "DAP12KO_1_singlets_PCA.rds")
DAP12KO_2 <- readRDS(file = "DAP12KO_2_singlets_PCA.rds")
DAP12KO_3 <- readRDS(file = "DAP12KO_3_singlets_PCA.rds")

DAP12WT_1 <- readRDS(file = "DAP12WT_1_singlets_PCA.rds")
DAP12WT_2 <- readRDS(file = "DAP12WT_2_singlets_PCA.rds")
DAP12WT_3 <- readRDS(file = "DAP12WT_3_singlets_PCA.rds")

setwd("/athena/ganlab/scratch/lif4001/Dap12/add_2more/integration")
WT <- c(WT_1, WT_2, WT_3)
anchors_WT <- FindIntegrationAnchors(object.list = WT, dims = 1:30)
WT_integrated <- IntegrateData(anchorset = anchors_WT, dims = 1:30)
rm(WT_1, WT_2, WT_3, WT)

DAP12KO <- c(DAP12KO_1, DAP12KO_2, DAP12KO_3)
anchors_DAP12KO <- FindIntegrationAnchors(object.list = DAP12KO, dims = 1:30)
DAP12KO_integrated <- IntegrateData(anchorset = anchors_DAP12KO, dims = 1:30)
rm(DAP12KO_1, DAP12KO_2, DAP12KO_3, DAP12KO)

DAP12WT <- c(DAP12WT_1, DAP12WT_2, DAP12WT_3)
anchors_DAP12WT <- FindIntegrationAnchors(object.list = DAP12WT, dims = 1:30)
DAP12WT_integrated <- IntegrateData(anchorset = anchors_DAP12WT, dims = 1:30)
rm(DAP12WT_1, DAP12WT_2, DAP12WT_3, DAP12WT)

Dap12 <- c(WT_integrated, DAP12KO_integrated, DAP12WT_integrated)
anchors_Dap12 <- FindIntegrationAnchors(object.list = Dap12, dims = 1:30)
Dap12_integrated <- IntegrateData(anchorset = anchors_Dap12, dims = 1:30)
rm(WT_integrated, DAP12KO_integrated, DAP12WT_integrated, Dap12)

#saveRDS(Dap12_integrated, file = "Dap12_integrated.rds")

DefaultAssay(Dap12_integrated) <- 'integrated'

# Dap12_integrated <- NormalizeData(Dap12_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
# Dap12_integrated <- FindVariableFeatures(Dap12_integrated, selection.method = "vst", nfeatures = 3000)

Dap12_integrated <- ScaleData(Dap12_integrated, verbose = FALSE)
Dap12_integrated <- RunPCA(Dap12_integrated, features = VariableFeatures(object = Dap12_integrated), verbose = FALSE)

Dap12_integrated <- FindNeighbors(Dap12_integrated, dims = 1:15)
Dap12_integrated <- FindClusters(Dap12_integrated, resolution = 0.1)
Dap12_integrated <- RunUMAP(Dap12_integrated, dims = 1: 15)

DefaultAssay(Dap12_integrated) <- 'RNA'
Dap12_integrated <- NormalizeData(Dap12_integrated, normalization.method = "LogNormalize", scale.factor = 10000)
Dap12_integrated <- ScaleData(Dap12_integrated, features = rownames(Dap12_integrated))

#saveRDS(Dap12_integrated, file = 'Dap12_integrated_PCA_0.1.rds')
#Dap12_integrated <- readRDS(file = "Dap12_integrated_PCA_0.1.rds")

Dap12_integrated$Condition <- factor(x = Dap12_integrated$Condition, levels = c("WT","DAP12WT","DAP12KO"))
Dap12_integrated$Sample_Name <- factor(x = Dap12_integrated$Sample_Name, levels = c("WT_1","WT_2","WT_3","DAP12WT_1","DAP12WT_2",
                                                                                    "DAP12WT_3","DAP12KO_1","DAP12KO_2","DAP12KO_3"))

pdf("Dap12_QC.pdf", width=9, height=4)
Idents(Dap12_integrated) <- "Condition"
VlnPlot(object = Dap12_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

pdf("Dap12_QC_Sample.pdf", width=12, height=4)
Idents(Dap12_integrated) <- "Sample_Name"
VlnPlot(object = Dap12_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()

Idents(Dap12_integrated) <- "seurat_clusters"
pdf("Dap12_integrated_umap.pdf", width=5, height=4)
DimPlot(Dap12_integrated, reduction = 'umap', label = T)
dev.off()
pdf("Dap12_integrated_umap_split_individual.pdf", width=8, height=7)
DimPlot(Dap12_integrated, reduction = "umap", split.by = "Sample_Name", label = T, ncol = 3)
dev.off()
pdf("Dap12_integrated_umap_split_Condition.pdf", width=10, height=3)
DimPlot(Dap12_integrated, reduction = "umap", split.by = "Condition", label = T, ncol = 3)
dev.off()

write.csv(table(Dap12_integrated$seurat_clusters, Dap12_integrated$Sample_Name), "Dap12_cell_counts_cluster_by_sample.csv")

saveRDS(Dap12_integrated, file = 'Dap12_integrated_PCA_0.1.rds')

DefaultAssay(Dap12_integrated) <- 'RNA'

Dap12_markers <- FindAllMarkers(Dap12_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1, test.use = "MAST")
write.csv(Dap12_markers, "Dap12_markers.csv")

Dap12_markers <- read.csv(file = "Dap12_markers.csv", header=T,row.names =1)
top5 <- Dap12_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top5$gene <- as.character(top5$gene)
pdf("Dap12_HeatMapTop5_0.1_new.pdf", width=24, height=16)
DoHeatmap(Dap12_integrated, features = top5$gene) + NoLegend()
dev.off()

#Add marker genes

sig_EN<-c("Snap25","Slc17a7", "Nrgn","Gad1", "Gad2","Plp1", "Mbp", "Mobp","Sntn","Aqp4", "Clu", "Aldoc", "Pla2g7","Cx3cr1", "P2ry12", "Csf1r",
          "Pdgfra", "Vcan", "Flt1","Vtn", "Igfbp7")
markers.to.plot <- as.matrix(sig_EN)
pdf("Dap12_annotation_combine.pdf", width=10, height=5)
DotPlot(object = Dap12_integrated, features = rev(x = markers.to.plot)) + RotatedAxis()
dev.off()

