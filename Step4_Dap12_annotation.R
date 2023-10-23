
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(dplyr)
library(cowplot)
library(reshape2)
library(MAST)
setwd("/athena/ganlab/.../integration")
Dap12_integrated <- readRDS("Dap12_integrated_PCA_0.1.rds")
#Remove cluster 5-unknown, mixed
Dap12_integrated <- subset(Dap12_integrated, idents=c("5"), invert=T)
Idents(Dap12_integrated) <- "seurat_clusters"
Dap12_integrated <- RenameIdents(Dap12_integrated,
                                 `0` = "oligodendrocytes", `1`="excitatory neurons", `2`="excitatory neurons", `3`="astrocytes",
                                 `4`="excitatory neurons", `6`="inhibitory neurons", `7`="microglia",
                                 `8`="excitatory neurons", `9`="excitatory neurons", `10`="vascular cells", `11`="OPCs",
                                `12`="CHOR", `13`="vascular cells"
)

pdf("Dap12_integrated_umap_annotation.pdf", width=6, height=3.8)
DimPlot(Dap12_integrated, reduction = 'umap', label = F)
dev.off()

Dap12_integrated$celltype.orig.ident <- paste(Idents(Dap12_integrated), Dap12_integrated$orig.ident, sep = "_")
Dap12_integrated$celltype <- Idents(Dap12_integrated)


saveRDS(Dap12_integrated, file = "Dap12_integrated_Annotation.rds")

pdf("Dap12_QC.pdf", width=9, height=4)
Idents(Dap12_integrated) <- "Condition"
VlnPlot(object = Dap12_integrated, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, idents=NULL)
dev.off()


data <- Dap12_integrated
# calculate ratio of each genotype in each cell type cluster
a<-as.data.frame(table(data$Condition,data$celltype))
colnames(a)<-c("clusters","cell.type","cell.no")
agg<-aggregate(cell.no~clusters,a,sum)
a$cluster.total <- agg$cell.no[match(a$clusters,agg$clusters)]
a$ratio<-a$cell.no/a$cluster.total

ggplot(a,aes(x=clusters, y=ratio, fill=cell.type))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Genotype")+
  ylab("Cell type ratio per genotype") + RotatedAxis()

ggsave("genotype_celltype_distribution.pdf",plot=last_plot(),path="/athena/ganlab/.../integration",
       width=4,height=4,units="in")

#markers for annotation
pdf("annotation_1.pdf", width=10.5, height=2.7)
DotPlot(data, features = c("Plp1", "Mbp", "Mobp","Slc17a7", "Nrgn", "Clu", "Plpp3",
                           "Pla2g7", "Cx3cr1", "P2ry12", "Csf1r","Gad1", "Gad2","Vcan", "Pdgfra", "Bmp6", "Adam12",
                           "Cped1","Clic6","Ttr")) + RotatedAxis()
dev.off()


data <- Dap12_integrated
# calculate ratio of each sample in each cell type cluster
a<-as.data.frame(table(data$Sample_Name,data$celltype))
colnames(a)<-c("clusters","cell.type","cell.no")
agg<-aggregate(cell.no~clusters,a,sum)
a$cluster.total <- agg$cell.no[match(a$clusters,agg$clusters)]
a$ratio<-a$cell.no/a$cluster.total

ggplot(a,aes(x=clusters, y=ratio, fill=cell.type))+
  geom_bar(stat="identity")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  xlab("Sample")+
  ylab("Cell type ratio per sample") + RotatedAxis()

ggsave("sample_celltype_distribution.pdf",plot=last_plot(),path="/athena/ganlab/.../integration",
       width=6,height=4,units="in")

Idents(Dap12_integrated) <- "celltype"
DefaultAssay(Dap12_integrated) <- 'RNA'
pdf("Dap12_integrated_umap_annotation_noLabel.pdf", width=6, height=4)
DimPlot(Dap12_integrated, reduction = 'umap', label = F)
dev.off()

Cluster_EN <- subset(Dap12_integrated, idents = "excitatory neurons")
Cluster_IN <- subset(Dap12_integrated, idents = "inhibitory neurons")
Cluster_MG <- subset(Dap12_integrated, idents = "microglia")
Cluster_AST <- subset(Dap12_integrated, idents = "astrocytes")
Cluster_OL <- subset(Dap12_integrated, idents = "oligodendrocytes")
Cluster_OPC <- subset(Dap12_integrated, idents = "OPCs")
Cluster_VC <- subset(Dap12_integrated, idents = "vascular cells")
Cluster_CHOR <- subset(Dap12_integrated, idents = "CHOR")

saveRDS(Cluster_EN, file = "Dap12_EN_subset.rds")
saveRDS(Cluster_IN, file = "Dap12_IN_subset.rds")
saveRDS(Cluster_MG, file = "Dap12_MG_subset.rds")
saveRDS(Cluster_AST, file = "Dap12_AST_subset.rds")
saveRDS(Cluster_OL, file = "Dap12_OL_subset.rds")
saveRDS(Cluster_OPC, file = "Dap12_OPC_subset.rds")
saveRDS(Cluster_VC, file = "Dap12_VC_subset.rds")
saveRDS(Cluster_CHOR, file = "Dap12_CHOR_subset.rds")

