#set working directory ====
setwd("/athena/ganlab/.../DF_1stRound")
#install.packages
library(Seurat)
library(ggplot2)
library(DoubletFinder)
library(SoupX)

#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/.../cellranger/B70/outs')
sc = autoEstCont(sc)
WT_1.counts = adjustCounts(sc)
WT_1 <- CreateSeuratObject(counts = WT_1.counts, project = "Dap12_B70", min.cells = 3, min.features = 200)
WT_1[["Condition"]] = c('WT')
WT_1[["Sample_Name"]] = c('WT_1')
rm(WT_1.counts)
#vizualize QC metrics and filtering====
WT_1[["percent.mt"]] <- PercentageFeatureSet(object = WT_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- WT_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("WT_1_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("WT_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("WT_1_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
saveRDS(all,"WT_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("WT_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/.../cellranger/B71/outs')
sc = autoEstCont(sc)
WT_2.counts = adjustCounts(sc)
WT_2 <- CreateSeuratObject(counts = WT_2.counts, project = "Dap12_B71", min.cells = 3, min.features = 200)
WT_2[["Condition"]] = c('WT')
WT_2[["Sample_Name"]] = c('WT_2')
rm(WT_2.counts)
#vizualize QC metrics and filtering====
WT_2[["percent.mt"]] <- PercentageFeatureSet(object = WT_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- WT_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("WT_2_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("WT_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("WT_2_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
saveRDS(all,"WT_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("WT_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/.../cellranger/LG180_503/outs')
sc = autoEstCont(sc)
WT_3.counts = adjustCounts(sc)
WT_3 <- CreateSeuratObject(counts = WT_3.counts, project = "LG180_503", min.cells = 3, min.features = 200)
WT_3[["Condition"]] = c('WT')
WT_3[["Sample_Name"]] = c('WT_3')
rm(WT_3.counts)
#vizualize QC metrics and filtering====
WT_3[["percent.mt"]] <- PercentageFeatureSet(object = WT_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- WT_3
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("WT_3_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("WT_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("WT_3_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
saveRDS(all,"WT_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("WT_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/.../cellranger/LG321_NY_53/outs')
sc = autoEstCont(sc)
DAP12KO_1.counts = adjustCounts(sc)
DAP12KO_1 <- CreateSeuratObject(counts = DAP12KO_1.counts, project = "Dap12_53", min.cells = 3, min.features = 200)
DAP12KO_1[["Condition"]] = c('DAP12KO')
DAP12KO_1[["Sample_Name"]] = c('DAP12KO_1')
rm(DAP12KO_1.counts)
#vizualize QC metrics and filtering====
DAP12KO_1[["percent.mt"]] <- PercentageFeatureSet(object = DAP12KO_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- DAP12KO_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("DAP12KO_1_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("DAP12KO_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("DAP12KO_1_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
saveRDS(all,"DAP12KO_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("DAP12KO_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/.../cellranger/LG321_NY_54/outs')
sc = autoEstCont(sc)
DAP12KO_2.counts = adjustCounts(sc)
DAP12KO_2 <- CreateSeuratObject(counts = DAP12KO_2.counts, project = "Dap12_54", min.cells = 3, min.features = 200)
DAP12KO_2[["Condition"]] = c('DAP12KO')
DAP12KO_2[["Sample_Name"]] = c('DAP12KO_2')
rm(DAP12KO_2.counts)
#vizualize QC metrics and filtering====
DAP12KO_2[["percent.mt"]] <- PercentageFeatureSet(object = DAP12KO_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- DAP12KO_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("DAP12KO_2_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("DAP12KO_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("DAP12KO_2_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
saveRDS(all,"DAP12KO_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("DAP12KO_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/.../cellranger/LG321_59/outs')
sc = autoEstCont(sc)
DAP12KO_3.counts = adjustCounts(sc)
DAP12KO_3 <- CreateSeuratObject(counts = DAP12KO_3.counts, project = "Dap12_59", min.cells = 3, min.features = 200)
DAP12KO_3[["Condition"]] = c('DAP12KO')
DAP12KO_3[["Sample_Name"]] = c('DAP12KO_3')
rm(DAP12KO_3.counts)
#vizualize QC metrics and filtering====
DAP12KO_3[["percent.mt"]] <- PercentageFeatureSet(object = DAP12KO_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- DAP12KO_3
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("DAP12KO_3_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("DAP12KO_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("DAP12KO_3_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
saveRDS(all,"DAP12KO_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("DAP12KO_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/.../cellranger/LG321_NY_60/outs')
sc = autoEstCont(sc)
DAP12WT_1.counts = adjustCounts(sc)
DAP12WT_1 <- CreateSeuratObject(counts = DAP12WT_1.counts, project = "Dap12_60", min.cells = 3, min.features = 200)
DAP12WT_1[["Condition"]] = c('DAP12WT')
DAP12WT_1[["Sample_Name"]] = c('DAP12WT_1')
rm(DAP12WT_1.counts)
#vizualize QC metrics and filtering====
DAP12WT_1[["percent.mt"]] <- PercentageFeatureSet(object = DAP12WT_1, pattern = "^mt-") #recognize mitochondrial transcripts
all <- DAP12WT_1
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("DAP12WT_1_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("DAP12WT_1_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("DAP12WT_1_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
saveRDS(all,"DAP12WT_1_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("DAP12WT_1_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/.../cellranger/LG321_NY_65/outs')
sc = autoEstCont(sc)
DAP12WT_2.counts = adjustCounts(sc)
DAP12WT_2 <- CreateSeuratObject(counts = DAP12WT_2.counts, project = "Dap12_65", min.cells = 3, min.features = 200)
DAP12WT_2[["Condition"]] = c('DAP12WT')
DAP12WT_2[["Sample_Name"]] = c('DAP12WT_2')
rm(DAP12WT_2.counts)
#vizualize QC metrics and filtering====
DAP12WT_2[["percent.mt"]] <- PercentageFeatureSet(object = DAP12WT_2, pattern = "^mt-") #recognize mitochondrial transcripts
all <- DAP12WT_2
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("DAP12WT_2_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("DAP12WT_2_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("DAP12WT_2_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
saveRDS(all,"DAP12WT_2_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("DAP12WT_2_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################
#for loading Cell Ranger counts:
sc = load10X('/athena/ganlab/.../cellranger/LG321_NY_68/outs')
sc = autoEstCont(sc)
DAP12WT_3.counts = adjustCounts(sc)
DAP12WT_3 <- CreateSeuratObject(counts = DAP12WT_3.counts, project = "Dap12_68", min.cells = 3, min.features = 200)
DAP12WT_3[["Condition"]] = c('DAP12WT')
DAP12WT_3[["Sample_Name"]] = c('DAP12WT_3')
rm(DAP12WT_3.counts)
#vizualize QC metrics and filtering====
DAP12WT_3[["percent.mt"]] <- PercentageFeatureSet(object = DAP12WT_3, pattern = "^mt-") #recognize mitochondrial transcripts
all <- DAP12WT_3
#plot of correlation between number of genes detected and number of transcripts detected - you generally want this to be 1:1
pdf("DAP12WT_3_FeatureScatter.pdf", width=12, height=4)
plot1 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(all, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()
#initial filtering step - usually you want cells with detected gene counts over 300 and mitochondrial transcripts below 5%
all <- subset(x = all, subset = nFeature_RNA > 300 & nFeature_RNA < 4000 & nCount_RNA < 20000 & percent.mt < 5)
#normalize counts=====
all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
all <- ScaleData(object = all)
all <- RunPCA(object = all, features = rownames(x = all), verbose = FALSE)
pdf("DAP12WT_3_Elbow_1.pdf", width=8, height=6)
ElbowPlot(all)
dev.off()
all <- FindNeighbors(all, dims = 1:15)
all <- FindClusters(all, resolution = 0.1)
all <- RunUMAP(all, dims = 1:15)
pdf("DAP12WT_3_UMAP.pdf", width=8, height=6)
DimPlot(all, reduction = "umap", label = T)
dev.off()
saveRDS(all,"DAP12WT_3_QC.rds")
#Doublet finder (no ground-truth) - please reference https://github.com/chris-mcginnis-ucsf/DoubletFinder for more information on parameters====
sweep.pbmc <- paramSweep_v3(all,PCs=1:15,sct=FALSE)
sweep.stats_pbmc <- summarizeSweep(sweep.pbmc,GT=FALSE)
bcmvn_pbmc <- find.pK(sweep.stats_pbmc)
pdf("DAP12WT_3_ggplot_pK.pdf", width=18, height=6)
ggplot(bcmvn_pbmc, aes(x=bcmvn_pbmc$pK, y=bcmvn_pbmc$BCmetric))+geom_bar(stat="identity") #look for pK at the initial peak
dev.off()
###############################################################################################
###############################################################################################




