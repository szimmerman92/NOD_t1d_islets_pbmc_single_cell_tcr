#Load in the dataset
require(Seurat)
require(ggplot2)
setwd()
obj <- readRDS(file = "./Integrated_allmergedsamples.rds")
obj
slotNames(obj)
head(obj$nCount_RNA)
head(obj$TCR)

#Import marker list
imm_mark <- read.csv(file = "./ImmunemarkersSheet2.csv")

#Weighted PCA based on marker importance
x10 <- imm_mark$Gene[c(rep(1:3))]
x5 <- imm_mark$Gene[!imm_mark$Gene %in% x10]
fnames <- VariableFeatures(obj)[1:200]
vec <- c(rep(x10,each = 10),rep(imm_mark$Gene[!imm_mark$Gene %in% c(x10)],each = 5),fnames[!fnames %in% imm_mark$Gene])\
cdpr <- prcomp(x = GetAssayData(obj)[vec,])
obj[["pca"]] <- CreateDimReducObject(embeddings = cdpr$rotation,stdev = cdpr$sdev,key = "PCA_", assay = DefaultAssay(obj))

#Run UMAP
obj <- RunUMAP(obj,reduction = "pca",dims = 1:10)

#Run Clustering with a resolution of 1.5
obj <- FindNeighbors(obj,reduction = "pca",dims = c(1:10))
obj <- FindClusters(obj,resolution = 1.5)

#Save Idents and plot clusters
Idents(obj) <- "RNA_snn_res.1.5"
DimPlot(obj,reduction = "umap",group.by = "RNA_snn_res.1.5",label = T,raster = F)

#Find Markers
df <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,features = imm_mark$Gene)


#Add Names
popnames1.5 <- c("CD4_Naive_Il7r+",#"CD4_Naive1",
                 "CD8_Naive_Il7r+",#"CD8_Tcm1",
                 "DN_Naive_Il7r",#"DN_Tcm", #cluster2
                 "CD4_Tem_Il7r+",#CD4_Tem",
                 "DN_Tem_Il7r+",#"DN_Tem/Trm1",
                 "CD8_Tem_Il7r-",#CD8_Tscm", #cluster5
                 "CD4_Effector_Treg",#CD4_Tem/Teff",
                 "DN_Effector_Treg",#DN_Tem/Trm2",
                 "CD4_Tcm_Il7r-",#CD4_Tcm", #cluster8
                 "CD4_Tscm",
                 "CD4_Tscm_Il7r-",#CD4_Naive2",
                 "CD4_Th1_Il7r-",#CD4_Treg1",#cluster11
                 "CD8_Naive_Il7r+",#CD8_Tcm1",
                 "DN_Tcm_Il7r-",#DN_Tcm1",
                 "CD8_Tc1_Gzma-",#CD8_Tem",#cluster14
                 "CD4_Tem_Il7rlo",#CD4_Tscm5",
                 "DN_Tscm_Il7r-",#DN_Naive1",
                 "CD8_Tcm_Il7r+",#CD8_Tcm2",
                 "DN_Tcm_Il7r+",#DN_Treg/Tcm",#cluster18
                 "CD4_RTE",#CD4_Naive3",
                 "CD8_RTE",#probably not change later
                 "DN_pre_iTreg_Il7r-",#DN_Treg1",
                 "CD4_Th1_Il7r+",#CD4_Treg2",#cluster22
                 "DN_Th1_Il7r+",#DN_Activated",
                 "DN_pre_nTreg_Il7r-",#DN_Tcm2",
                 "CD8_Trm/Tc1_Gzma+",#CD8_Tem/Tc1",
                 "DN_Naive_Il7r-",#DN_Naive2",#cluster26
                 "CD8_Trm/Tc1_Gzma+",#"CD8_Tem/Tc1w32",
                 "CD8_Trm/Tc1_Gzma+",#"CD8_Tem/Tc134",
                 "CD8_Naive_Il7r-",#CD8_Naive",
                 "CD4_preTreg_Il7r-",#CD4_Treg/Tcm",
                 "DP_naive",#cluster31
                 "DN_Tscm",
                 "DN_Treg",
                 "DN_Trm/Tc1_Il7r-")#DN_Tem")
names(popnames1.5) <- levels(obj)
obj <- RenameIdents(obj, popnames1.5)


