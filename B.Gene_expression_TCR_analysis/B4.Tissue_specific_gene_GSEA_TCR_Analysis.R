#!/usr/bin/env Rscript
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(ggplot2)
library(cowplot)
library(viridis)
library(gridExtra)
library(RColorBrewer)
library(tibble)


############ A. Analysis for T1D PBMC ###########

#### A1. load integrated rds ####
Integratedt1dpbmc <- readRDS("/Users/zahir/Documents/Tcell/seurat/1.DataProcessing/Integrated_rds/Integratedt1dpbmc.rds")

#For performing differential expression after integration, we switch back to the original data
DefaultAssay(Integratedt1dpbmc) <- "RNA"

#### A2. #Make UMAP plot ####
pdf(file = "t1dpbmc_integrated_unlabeled_umap.pdf", width = 8, height = 6)
pa1 <- DimPlot(Integratedt1dpbmc, reduction = "umap", label = TRUE, repel = TRUE)
pa1
dev.off()


#### A3. Cluster cell population identification ####

#First we will use Findconserved markers function to identified conserved markers across different clusters. 
#using this loop for all 0-10 clusters
for (i in 0:10){
  marker_i <- FindConservedMarkers(Integratedt1dpbmc, ident.1 = i, grouping.var = "sampletype", verbose =TRUE)
  filename <- paste0("t1dpbmc_cluster.", i,".csv")
  write.csv(marker_i, filename)
}

#Now we will identify cluster specific marker by differential expression of markers in a specific cluster
t1dpbmc_findall_markers <-FindAllMarkers(Integratedt1dpbmc, min.pct = 0.25, logfc.threshold = 0.25)

write.csv(t1dpbmc_findall_markers, "T1D_PBMC_Find_cluster_markers/t1dpbmc_findall_markers.csv", quote = F)

#if there is any confusion about expression of any canonical marker genes then try to solve it by violin ploting or featureploting 

#plot by simple violin 
plots <- VlnPlot(Integratedt1dpbmc, features = c("Tmem176b", "Tmem176a", "Il12rb1", "Cd4"), group.by = "seurat_clusters",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)




#FeaturePlot(object = Integratedt1dpbmc, features = "Cd4")

#staked violin plot in SeuratTheme
## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj,
                          feature,
                          pt.size = 0,
                          plot.margin = margin(0, 0, 0, 0, "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, group.by = "seurat_clusters", pt.size = pt.size, ... )+
    ylab(feature) +
    theme(legend.position = "none",
          plot.title= element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = rel(1), angle = 0),
          axis.text.y = element_text(size = rel(1)),
          plot.margin = plot.margin )
  return(p)
}

# extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = margin(0, 0, 0, 0, "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5), axis.ticks.x = element_line(), legend.position = "right")
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


features<- c("Ccl5", "Ccr10", "Ccr6","Ccr7", "Cd28", "Cd38", "Cd4", "Cd44", "Cd69", "Cd7", "Cd74", "Cd8a", "Cd8b1", "Ctla2a", "Ctla4", "Cxcr5", "Fasl", "Foxp3", "Gata3", "Ifng", "Igfbp4", "Il12rb1", "Il17a", "Il1r1", "Il2ra", "Il2rb", "Il7r", "Itgb1", "Klrd1", "Lef1", "Ly6a", "Nkg7", "Rora", "Sell", "Tbx21", "Tmem176a", "Tmem176b")

pdf(file = "t1dpbmc_canonical_violin_staked.pdf", width = 12, height = 16)
stackedplot_violin <- StackedVlnPlot(obj = Integratedt1dpbmc, features = features)
stackedplot_violin
dev.off()

#rename cluster ID 
DefaultAssay(Integratedt1dpbmc) <- "RNA"

Integratedt1dpbmc_labelled <- RenameIdents(Integratedt1dpbmc, `0` = "CD4 naive T cells", `1` = "CD8 central memory CTLs", `2` = "CD4 central memory T cells", `3` = "CD4 effector memory T cells", `4` = "CD4 CD25 Foxp3 Tregs", `5` = "CD8 effector memory CTLs", `6` = "CD8 CTLA2a Trm CTLs", `7` = "DN CD7 NKT cells", `8` = "CD4 CD38 CD74 T cells", `9` = "CD4 CD74 T cells",`10` = "DN CD25 Th1/Th17 like cells")

t1d_pbmcs_colors = c('CD4 naive T cells'='#F0A0FF','CD8 central memory CTLs'='#0075DC','CD4 central memory T cells'='#993F00','CD4 CD25 Foxp3 Tregs'='#191919','CD8 effector memory CTLs'='#005C31','CD8 CTLA2a Trm CTLs'='#2BCE48','DN CD7 NKT cells'='#FFCC99','DN CD25 Th1/Th17 like cells'='#8F7C00','CD4 effector memory T cells'='#4C005C','CD4 CD38 CD74 T cells'='#808080','CD4 CD74 T cells'='#94FFB5')


#colored umap
pdf(file = "t1dpbmc_labelled_umap_colored.pdf", width = 9, height = 7)
pa2 <- DimPlot(Integratedt1dpbmc_labelled,cols=t1d_pbmcs_colors, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 0.2, label.size = 5) + theme(axis.text = element_text(size = 20)) + theme(legend.text=element_text(size=16))
pa2
dev.off()



# Supplementary Dotplot for canonical markers for supplementary with UMAP 
Idents(Integratedt1dpbmc_labelled) <- factor(Idents(Integratedt1dpbmc_labelled), levels = c("CD4 naive T cells", "CD8 central memory CTLs", "CD4 central memory T cells", "CD4 effector memory T cells", "CD4 CD25 Foxp3 Tregs", "CD8 effector memory CTLs", "CD8 CTLA2a Trm CTLs", "DN CD7 NKT cells", "CD4 CD38 CD74 T cells", "CD4 CD74 T cells", "DN CD25 Th1/Th17 like cells"))


canonicalmarker <- c("Ccl5", "Ccr10", "Ccr6","Ccr7", "Cd28", "Cd38", "Cd4", "Cd44", "Cd69", "Cd7", "Cd74", "Cd8a", "Cd8b1", "Ctla2a", "Ctla4", "Cxcr5", "Fasl", "Foxp3", "Gata3", "Ifng", "Igfbp4", "Il12rb1", "Il17a", "Il1r1", "Il2ra", "Il2rb", "Il7r", "Itgb1", "Klrd1", "Lef1", "Ly6a", "Nkg7", "Rora", "Sell", "Tbx21", "Tmem176a", "Tmem176b")


pdf(file = "t1dpbmc_canonicaldotplot.pdf", width = 12, height = 6)
canonicaldotplot <- DotPlot(Integratedt1dpbmc_labelled, features = canonicalmarker, cols = c("blue", "red", "green")) + RotatedAxis()
canonicaldotplot
dev.off()


#Cell count in different cluster
cell.num <- table(Idents(Integratedt1dpbmc_labelled))

write.csv(cell.num, "t1dpbmc_cellcountbycluster.csv", quote = F)



#### A4. T cell clone matching ####
#plot islets matching cells in blood

pdf(file = "t1dpbmc_isletsmatching_umap.pdf", width = 8, height = 6)
t1dpbmc_isletsmatching_umap <- DimPlot(object = Integratedt1dpbmc_labelled, sizes.highlight=0.08, cells.highlight = WhichCells(object = Integratedt1dpbmc_labelled, expression = Matching_pre_filter == "matching")) + scale_color_manual(labels = c("Non-matching", "Islets matching"), values = c("grey", "red"))+ theme(axis.text = element_text(size = 20)) + theme(legend.text=element_text(size=16))
t1dpbmc_isletsmatching_umap
dev.off()

                                                                                                                                                                                                                                                                                

pdf(file = "t1dpbmc_isletsmatching_umap_bysampleid.pdf", width = 16, height = 4)
t1dpbmc_isletsmatching_umap_bysampleid <- DimPlot(object = Integratedt1dpbmc_labelled, cells.highlight = WhichCells(object = Integratedt1dpbmc_labelled, expression = Matching_pre_filter == "matching"), split.by = "sampleid") + scale_color_manual(labels = c("Non-matching", "Islets matching"), values = c("grey", "red"))
t1dpbmc_isletsmatching_umap_bysampleid
dev.off()


#Plot clonal expansion
pdf(file = "t1dpbmc_isletsmatching_clinalexpansion.pdf", width = 8, height = 6)
t1dpbmc_isletsmatching_clinalexpansion <- FeaturePlot(object = Integratedt1dpbmc_labelled, features = "Clone.size", pt.size = 0.08)+ theme(axis.text = element_text(size = 20)) + theme(legend.text=element_text(size=16))
t1dpbmc_isletsmatching_clinalexpansion
dev.off()


####plot and significance test of clonl expansion###
Integratedt1dpbmc_metadata= Integratedt1dpbmc@meta.data
Integratedt1dpbmc_metadata$Matching[grep("_pbmc_matching",Integratedt1dpbmc_metadata$Matching)] = 'Islets-matching'
Integratedt1dpbmc_metadata$Matching[Integratedt1dpbmc_metadata$Matching == "alpha_chain_only"] = "Not-matching"
Integratedt1dpbmc_metadata$Matching[Integratedt1dpbmc_metadata$Matching == "beta_chain_only"] = "Not-matching"
Integratedt1dpbmc_metadata$Matching[Integratedt1dpbmc_metadata$Matching == "notcr"] = "Not-matching"
Integratedt1dpbmc_metadata$Matching[Integratedt1dpbmc_metadata$Matching == "not_matching"] = "Not-matching"

Integratedt1dpbmc_metadata$Clone.size = log10(Integratedt1dpbmc_metadata$Clone.size+1)

# do wicoin rank sum test/Mann-whitney test
clone_size_matching = Integratedt1dpbmc_metadata[Integratedt1dpbmc_metadata$Matching == "Islets-matching","Clone.size"]
clone_size_non_matching = Integratedt1dpbmc_metadata[Integratedt1dpbmc_metadata$Matching == "Not-matching","Clone.size"]
wilcox_res_matching_vs_nonmatching = wilcox.test(clone_size_matching,clone_size_non_matching,alternative="greater")
wilcox_res_matching_vs_nonmatching


# make boxplot
pdf(file = "t1dpbmc_isletsmatching_clinalexpansion_significance.pdf", width = 4, height = 4)
p1 <- ggplot(Integratedt1dpbmc_metadata,aes(x=Matching,y=Clone.size, fill = Matching)) + geom_boxplot () + scale_fill_manual(values=c("red","gray"))+labs(y=expression("Log"[10]*"(Clone size)"))+ 
  theme(axis.title.x=element_blank(),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border= element_blank(),
        legend.position = "none",
        text = element_text(size = 16))+
  annotate(geom="text", x=1.7, y=1.7, label="Wilcoxon, p-value < 2.2e-16", color="black", size = 4)
p1
dev.off()



# cell count by group
library(data.table)
library(magrittr)
# extract meta data
metadata <- Integratedt1dpbmc_labelled@meta.data %>% as.data.table
# count the number of cells per unique combinations of "matcing clone number" and "seurat_clusters"

cell.num.by.matchingclone.seuratcluster <- metadata[, .N, by = c("Matching_pre_filter", "seurat_clusters")] %>% dcast(., Matching_pre_filter ~ seurat_clusters, value.var = "N")
write.csv(cell.num.by.matchingclone.seuratcluster, "t1dpbmc_cell.num.by.matchingclone.seuratcluster.csv", quote = F)


#count by clone size
cell.num.by.clonesize <- metadata[, .N, by = c("Matching_pre_filter", "Clone.size")] %>% dcast(., Matching_pre_filter ~ Clone.size, value.var = "N")
write.csv(cell.num.by.clonesize, "t1dpbmc_cell.num.by.clonesize.csv", quote = F)


#count by sample id
cell.num.by.sample <- metadata[, .N, by = c("sampleid", "seurat_clusters")] %>% dcast(., sampleid ~ seurat_clusters, value.var = "N")
write.csv(cell.num.by.sample, "t1dpbmc_cell.num.by.sample.csv", quote = F)


#### A5. DE islets matching T cells and non-matching T cells ####

#compare all cluster together
allcluster_matching_nonmatching <- FindMarkers(Integratedt1dpbmc_labelled, ident.1 = "matching", ident.2 = "not_matching", group.by = 'Matching_pre_filter')

#to write the results in csv
write.csv(allcluster_matching_nonmatching, file ="DE_t1dpbmc_allcluster_matching_nonmatching.csv")

#We can also perform DE in a specific cluster
cluster3_CD4_effector_memory_matching_nonmatching <- FindMarkers(Integratedt1dpbmc_labelled, ident.1 = "matching", ident.2 = "not_matching", group.by = 'Matching_pre_filter', subset.ident = "CD4 effector memory T cells")

#write.csv(cluster3_CD4_effector_memory_matching_nonmatching, file ="DE_t1dpbmc_CD4)_effector_memory_T_cells_matching_nonmatching.csv")


cluster5_CD8_effector_memory_CTLs_matching_nonmatching <- FindMarkers(Integratedt1dpbmc_labelled, ident.1 = "matching", ident.2 = "not_matching", group.by = 'Matching_pre_filter', subset.ident = "CD8 effector memory CTLs")

#write.csv(cluster5_CD8_effector_memory_CTLs_matching_nonmatching, file ="DE_t1dpbmc_CD8_effector_memory_T_cells_matching_nonmatching.csv")


####making volkano plot with the differentially expressed genes
#https://bioconductor.org/packages/release/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
library(EnhancedVolcano)

pdf(file = "volcano_DE_islets-matching_vs_non-mathcing_in_PBMC.pdf", width = 12, height = 10)
pbmc_volcano <- EnhancedVolcano(allcluster_matching_nonmatching,
                lab = rownames(allcluster_matching_nonmatching),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'DE islets-matching vs non-mathcing in PBMC',
                pCutoff = 10e-10,
                FCcutoff = 1.5,
                pointSize = 3.0,
                labSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5)+theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border= element_blank())

pbmc_volcano
dev.off()


pdf(file = "DE_islets-matching_vs_non-mathcing_in_t1dPBMC_CD4_effector_memory.pdf", width = 8, height = 8)
pbmc_cluster3_CD4_effector_memory_volcano <- EnhancedVolcano(cluster3_CD4_effector_memory_matching_nonmatching,
                                lab = rownames(cluster3_CD4_effector_memory_matching_nonmatching),
                                x = 'avg_log2FC',
                                y = 'p_val_adj',
                                title = 'volcano_DE_islets-matching_vs_non-mathcing_in_t1dPBMC_cluster3_CD4_effector_memory',
                                pCutoff = 10e-10,
                                FCcutoff = .5,
                                pointSize = 3.0,
                                labSize = 5.0,
                                drawConnectors = TRUE,
                                widthConnectors = 0.5)+theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border= element_blank())

pbmc_cluster3_CD4_effector_memory_volcano
dev.off()


pdf(file = "DE_islets-matching_vs_non-mathcing_t1dPBMC_CD8_effector_memory.pdf", width = 6, height = 6)
cluster5_CD8_effector_memory_CTLs_volcano <- EnhancedVolcano(cluster5_CD8_effector_memory_CTLs_matching_nonmatching,
                                                             lab = rownames(cluster5_CD8_effector_memory_CTLs_matching_nonmatching),
                                                             x = 'avg_log2FC',
                                                             y = 'p_val_adj',
                                                             title = 'volcano_DE_islets-matching_vs_non-mathcing_in_t1dPBMC_cluster5_CD8_effector_memory_CTLs',
                                                             pCutoff = 10e-3,
                                                             FCcutoff = .5,
                                                             pointSize = 3.0,
                                                             labSize = 5.0,
                                                             drawConnectors = TRUE,
                                                             widthConnectors = 0.5)+theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border= element_blank())

cluster5_CD8_effector_memory_CTLs_volcano
dev.off()




#plot Il17a expression in all population of t1dpbmc, t1dislets and nont1dpbmc
#plot Il17a in all sample of cluster 10 in t1d pbmc 
#the follwing command needs to run the stacked violin loop function
features_Il17a<- c("Il17a")

pdf(file = "t1dpbmc_Il17a_violin.pdf", width = 4, height = 3)
stackedplot_violin_il17a <- StackedVlnPlot(obj = Integratedt1dpbmc, features = features_Il17a)
stackedplot_violin_il17a
dev.off()

pdf(file = "t1dislets_Il17a_violin.pdf", width = 4, height = 3)
stackedplot_violin_il17a <- StackedVlnPlot(obj = Integratedt1dislets, features = features_Il17a)
stackedplot_violin_il17a
dev.off()

pdf(file = "nont1dpbmc_Il17a_violin.pdf", width = 4, height = 3)
stackedplot_violin_il17a <- StackedVlnPlot(obj = Integratednont1dpbmc, features = features_Il17a)
stackedplot_violin_il17a
dev.off()

#by sampleid
pdf(file = "t1dpbmc_Il17a_violin_bysampleid.pdf", width = 4, height = 3)
v1 <- VlnPlot(subset(Integratedt1dpbmc_labelled, idents = "DN CD25 Th1/Th17 like cells"), features = "Il17a", group.by = "sampleid", cols=c("brown4", "brown4", "brown4", "brown4", "brown4")) + theme(legend.position = 'none')
v1
dev.off()


pdf(file = "t1dpbmc_Il17a_violin_bysampleid_matchingclone.pdf", width = 4, height = 3)
v2 <- VlnPlot(subset(Integratedt1dpbmc_labelled, idents = "DN CD25 Th1/Th17 like cells"), features = "Il17a", group.by = "Matching_pre_filter")
v2
dev.off()


pdf(file = "t1dpbmc_vs_nont1dpbmc_Il17a_violin.pdf", width = 4, height = 3)
v3 <- VlnPlot(Integrated_allmergedsamples, features = "Il17a", group.by = "sampletype")
v3
dev.off()



VlnPlot(Integratedt1dpbmc_labelled, features = "Il17a")






############ B. Analysis for T1D islets ###########

#### B1. load integrated rds ####
Integratedt1dislets <- readRDS("/Users/zahir/Documents/Tcell/seurat/1.DataProcessing/Integrated_rds/Integratedislets.rds")

#For performing differential expression after integration, we switch back to the original data
DefaultAssay(Integratedt1dislets) <- "RNA"

#### B2. #Make UMAP plot ####
pdf(file = "t1dislets_integrated_unlabeled_umap.pdf", width = 8, height = 6)
pb1 <- DimPlot(Integratedt1dislets, reduction = "umap", label = TRUE, repel = TRUE)
pb1
dev.off()


#### B3. Cluster cell population identification ####

#First we will use Findconserved markers function to identified conserved markers across different clusters. 
#using this loop for all 0-10 clusters
for (i in 0:16){
  marker_i <- FindConservedMarkers(Integratedt1dislets, ident.1 = i, grouping.var = "sampletype", verbose =TRUE)
  filename <- paste0("t1dislets_cluster.", i,".csv")
  write.csv(marker_i, filename)
}

#Now we will identify cluster specific marker by differential expression of markers in a specific cluster
t1dislets_findall_markers <-FindAllMarkers(Integratedt1dislets, min.pct = 0.25, logfc.threshold = 0.25)

write.csv(t1dislets_findall_markers, "t1d_islets_Find_cluster_markers/t1dislets_findall_markers.csv", quote = F)

#if there is any confusion about expression of any canonical marker genes then try to solve it by violin ploting or featureploting 

#plot by simple violin 
plots <- VlnPlot(Integratedt1dislets, features = c("Tmem176b", "Tmem176a", "Il12rb1", "Cd4"), group.by = "seurat_clusters",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

#FeaturePlot(object = Integratedt1dpbmc, features = "Cd4")

#staked violin plot in SeuratTheme
## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj,
                          feature,
                          pt.size = 0,
                          plot.margin = margin(0, 0, 0, 0, "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, group.by = "seurat_clusters", pt.size = pt.size, ... )+
    ylab(feature) +
    theme(legend.position = "none",
          plot.title= element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = rel(1), angle = 0),
          axis.text.y = element_text(size = rel(1)),
          plot.margin = plot.margin )
  return(p)
}

# extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = margin(0, 0, 0, 0, "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5), axis.ticks.x = element_line(), legend.position = "right")
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


features<- c("Ccl5", "Ccr10", "Ccr6","Ccr7", "Cd28", "Cd38", "Cd4", "Cd44", "Cd69", "Cd7", "Cd74", "Cd8a", "Cd8b1", "Ctla2a", "Ctla4", "Cxcr5", "Fasl", "Foxp3", "Gata3", "Ifng", "Igfbp4", "Il12rb1", "Il17a", "Il1r1", "Il2ra", "Il2rb", "Il7r", "Itgb1", "Klrd1", "Lef1", "Ly6a", "Nkg7", "Rora", "Sell", "Tbx21", "Tmem176a", "Tmem176b")

pdf(file = "t1dislets_canonical_violin_staked.pdf", width = 12, height = 16)
stackedplot_violin_islets <- StackedVlnPlot(obj = Integratedt1dislets, features = features)
stackedplot_violin_islets
dev.off()


#rename cluster ID 
DefaultAssay(Integratedt1dislets) <- "RNA"

Integratedt1dislets_labelled <- RenameIdents(Integratedt1dislets, `0` = "CD4 Lef1 memory T cells", `1` = "CD4 CD29 CD25 T cells", `2` = "CD8 Ly6a central memory CTLs", `3` = "CD4 central memory T cells", `4` = "CD4 CD25 CTLA4 Foxp3- Tregs", `5` = "CD8 effector memory CTLs", `6` = "CD4 CD25 Th1 cells", `7` = "CD4 CD25 CTLA4 Foxp3 Tregs", `8` = "CD8 central memory CTLs", `9` = "CD4 CTLA4 memory T cells",`10` = "DN CD7 NKT cells", `11` = "CD8 CTLA2a Trm CTLs", `12` = "DN CD74 memory T cells", `13` = "CD4 Th17 like cells", `14` = "DP CD25 CTLA4 T cells", `15` = "CD4 Th17 memory like cells", `16` = "DN CTLA4 Th17 like cells")


t1d_islets_colors =c('CD4 Lef1 memory T cells'='#F0A0FF','CD4 CD29 CD25 T cells'='#854c65','CD8 Ly6a central memory CTLs'= '#00FF00','CD4 central memory T cells'='#993F00','CD4 CD25 CTLA4 Foxp3- Tregs'='#BDB76B','CD8 effector memory CTLs'='#005C31','CD4 CD25 Th1 cells'='#a19d94','CD4 CD25 CTLA4 Foxp3 Tregs'='#191919','CD8 central memory CTLs'='#0075DC','CD4 CTLA4 memory T cells'='#FF0000','DN CD7 NKT cells'='#87CEEB','CD8 CTLA2a Trm CTLs'='#FFA500','DN CD74 memory T cells'='#E0FF66','CD4 Th17 like cells'='violet','DP CD25 CTLA4 T cells'='#722F37','CD4 Th17 memory like cells'='#FFFF00','DN CTLA4 Th17 like cells'= '#ffa010')

#library("pals")
#pal.bands(alphabet, main="Discrete", show.names=TRUE)


#colored umap
pdf(file = "t1dislets_labelled_umap_colored.pdf", width = 10, height = 6)
pb2 <- DimPlot(Integratedt1dislets_labelled,cols=t1d_islets_colors, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 0.2, label.size = 5) + theme(axis.text = element_text(size = 20)) + theme(legend.text=element_text(size=14))
pb2
dev.off()



# Supplementary Dotplot for canonical markers for supplementary with UMAP 
Idents(Integratedt1dislets_labelled) <- factor(Idents(Integratedt1dislets_labelled), levels = c("CD4 Lef1 memory T cells", "CD4 CD29 CD25 T cells", "CD8 Ly6a central memory CTLs", "CD4 central memory T cells", "CD4 CD25 CTLA4 Foxp3- Tregs", "CD8 effector memory CTLs", "CD4 CD25 Th1 cells", "CD4 CD25 CTLA4 Foxp3 Tregs", "CD8 central memory CTLs", "CD4 CTLA4 memory T cells", "DN CD7 NKT cells", "CD8 CTLA2a Trm CTLs", "DN CD74 memory T cells", "CD4 Th17 like cells", "DP CD25 CTLA4 T cells", "CD4 Th17 memory like cells", "DN CTLA4 Th17 like cells"))


canonicalmarker <- c("Ccl5", "Ccr10", "Ccr6","Ccr7", "Cd28", "Cd38", "Cd4", "Cd44", "Cd69", "Cd7", "Cd74", "Cd8a", "Cd8b1", "Ctla2a", "Ctla4", "Cxcr5", "Fasl", "Foxp3", "Gata3", "Ifng", "Igfbp4", "Il12rb1", "Il17a", "Il1r1", "Il2ra", "Il2rb", "Il7r", "Itgb1", "Klrd1", "Lef1", "Ly6a", "Nkg7", "Rora", "Sell", "Tbx21", "Tmem176a", "Tmem176b")


pdf(file = "t1dislets_canonicaldotplot.pdf", width = 12, height = 6)
canonicaldotplot_islets <- DotPlot(Integratedt1dislets_labelled, features = canonicalmarker, cols = c("blue", "red", "green")) + RotatedAxis()
canonicaldotplot_islets
dev.off()


#Cell count in different cluster
cell.num <- table(Idents(Integratedt1dislets_labelled))

write.csv(cell.num, "t1dislets_cellcountbycluster.csv", quote = F)



#### B4. T cell clone matching ####
#plot PBMC matching cells in islets

pdf(file = "t1dislets_pbmcmatching_umap.pdf", width = 8, height = 6)
t1dislets_pbmcmatching_umap <- DimPlot(object = Integratedt1dislets_labelled, sizes.highlight=0.08, cells.highlight = WhichCells(object = Integratedt1dislets_labelled, expression = Matching_pre_filter == "matching")) + scale_color_manual(labels = c("Non-matching", "PBMC matching"), values = c("grey", "red"))+ theme(axis.text = element_text(size = 20)) + theme(legend.text=element_text(size=16))
t1dislets_pbmcmatching_umap
dev.off()



pdf(file = "t1dislets_pbmcmatching_umap_bysampleid.pdf", width = 16, height = 4)
t1dislets_pbmcmatching_umap_bysampleid<- DimPlot(object = Integratedt1dislets_labelled, cells.highlight = WhichCells(object = Integratedt1dislets_labelled, expression = Matching_pre_filter == "matching"), split.by = "sampleid") + scale_color_manual(labels = c("Non-matching", "PBMC matching"), values = c("grey", "red"))
t1dislets_pbmcmatching_umap_bysampleid
dev.off()


#Plot clonal expansion
pdf(file = "t1dislets_pbmcmatching_clinalexpansion.pdf", width = 7, height = 6)
t1dislets_pbmcmatching_clinalexpansion <- FeaturePlot(object = Integratedt1dislets_labelled, features = "Clone.size", pt.size = 0.08)+ theme(axis.text = element_text(size = 20)) + theme(legend.text=element_text(size=16))
t1dislets_pbmcmatching_clinalexpansion
dev.off()



####plot and significance test of clonl expansion###
Integratedt1dislets_metadata= Integratedt1dislets@meta.data
Integratedt1dislets_metadata$Matching[grep("_islets_matching",Integratedt1dislets_metadata$Matching)] = 'PBMC-matching'
Integratedt1dislets_metadata$Matching[Integratedt1dislets_metadata$Matching == "alpha_chain_only"] = "Not-matching"
Integratedt1dislets_metadata$Matching[Integratedt1dislets_metadata$Matching == "beta_chain_only"] = "Not-matching"
Integratedt1dislets_metadata$Matching[Integratedt1dislets_metadata$Matching == "notcr"] = "Not-matching"
Integratedt1dislets_metadata$Matching[Integratedt1dislets_metadata$Matching == "not_matching"] = "Not-matching"

Integratedt1dislets_metadata$Clone.size = log10(Integratedt1dislets_metadata$Clone.size+1)

# do wicoin rank sum test/Mann-whitney test
islets_clone_size_matching = Integratedt1dislets_metadata[Integratedt1dislets_metadata$Matching == "PBMC-matching","Clone.size"]
islets_clone_size_non_matching = Integratedt1dislets_metadata[Integratedt1dislets_metadata$Matching == "Not-matching","Clone.size"]
islets_wilcox_res_matching_vs_nonmatching = wilcox.test(islets_clone_size_matching,islets_clone_size_non_matching,alternative="greater")
islets_wilcox_res_matching_vs_nonmatching


# make boxplot
pdf(file = "Islets_pbmcmatching_clinalexpansion_significance.pdf", width = 4, height = 4)
p2 <- ggplot(Integratedt1dislets_metadata,aes(x=Matching,y=Clone.size, fill = Matching)) + geom_boxplot () + scale_fill_manual(values=c("gray","red"))+labs(y=expression("Log"[10]*"(Clone size)"))+ 
  theme(axis.title.x=element_blank(),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border= element_blank(),
        legend.position = "none",
        text = element_text(size = 16))+
  annotate(geom="text", x=1.2, y=1.7, label="Wilcoxon, p-value < 2.2e-16", color="black", size = 4)
p2
dev.off()




# cell count by group
library(data.table)
library(magrittr)
# extract meta data
metadata <- Integratedt1dislets_labelled@meta.data %>% as.data.table
# count the number of cells per unique combinations of "matcing clone number" and "seurat_clusters"

cell.num.by.matchingclone.seuratcluster <- metadata[, .N, by = c("Matching_pre_filter", "seurat_clusters")] %>% dcast(., Matching_pre_filter ~ seurat_clusters, value.var = "N")
write.csv(cell.num.by.matchingclone.seuratcluster, "t1dislets_cell.num.by.matchingclone.seuratcluster.csv", quote = F)


#count by clone size
cell.num.by.clonesize <- metadata[, .N, by = c("Matching_pre_filter", "Clone.size")] %>% dcast(., Matching_pre_filter ~ Clone.size, value.var = "N")
write.csv(cell.num.by.clonesize, "t1dislets_cell.num.by.clonesize.csv", quote = F)


#count by sample id
cell.num.by.sample <- metadata[, .N, by = c("sampleid", "seurat_clusters")] %>% dcast(., sampleid ~ seurat_clusters, value.var = "N")
write.csv(cell.num.by.sample, "t1dislets_cell.num.by.sample.csv", quote = F)


#### B5. DE pbmc matching T cells and non-matching T cells ####

#compare all cluster together
allcluster_matching_nonmatching_islets <- FindMarkers(Integratedt1dislets_labelled, ident.1 = "matching", ident.2 = "not_matching", group.by = 'Matching_pre_filter')

#write.csv(allcluster_matching_nonmatching_islets, file ="DE_t1dislets_allcluster_matching_nonmatching.csv")

#compare specific cluster
CD8_effector_memory_CTLs_matching_nonmatching_in_islets <- FindMarkers(Integratedt1dislets_labelled, ident.1 = "matching", ident.2 = "not_matching", group.by = 'Matching_pre_filter', subset.ident = "CD8 effector memory CTLs")

#write.csv(CD8_effector_memory_CTLs_matching_nonmatching_in_islets, file ="DE_t1dislets_CD8_effector_memory_T_cells_matching_nonmatching.csv")


CD4_CD25_Th1_cells_matching_nonmatching_in_islets <- FindMarkers(Integratedt1dislets_labelled, ident.1 = "matching", ident.2 = "not_matching", group.by = 'Matching_pre_filter', subset.ident = "CD4 CD25 Th1 cells")

#write.csv(CD4_CD25_Th1_cells_matching_nonmatching_in_islets, file ="DE_t1dislets_CD4_CD25_Th1_cells_matching_nonmatching.csv")

#plot violin 
pdf(file = "volcano_DE_pbmc-matching_vs_non-mathcing_in_islets.pdf", width = 12, height = 10)
islets_volcano <- EnhancedVolcano(allcluster_matching_nonmatching_islets,
                                lab = rownames(allcluster_matching_nonmatching_islets),
                                x = 'avg_log2FC',
                                y = 'p_val_adj',
                                title = 'DE_pbmc-matching_vs_non-mathcing_in_islets',
                                pCutoff = 10e-10,
                                FCcutoff = 1,
                                pointSize = 3.0,
                                labSize = 5.0,
                                drawConnectors = TRUE,
                                widthConnectors = 0.5)+theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border= element_blank())

islets_volcano
dev.off()


pdf(file = "volcano_DE_pbmc-matching_vs_non-mathcing_in_islets_CD8_effector_memory.pdf", width = 8, height = 8)
CD8_effector_memory_islets_volcano <- EnhancedVolcano(CD8_effector_memory_CTLs_matching_nonmatching_in_islets,
                                  lab = rownames(CD8_effector_memory_CTLs_matching_nonmatching_in_islets),
                                  x = 'avg_log2FC',
                                  y = 'p_val_adj',
                                  title = 'DE_pbmc-matching_vs_non-mathcing_in_islets_CD8_effector_memory',
                                  pCutoff = 10e-5,
                                  FCcutoff = .5,
                                  pointSize = 3.0,
                                  labSize = 5.0,
                                  drawConnectors = TRUE,
                                  widthConnectors = 0.5)+theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border= element_blank())

CD8_effector_memory_islets_volcano
dev.off()


pdf(file = "volcano_DE_pbmc-matching_vs_non-mathcing_in_islets_CD4_CD25_Th1.pdf", width = 8, height = 8)
CD4_CD25_Th1_islets_volcano <- EnhancedVolcano(CD4_CD25_Th1_cells_matching_nonmatching_in_islets,
                                                      lab = rownames(CD4_CD25_Th1_cells_matching_nonmatching_in_islets),
                                                      x = 'avg_log2FC',
                                                      y = 'p_val_adj',
                                                      title = 'DE_pbmc-matching_vs_non-mathcing_in_islets_CD8_effector_memory',
                                                      pCutoff = 10e-2,
                                                      FCcutoff = .5,
                                                      pointSize = 3.0,
                                                      labSize = 5.0,
                                                      drawConnectors = TRUE,
                                                      widthConnectors = 0.5)+theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border= element_blank())

CD4_CD25_Th1_islets_volcano
dev.off()


############ C. Analysis for non-T1D PBMC ###########

#### C1. load integrated rds ####
Integratednont1dpbmc <- readRDS("/Users/zahir/Documents/Tcell/seurat/1.DataProcessing/Integrated_rds/Integratednont1dpbmc.rds")

#For performing differential expression after integration, we switch back to the original data
DefaultAssay(Integratednont1dpbmc) <- "RNA"

#### C2. #Make UMAP plot ####
pdf(file = "nont1dpbmc_integrated_unlabeled_umap.pdf", width = 8, height = 6)
pc1 <- DimPlot(Integratednont1dpbmc, reduction = "umap", label = TRUE, repel = TRUE)
pc1
dev.off()


#### C3. Cluster cell population identification ####

#First we will use Findconserved markers function to identified conserved markers across different clusters. 
#using this loop for all 0-10 clusters
for (i in 0:11){
  marker_i <- FindConservedMarkers(Integratednont1dpbmc, ident.1 = i, grouping.var = "sampletype", verbose =TRUE)
  filename <- paste0("nont1dpbmc_cluster.", i,".csv")
  write.csv(marker_i, filename)
}

#Now we will identify cluster specific marker by differential expression of markers in a specific cluster
nont1dpbmc_findall_markers <-FindAllMarkers(Integratednont1dpbmc, min.pct = 0.25, logfc.threshold = 0.25)

write.csv(nont1dpbmc_findall_markers, "nont1d_pbmc_find_clusters_markers/nont1dpbmc_findall_markers.csv", quote = F)

#if there is any confusion about expression of any canonical marker genes then try to solve it by violin ploting or featureploting 

#plot by simple violin 
plots <- VlnPlot(Integratednont1dpbmc, features = c("Tmem176b", "Tmem176a", "Il12rb1", "Cd4"), group.by = "seurat_clusters",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

#FeaturePlot(object = Integratednont1dpbmc, features = "Cd4")

#staked violin plot in SeuratTheme
## remove the x-axis text and tick
## plot.margin to adjust the white space between each plot.
## ... pass any arguments to VlnPlot in Seurat
modify_vlnplot<- function(obj,
                          feature,
                          pt.size = 0,
                          plot.margin = margin(0, 0, 0, 0, "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, group.by = "seurat_clusters", pt.size = pt.size, ... )+
    ylab(feature) +
    theme(legend.position = "none",
          plot.title= element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_text(size = rel(1), angle = 0),
          axis.text.y = element_text(size = rel(1)),
          plot.margin = plot.margin )
  return(p)
}

# extract the max value of the y axis
extract_max<- function(p){
  ymax<- max(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range)
  return(ceiling(ymax))
}

## main function
StackedVlnPlot<- function(obj, features,
                          pt.size = 0, 
                          plot.margin = margin(0, 0, 0, 0, "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  
  # Add back x-axis title to bottom plot. patchwork is going to support this?
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5), axis.ticks.x = element_line(), legend.position = "right")
  
  # change the y-axis tick to only max value 
  ymaxs<- purrr::map_dbl(plot_list, extract_max)
  plot_list<- purrr::map2(plot_list, ymaxs, function(x,y) x + 
                            scale_y_continuous(breaks = c(y)) + 
                            expand_limits(y = y))
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


features<- c("Ccl5", "Ccr10", "Ccr6","Ccr7", "Cd28", "Cd38", "Cd4", "Cd44", "Cd69", "Cd7", "Cd74", "Cd8a", "Cd8b1", "Ctla2a", "Ctla4", "Cxcr5", "Fasl", "Foxp3", "Gata3", "Ifng", "Igfbp4", "Il12rb1", "Il17a", "Il1r1", "Il2ra", "Il2rb", "Il7r", "Itgb1", "Klrd1", "Lef1", "Ly6a", "Nkg7", "Rora", "Sell", "Tbx21", "Tmem176a", "Tmem176b")

pdf(file = "nont1dpbmc_canonical_violin_staked.pdf", width = 12, height = 16)
stackedplot_violin_nont1dpbmc <- StackedVlnPlot(obj = Integratednont1dpbmc, features = features)
stackedplot_violin_nont1dpbmc
dev.off()


#rename cluster ID 
DefaultAssay(Integratednont1dpbmc) <- "RNA"

Integratednont1dpbmc_labelled <- RenameIdents(Integratednont1dpbmc, `0` = "CD4 naive T cells", `1` = "CD8 central memory CTLs", `2` = "CD4 CTLA2 central memory T cells", `3` = "CD4 central memory T cells", `4` = "CD8 effector memory CTLs", `5` = "CD4 CD25 Foxp3 Tregs", `6` = "CD8 CTLA2a Trm CTLs", `7` = "CD4 CD74 central memory T cells", `8` = "DN CD7 NKT cells", `9` = "DN unclassified",`10` = "DP NKT cells", `11` = "DP Trm like")


nont1D_colors = c('CD4 naive T cells'='#F0A0FF','CD8 central memory CTLs'='#0075DC','CD4 central memory T cells'='#993F00','CD4 CD25 Foxp3 Tregs'='#191919','CD8 effector memory CTLs'='#005C31','CD8 CTLA2a Trm CTLs'='#2BCE48','DN CD7 NKT cells'='#FFCC99','CD4 CTLA2 central memory T cells'='#9DCC00','CD4 CD74 central memory T cells'='#C20088','DN unclassified'='#003380','DP NKT cells'='#FFA405','DP Trm like'='#FFA8BB')



#colored umap
pdf(file = "nont1dpbmc_labelled_umap_colored.pdf", width = 10, height = 6)
pc2 <- DimPlot(Integratednont1dpbmc_labelled,cols=nont1D_colors, reduction = "umap", label = FALSE, repel = TRUE, pt.size = 0.2, label.size = 5) + theme(axis.text = element_text(size = 20)) + theme(legend.text=element_text(size=16))
pc2
dev.off()




#Labeled umap plot 
pdf(file = "nont1dpbmc_labelled_umap.pdf", width = 10, height = 6)
pc2 <- DimPlot(Integratednont1dpbmc_labelled, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.3, label.size = 3)
pc2
dev.off()


# Supplementary Dotplot for canonical markers for supplementary with UMAP 
Idents(Integratednont1dpbmc_labelled) <- factor(Idents(Integratednont1dpbmc_labelled), levels = c("CD4 naive T cells", "CD8 central memory CTLs", "CD4 CTLA2 central memory T cells", "CD4 central memory T cells", "CD8 effector memory CTLs", "CD4 CD25 Foxp3 Tregs", "CD8 CTLA2a Trm CTLs", "CD4 CD74 central memory T cells", "DN CD7 NKT cells", "DN unclassified", "DP NKT cells", "DP Trm like"))


canonicalmarker <- c("Ccl5", "Ccr10", "Ccr6","Ccr7", "Cd28", "Cd38", "Cd4", "Cd44", "Cd69", "Cd7", "Cd74", "Cd8a", "Cd8b1", "Ctla2a", "Ctla4", "Cxcr5", "Fasl", "Foxp3", "Gata3", "Ifng", "Igfbp4", "Il12rb1", "Il17a", "Il1r1", "Il2ra", "Il2rb", "Il7r", "Itgb1", "Klrd1", "Lef1", "Ly6a", "Nkg7", "Rora", "Sell", "Tbx21", "Tmem176a", "Tmem176b")


pdf(file = "non-t1dpbmc_canonicaldotplot.pdf", width = 12, height = 6)
canonicaldotplot_nont1dpbmc <- DotPlot(Integratednont1dpbmc_labelled, features = canonicalmarker, cols = c("blue", "red", "green")) + RotatedAxis()
canonicaldotplot_nont1dpbmc
dev.off()


#Cell count in different cluster
cell.num <- table(Idents(Integratednont1dpbmc_labelled))

write.csv(cell.num, "nont1dpbmc_cellcountbycluster.csv", quote = F)


# cell count by group
library(data.table)
library(magrittr)
# extract meta data
metadata <- Integratednont1dpbmc_labelled@meta.data %>% as.data.table


#count by sample id
cell.num.by.sample <- metadata[, .N, by = c("sampleid", "seurat_clusters")] %>% dcast(., sampleid ~ seurat_clusters, value.var = "N")
write.csv(cell.num.by.sample, "nont1dpbmc_cell.num.by.sample.csv", quote = F)




### DIfferential expression ofislets  matching T cells in PBMC vs obmc matching T cells in islets 
#Sam did this analysis and found fowllwing DE genes 

DE_isletsmatchinginpbmc_vs_pbmcmatchinginislets <- read.csv(file = "MousePseudoBulk_DEgenes_unsorted.csv", header = TRUE)

#change the rowname with gene name that is located in the column1 
rownames(DE_isletsmatchinginpbmc_vs_pbmcmatchinginislets) = DE_isletsmatchinginpbmc_vs_pbmcmatchinginislets[,1]

#Change column name to make it same as previous plot 

colnames(DE_isletsmatchinginpbmc_vs_pbmcmatchinginislets)[2] <- "avg_log2FC"

colnames(DE_isletsmatchinginpbmc_vs_pbmcmatchinginislets)[6] <- "p_val_adj"



pdf(file = "volcano_DE_isletsmatchinginpbmc_vs_pbmcmatchinginislets.pdf", width = 12, height = 10)
DE_isletsmatchinginpbmc_vs_pbmcmatchinginislets_volcano <- EnhancedVolcano(DE_isletsmatchinginpbmc_vs_pbmcmatchinginislets,
                                                                           lab = rownames(DE_isletsmatchinginpbmc_vs_pbmcmatchinginislets),
                                                                           x = 'avg_log2FC',
                                                                           y = 'p_val_adj',
                                                                           title = 'DE_isletsmatchinginpbmc_vs_pbmcmatchinginislets',
                                                                           pCutoff = 10e-10,
                                                                           FCcutoff = 1.5,
                                                                           pointSize = 3.0,
                                                                           labSize = 5.0,
                                                                           drawConnectors = TRUE,
                                                                           widthConnectors = 0.5)+theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border= element_blank())

DE_isletsmatchinginpbmc_vs_pbmcmatchinginislets_volcano
dev.off()


#It seems lot of the significantly associated genes are tissue-specific gene contamination 
#so we will now identify whether they are contamination or not 
#violin plot for islets enriched feature

islets_contaminant_genes <- c("Ctrb1", "Try5", "Cela2a", "Gcg", "Clps", "Try4", "Pnlip", "Cela3b", "Cpa1", "2210010C04Rik", "Sycn", "Ins2", "Rnase1", "Reg1", "Ctrc", "Sst", "Cel", "Cpb1", "Ins1", "Pnliprp1", "Zg16", "Iapp", "Ppy", "Cpa2", "Pla2g1b", "Prss2", "Cela1", "Ctrl", "Gm20481", "Hspa1b", "Lag3", "Pdcd1", "Hspa1a", "Tnfsf11")


DE_isletsmatchinginpbmc_vs_pbmcmatchinginislets_contaminat_removed <- DE_isletsmatchinginpbmc_vs_pbmcmatchinginislets[!(row.names(DE_isletsmatchinginpbmc_vs_pbmcmatchinginislets) %in% islets_contaminant_genes),]

#Now plot after removing contaminants
pdf(file = "volcano_DE_isletsmatchinginpbmc_vs_pbmcmatchinginislets_contaminant_removed.pdf", width = 12, height = 10)
DE_isletsmatchinginpbmc_vs_pbmcmatchinginislets_contaminant_removed_volcano <- EnhancedVolcano(DE_isletsmatchinginpbmc_vs_pbmcmatchinginislets_contaminat_removed, lab = rownames(DE_isletsmatchinginpbmc_vs_pbmcmatchinginislets_contaminat_removed), x = 'avg_log2FC', y = 'p_val_adj', title = 'DE_isletsmatchinginpbmc_vs_pbmcmatchinginislets_contaminat_removed', pCutoff = 10e-10, FCcutoff = 1, pointSize = 3.0, labSize = 5.0, drawConnectors = TRUE, widthConnectors = 0.5)+theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border= element_blank())

DE_isletsmatchinginpbmc_vs_pbmcmatchinginislets_contaminant_removed_volcano
dev.off()





#gene set enrichment analysis
#https://github.com/erilu/single-cell-rnaseq-analysis#gene-set-enrichment-analysis-gsea-across-groups


#### load integrated rds ####
Integratedt1dpbmc <- readRDS("/Users/zahir/Documents/Tcell/seurat/1.DataProcessing/Integrated_rds/Integratedt1dpbmc.rds")

#For performing differential expression after integration, we switch back to the original data
DefaultAssay(Integratedt1dpbmc) <- "RNA"

library(fgsea)
library(tidyverse)

#first compare all cluster together
# use FindMarkers with min.pct = 0.1, and no logfc threshold to rank all the genes available
GSEA_list_pbmc_allcluster_matching_nonmatching <- FindMarkers(Integratedt1dpbmc_labelled, ident.1 = "matching", ident.2 = "not_matching", group.by = 'Matching_pre_filter', min.pct = 0.1, logfc.threshold = 0)

#to write the results in csv
write.csv(GSEA_list_pbmc_allcluster_matching_nonmatching, file ="GSEA_list_pbmc_allcluster_matching_nonmatching.csv")

#fix row name issue
Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching <- cbind(Gene.name = rownames(GSEA_list_pbmc_allcluster_matching_nonmatching), GSEA_list_pbmc_allcluster_matching_nonmatching)


# order list, pull out gene name and log2fc, and convert genes to uppercase
Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching <- Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching[order(Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching$avg_log2FC, decreasing = T),]


Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching$Gene.name = toupper(Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching$Gene.name)


Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching <- Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching[,c("Gene.name", "avg_log2FC")]
rownames(Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching) <- NULL

# write to table if want to save for GSEA application
save(Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching, file = "GSEA_list_pbmc_allcluster_matching_nonmatching_ranked.Robj")
write.table(Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching, file = "GSEA_list_pbmc_allcluster_matching_nonmatching_ranked.rnk", sep = "\t", row.names = F, quote = F)

#First, we have to download the pathway files from the GSEA MSigDB webpage at: https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp and load them in using fgsea::gmtPathways(). A commonly-used set of pathways is the HALLMARK pathway set, which I will use below.

# read in file containing lists of genes for each pathway
hallmark_pathway <- gmtPathways("GSEA/h.all.v7.5.1.symbols.gmt.txt")
head(names(hallmark_pathway))


#then, we have to turn our ranked list into a vector in which the avg_logFC is named with the gene name, as well as get rid of any NA values or duplicate gene entries.

# formats the ranked list for the fgsea() function
prepare_ranked_list <- function(ranked_list) {
  # if duplicate gene names present, average the values
  if( sum(duplicated(ranked_list$Gene.name)) > 0) {
    ranked_list <- aggregate(.~Gene.name, FUN = mean, data = ranked_list)
    ranked_list <- ranked_list[order(ranked_list$avg_log2FC, decreasing = T),]
  }
  # omit rows with NA values
  ranked_list <- na.omit(ranked_list)
  # turn the dataframe into a named vector
  ranked_list <- tibble::deframe(ranked_list)
  ranked_list
}

Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching <- prepare_ranked_list(Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching)
head(Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching)


#Now that we have a cleaned-up vector of avg_logFC named by gene, we can plug it into the fgsea::fgsea() function along with the hallmark pathway object, which generates a table of results containing the enrichment scores associated with each pathway:
  
  # generate GSEA result table using fgsea() by inputting the pathway list and ranked list
fgsea_results_Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching <- fgsea(pathways = hallmark_pathway,
                         stats = Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching,
                         minSize = 15,
                         maxSize = 500,
                         nperm= 1000)

fgsea_results_Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching %>% arrange (desc(NES)) %>% select (pathway, padj, NES) %>% head()




plot_enrichment <- function (geneset, pathway, ranked_list) {
  plotEnrichment(geneset[[pathway]], ranked_list)+labs (title = pathway)
}

# download fgsea result table and plot in r or other program like prism

table1 <- data.frame(fgsea_results_Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching)

write.csv(table1,"fgsea_results_Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching.csv")

fwrite(table1, file ="fgsea_results_Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching.csv")


# example of pathway highly enriched in ADAR
plot_enrichment(hallmark_pathway, "HALLMARK_IL2_STAT5_SIGNALING" , Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching)

#save 5 by 6 

plot_enrichment(hallmark_pathway, "HALLMARK_INTERFERON_GAMMA_RESPONSE" , Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching)


#The fsgea_results table above is useful, but it would be nicer to view in a more graphical format. Below, we construct a "waterfall" plot using the function waterfall_plot(), which is basically a sideways bar plot of enrichment scores for each of the pathways. The significantly enriched pathways are highlighted in a different color.
library(ggplot2)
library(tidyverse)

waterfall_plot <- function (fgsea_results_Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching, graph_title) {
  fgsea_results_Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching %>%
    mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>%
    ggplot( aes(reorder(short_name,NES), NES)) +
    geom_bar(stat= "identity", aes(fill = padj<0.05))+
    coord_flip()+
    labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = graph_title)+
    theme(axis.text.y = element_text(size = 7),
          plot.title = element_text(hjust = 1))
}


waterfall_plot(fgsea_results_Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching, "Pathways enriched in islets-matching cells vs non-matching cells in blood")

#save 5 height and 6 weidth 


#We can also perform the analysis in a specific cluster
cluster5_CD8_EM_CTL_GSEA_list_pbmc_allcluster_matching_nonmatching <- FindMarkers(Integratedt1dpbmc_labelled, ident.1 = "matching", ident.2 = "not_matching", group.by = 'Matching_pre_filter', subset.ident = "CD8 effector memory CTLs", min.pct = 0.1, logfc.threshold = 0)

#to write files
write.csv(cluster5_CD8_EM_CTL_GSEA_list_pbmc_allcluster_matching_nonmatching, file ="cluster5_CD8_EM_CTL_GSEA_list_pbmc_allcluster_matching_nonmatching.csv")





##GSEA for islets 

#### load integrated rds ####
Integratedt1dislets <- readRDS("/Users/zahir/Documents/Tcell/seurat/1.DataProcessing/Integrated_rds/Integratedislets.rds")

#For performing differential expression after integration, we switch back to the original data
DefaultAssay(Integratedt1dislets) <- "RNA"

library(fgsea)
library(tidyverse)

#first compare all cluster together
# use FindMarkers with min.pct = 0.1, and no logfc threshold to rank all the genes available
GSEA_list_islets_allcluster_matching_nonmatching <- FindMarkers(Integratedt1dislets_labelled, ident.1 = "matching", ident.2 = "not_matching", group.by = 'Matching_pre_filter', min.pct = 0.1, logfc.threshold = 0)

#to write the results in csv
write.csv(GSEA_list_islets_allcluster_matching_nonmatching, file ="GSEA_list_islets_allcluster_matching_nonmatching.csv")

#fix row name issue
Rownames_GSEA_list_islets_allcluster_matching_nonmatching <- cbind(Gene.name = rownames(GSEA_list_islets_allcluster_matching_nonmatching), GSEA_list_islets_allcluster_matching_nonmatching)


# order list, pull out gene name and log2fc, and convert genes to uppercase
Rownames_GSEA_list_islets_allcluster_matching_nonmatching <- Rownames_GSEA_list_islets_allcluster_matching_nonmatching[order(Rownames_GSEA_list_islets_allcluster_matching_nonmatching$avg_log2FC, decreasing = T),]


Rownames_GSEA_list_islets_allcluster_matching_nonmatching$Gene.name = toupper(Rownames_GSEA_list_islets_allcluster_matching_nonmatching$Gene.name)


Rownames_GSEA_list_islets_allcluster_matching_nonmatching <- Rownames_GSEA_list_islets_allcluster_matching_nonmatching[,c("Gene.name", "avg_log2FC")]
rownames(Rownames_GSEA_list_islets_allcluster_matching_nonmatching) <- NULL


# write to table if want to save for GSEA application
save(Rownames_GSEA_list_islets_allcluster_matching_nonmatching, file = "Rownames_GSEA_list_islets_allcluster_matching_nonmatching.Robj")
write.table(Rownames_GSEA_list_pbmc_allcluster_matching_nonmatching, file = "Rownames_GSEA_list_islets_allcluster_matching_nonmatching.rnk", sep = "\t", row.names = F, quote = F)

#First, we have to download the pathway files from the GSEA MSigDB webpage at: https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp and load them in using fgsea::gmtPathways(). A commonly-used set of pathways is the HALLMARK pathway set, which I will use below.

# read in file containing lists of genes for each pathway
hallmark_pathway <- gmtPathways("GSEA/h.all.v7.5.1.symbols.gmt.txt")
head(names(hallmark_pathway))


#then, we have to turn our ranked list into a vector in which the avg_logFC is named with the gene name, as well as get rid of any NA values or duplicate gene entries.

# formats the ranked list for the fgsea() function
prepare_ranked_list <- function(ranked_list) {
  # if duplicate gene names present, average the values
  if( sum(duplicated(ranked_list$Gene.name)) > 0) {
    ranked_list <- aggregate(.~Gene.name, FUN = mean, data = ranked_list)
    ranked_list <- ranked_list[order(ranked_list$avg_log2FC, decreasing = T),]
  }
  # omit rows with NA values
  ranked_list <- na.omit(ranked_list)
  # turn the dataframe into a named vector
  ranked_list <- tibble::deframe(ranked_list)
  ranked_list
}

Rownames_GSEA_list_islets_allcluster_matching_nonmatching <- prepare_ranked_list(Rownames_GSEA_list_islets_allcluster_matching_nonmatching)
head(Rownames_GSEA_list_islets_allcluster_matching_nonmatching)


#Now that we have a cleaned-up vector of avg_logFC named by gene, we can plug it into the fgsea::fgsea() function along with the hallmark pathway object, which generates a table of results containing the enrichment scores associated with each pathway:

# generate GSEA result table using fgsea() by inputting the pathway list and ranked list
fgsea_results_Rownames_GSEA_list_islets_allcluster_matching_nonmatching <- fgsea(pathways = hallmark_pathway,
                                                                               stats = Rownames_GSEA_list_islets_allcluster_matching_nonmatching,
                                                                               minSize = 15,
                                                                               maxSize = 500,
                                                                               nperm= 1000)

fgsea_results_Rownames_GSEA_list_islets_allcluster_matching_nonmatching %>% arrange (desc(NES)) %>% select (pathway, padj, NES) %>% head()




plot_enrichment <- function (geneset, pathway, ranked_list) {
  plotEnrichment(geneset[[pathway]], ranked_list)+labs (title = pathway)
}

# download fgsea result table and plot in r or other program like prism

table2 <- data.frame(fgsea_results_Rownames_GSEA_list_islets_allcluster_matching_nonmatching)

fwrite(table2, file ="fgsea_results_Rownames_GSEA_list_islets_allcluster_matching_nonmatching.csv")


# example of pathway highly enriched in ADAR
plot_enrichment(hallmark_pathway, "HALLMARK_IL2_STAT5_SIGNALING" , Rownames_GSEA_list_islets_allcluster_matching_nonmatching)

#save 5 by 6 

plot_enrichment(hallmark_pathway, "HALLMARK_INTERFERON_GAMMA_RESPONSE" , Rownames_GSEA_list_islets_allcluster_matching_nonmatching)


#The fsgea_results table above is useful, but it would be nicer to view in a more graphical format. Below, we construct a "waterfall" plot using the function waterfall_plot(), which is basically a sideways bar plot of enrichment scores for each of the pathways. The significantly enriched pathways are highlighted in a different color.
library(ggplot2)
library(tidyverse)

waterfall_plot <- function (fgsea_results_Rownames_GSEA_list_islets_allcluster_matching_nonmatching, graph_title) {
  fgsea_results_Rownames_GSEA_list_islets_allcluster_matching_nonmatching %>%
    mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>%
    ggplot( aes(reorder(short_name,NES), NES)) +
    geom_bar(stat= "identity", aes(fill = padj<0.05))+
    coord_flip()+
    labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = graph_title)+
    theme(axis.text.y = element_text(size = 7),
          plot.title = element_text(hjust = 1))
}


waterfall_plot(fgsea_results_Rownames_GSEA_list_islets_allcluster_matching_nonmatching, "Pathways enriched in blood-matching cells vs non-matching cells in islets")

#save 5 height and 6 weidth 





## GSEA analysis from psudobulk output 
library(fgsea)
library(tidyverse)

#### read csv file 
GSEA_list_between_blood_islets_matching_clone <- read.csv(file = 'MousePseudoBulk_DEgenes_unsorted.csv')

#rename gene column
names(GSEA_list_between_blood_islets_matching_clone)[1] <- 'Gene.name'

names(GSEA_list_between_blood_islets_matching_clone)[2] <- 'avg_log2FC'

# order list, pull out gene name and log2fc, and convert genes to uppercase
GSEA_list_between_blood_islets_matching_clone <- GSEA_list_between_blood_islets_matching_clone[order(GSEA_list_between_blood_islets_matching_clone$avg_log2FC, decreasing = T),]


GSEA_list_between_blood_islets_matching_clone$Gene.name = toupper(GSEA_list_between_blood_islets_matching_clone$Gene.name)


GSEA_list_between_blood_islets_matching_clone <- GSEA_list_between_blood_islets_matching_clone[,c("Gene.name", "avg_log2FC")]
rownames(GSEA_list_between_blood_islets_matching_clone) <- NULL


# write to table if want to save for GSEA application
save(GSEA_list_between_blood_islets_matching_clone, file = "GSEA_list_between_blood_islets_matching_clone.Robj")
write.table(GSEA_list_between_blood_islets_matching_clone, file = "GSEA_list_between_blood_islets_matching_clone.rnk", sep = "\t", row.names = F, quote = F)

#First, we have to download the pathway files from the GSEA MSigDB webpage at: https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp and load them in using fgsea::gmtPathways(). A commonly-used set of pathways is the HALLMARK pathway set, which I will use below.

# read in file containing lists of genes for each pathway
hallmark_pathway <- gmtPathways("GSEA/h.all.v7.5.1.symbols.gmt.txt")
head(names(hallmark_pathway))


#then, we have to turn our ranked list into a vector in which the avg_logFC is named with the gene name, as well as get rid of any NA values or duplicate gene entries.

# formats the ranked list for the fgsea() function
prepare_ranked_list <- function(ranked_list) {
  # if duplicate gene names present, average the values
  if( sum(duplicated(ranked_list$Gene.name)) > 0) {
    ranked_list <- aggregate(.~Gene.name, FUN = mean, data = ranked_list)
    ranked_list <- ranked_list[order(ranked_list$avg_log2FC, decreasing = T),]
  }
  # omit rows with NA values
  ranked_list <- na.omit(ranked_list)
  # turn the dataframe into a named vector
  ranked_list <- tibble::deframe(ranked_list)
  ranked_list
}

GSEA_list_between_blood_islets_matching_clone <- prepare_ranked_list(GSEA_list_between_blood_islets_matching_clone)
head(GSEA_list_between_blood_islets_matching_clone)


#Now that we have a cleaned-up vector of avg_logFC named by gene, we can plug it into the fgsea::fgsea() function along with the hallmark pathway object, which generates a table of results containing the enrichment scores associated with each pathway:

# generate GSEA result table using fgsea() by inputting the pathway list and ranked list
fgsea_results_GSEA_list_between_blood_islets_matching_clone <- fgsea(pathways = hallmark_pathway,
                                                                                 stats = GSEA_list_between_blood_islets_matching_clone,
                                                                                 minSize = 15,
                                                                                 maxSize = 500,
                                                                                 nperm= 1000)

fgsea_results_GSEA_list_between_blood_islets_matching_clone %>% arrange (desc(NES)) %>% select (pathway, padj, NES) %>% head()




plot_enrichment <- function (geneset, pathway, ranked_list) {
  plotEnrichment(geneset[[pathway]], ranked_list)+labs (title = pathway)
}

# download fgsea result table and plot in r or other program like prism

table3 <- data.frame(fgsea_results_GSEA_list_between_blood_islets_matching_clone)

fwrite(table3, file ="fgsea_results_GSEA_list_between_blood_islets_matching_clone.csv")


# example of pathway highly enriched in ADAR
plot_enrichment(hallmark_pathway, "HALLMARK_IL2_STAT5_SIGNALING" , GSEA_list_between_blood_islets_matching_clone)

#save 5 by 6 

plot_enrichment(hallmark_pathway, "HALLMARK_INTERFERON_GAMMA_RESPONSE" , GSEA_list_between_blood_islets_matching_clone)


#The fsgea_results table above is useful, but it would be nicer to view in a more graphical format. Below, we construct a "waterfall" plot using the function waterfall_plot(), which is basically a sideways bar plot of enrichment scores for each of the pathways. The significantly enriched pathways are highlighted in a different color.
library(ggplot2)
library(tidyverse)

waterfall_plot <- function (fgsea_results_GSEA_list_between_blood_islets_matching_clone, graph_title) {
  fgsea_results_GSEA_list_between_blood_islets_matching_clone %>%
    mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>%
    ggplot( aes(reorder(short_name,NES), NES)) +
    geom_bar(stat= "identity", aes(fill = padj<0.05))+
    coord_flip()+
    labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = graph_title)+
    theme(axis.text.y = element_text(size = 7),
          plot.title = element_text(hjust = 1))
}


waterfall_plot(fgsea_results_GSEA_list_between_blood_islets_matching_clone, "Pathways enriched in blood-matching  vs islets-matching")


#save 5 height and 6 weidth 


##GSEA for all mergged sample 


Integrated_allmergedsamples <- readRDS("/Users/zahir/Documents/Tcell/seurat/1.DataProcessing/Integrated_rds/Integrated_allmergedsamples.rds")

#For performing differential expression after integration, we switch back to the original data
DefaultAssay(Integrated_allmergedsamples) <- "RNA"

library(fgsea)
library(tidyverse)

#first compare all cluster together
# use FindMarkers with min.pct = 0.1, and no logfc threshold to rank all the genes available
GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc <- FindMarkers(Integrated_allmergedsamples, ident.1 = "t1dpbmc", ident.2 = "nont1dpbmc", group.by = 'sampletype', min.pct = 0.1, logfc.threshold = 0)



#to write the results in csv
write.csv(GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc, file ="GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc.csv")

#fix row name issue
Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc <- cbind(Gene.name = rownames(GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc), GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc)


# order list, pull out gene name and log2fc, and convert genes to uppercase
Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc <- Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc[order(Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc$avg_log2FC, decreasing = T),]


Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc$Gene.name = toupper(Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc$Gene.name)


Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc <- Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc[,c("Gene.name", "avg_log2FC")]
rownames(Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc) <- NULL

# write to table if want to save for GSEA application
save(Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc, file = "GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc_ranked.Robj")
write.table(Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc, file = "GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc_ranked.rnk", sep = "\t", row.names = F, quote = F)

#First, we have to download the pathway files from the GSEA MSigDB webpage at: https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp and load them in using fgsea::gmtPathways(). A commonly-used set of pathways is the HALLMARK pathway set, which I will use below.

# read in file containing lists of genes for each pathway
hallmark_pathway <- gmtPathways("GSEA/h.all.v7.5.1.symbols.gmt.txt")
head(names(hallmark_pathway))


#then, we have to turn our ranked list into a vector in which the avg_logFC is named with the gene name, as well as get rid of any NA values or duplicate gene entries.

# formats the ranked list for the fgsea() function
prepare_ranked_list <- function(ranked_list) {
  # if duplicate gene names present, average the values
  if( sum(duplicated(ranked_list$Gene.name)) > 0) {
    ranked_list <- aggregate(.~Gene.name, FUN = mean, data = ranked_list)
    ranked_list <- ranked_list[order(ranked_list$avg_log2FC, decreasing = T),]
  }
  # omit rows with NA values
  ranked_list <- na.omit(ranked_list)
  # turn the dataframe into a named vector
  ranked_list <- tibble::deframe(ranked_list)
  ranked_list
}

Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc <- prepare_ranked_list(Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc)
head(Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc)


#Now that we have a cleaned-up vector of avg_logFC named by gene, we can plug it into the fgsea::fgsea() function along with the hallmark pathway object, which generates a table of results containing the enrichment scores associated with each pathway:

# generate GSEA result table using fgsea() by inputting the pathway list and ranked list
fgsea_results_Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc <- fgsea(pathways = hallmark_pathway,
                                                                               stats = Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc,
                                                                               minSize = 15,
                                                                               maxSize = 500,
                                                                               nperm= 1000)

fgsea_results_Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc %>% arrange (desc(NES)) %>% select (pathway, padj, NES) %>% head()




plot_enrichment <- function (geneset, pathway, ranked_list) {
  plotEnrichment(geneset[[pathway]], ranked_list)+labs (title = pathway)
}

# download fgsea result table and plot in r or other program like prism

table1 <- data.frame(fgsea_results_Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc)

fwrite(table1, file ="fgsea_results_Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc.csv")


# example of pathway highly enriched in ADAR
plot_enrichment(hallmark_pathway, "HALLMARK_TNFA_SIGNALING_VIA_NFKB" , Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc)

#save 5 by 6 

plot_enrichment(hallmark_pathway, "HALLMARK_TGF_BETA_SIGNALING" , Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc)


plot_enrichment(hallmark_pathway, "HALLMARK_INTERFERON_ALPHA_RESPONSE" , Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc)



plot_enrichment(hallmark_pathway, "HALLMARK_INTERFERON_GAMMA_RESPONSE" , Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc)






#The fsgea_results table above is useful, but it would be nicer to view in a more graphical format. Below, we construct a "waterfall" plot using the function waterfall_plot(), which is basically a sideways bar plot of enrichment scores for each of the pathways. The significantly enriched pathways are highlighted in a different color.
library(ggplot2)
library(tidyverse)

waterfall_plot <- function (fgsea_results_Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc, graph_title) {
  fgsea_results_Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc %>%
    mutate(short_name = str_split_fixed(pathway, "_",2)[,2])%>%
    ggplot( aes(reorder(short_name,NES), NES)) +
    geom_bar(stat= "identity", aes(fill = padj<0.05))+
    coord_flip()+
    labs(x = "Hallmark Pathway", y = "Normalized Enrichment Score", title = graph_title)+
    theme(axis.text.y = element_text(size = 7),
          plot.title = element_text(hjust = 1))
}


waterfall_plot(fgsea_results_Rownames_GSEA_list_Integrated_allmergedsamples_t1d_nonT1D_pbmc, "Pathways enriched in non-T1D-blood vs T1D-blood")

#save 5 height and 6 weidth 

