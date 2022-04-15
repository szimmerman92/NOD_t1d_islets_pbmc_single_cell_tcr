library(monocle3)
library(Seurat)
library(SeuratData)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(magrittr)

####1.load integrated rds ####
Integratedt1dislets <- readRDS("/Users/zahir/Documents/Tcell/seurat/1.DataProcessing/Integrated_rds/Integratedislets.rds")

#For performing differential expression after integration, we switch back to the original data
DefaultAssay(Integratedt1dislets) <- "RNA"

####2.identify the pseudotime root cluster based on lowest entropy value in a particular cluster ####
#install 
#library(remotes)
#remotes::install_github("rnabioco/scbp")
library(scbp)

entropy <- calc_diversity(Integratedt1dislets, sample_id = "sampleid", group_id = "seurat_clusters")
# extract entropy for each cluster
entropy_cluster_df = entropy@meta.data[,c("seurat_clusters","entropy")]
entropy_cluster_df = unique(entropy_cluster_df)
entropy_cluster_df = entropy_cluster_df[order(entropy_cluster_df$seurat_clusters),]

pdf(file = "Entropy.pdf", width = 6, height = 4)
entropy <- ggplot(entropy_cluster_df,aes(x=reorder(seurat_clusters,entropy),y=entropy)) + geom_bar(stat="identity", fill = "brown2")+ 
  theme(axis.title.x=element_blank(),
        axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border= element_blank(),
        legend.position = "none",
        text = element_text(size = 20),
        axis.text.y = element_text(color = "black", size =16),
        axis.text.x = element_text(color = "black", size =16))
entropy
dev.off()



####3.Subset CD4 and Cd8 clusters into two seperate seurat objects ####
cd4_islets_cluster <- subset(Integratedt1dislets, subset=seurat_clusters %in% c(0, 1, 3, 4, 6, 7, 9,13,15))
DefaultAssay(cd4_islets_cluster) <- "RNA"

#label cluster
cd4_islets_cluster <- RenameIdents(cd4_islets_cluster, `0` = "Lef1 memory", `1` = "CD29 CD25", `3` = "Central memory", `4` = "CD25 CTLA4 Foxp3-", `6` = "CD25 Th1", `7` = "CD25 CTLA4 Foxp3+", `9` = "CTLA4 memory", `13` = "Th17 like", `15` = "Th17 memory like")

pdf(file = "cd4_seurat_cluster.pdf", width = 6, height = 4)
cd4_seurat_cluster <- DimPlot(cd4_islets_cluster, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend() + theme(axis.text = element_text(size = 16)) + theme(legend.text=element_text(size=8))
cd4_seurat_cluster
dev.off()


cd8_islets_cluster <- subset(Integratedt1dislets, subset=seurat_clusters %in% c(2, 5, 8, 11))
DefaultAssay(cd8_islets_cluster) <- "RNA"

#label cluster
cd8_islets_cluster <- RenameIdents(cd8_islets_cluster, `2` = "Ly6a central memory", `5` = "Effector memory", `8` = "Central memory", `11` = "CTLA2a Trm")

pdf(file = "cd8_seurat_cluster.pdf", width = 5, height = 4)
cd8_seurat_cluster <- DimPlot(cd8_islets_cluster, reduction = "umap", label = TRUE, repel = TRUE) + NoLegend() + theme(axis.text = element_text(size = 16)) + theme(legend.text=element_text(size=8))
cd8_seurat_cluster
dev.off()



####4. perform monocle clustering and learn graph ####
#CD4
cd4.cds <- as.cell_data_set(cd4_islets_cluster)
cd4.cds <- preprocess_cds(cd4.cds, method = "PCA")
cd4.cds <- cluster_cells(cds = cd4.cds, reduction_method = "UMAP")

#!#!!#! NOTE YOU WANT TO CHOOSE A SINGLE PARTITION HERE BEFORE MOVING ON!!!
cd4.cds <- learn_graph(cd4.cds, use_partition = TRUE)

pdf("cd4_monocle_clusters.pdf",width=4,height=5)
plot_cells(cd4.cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
dev.off()
pdf("cd4_monocle_partitions.pdf",width=4,height=5)
plot_cells(cd4.cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
dev.off()

pdf(file = "cd4_monocl_cluster_trajectory.pdf", width = 4, height = 3)
cd4_monocl_cluster_trajectory <- plot_cells(cd4.cds, color_cells_by = "ident", show_trajectory_graph = TRUE, label_cell_groups=TRUE, label_leaves=FALSE,label_branch_points=FALSE, label_roots = FALSE,trajectory_graph_color = "gray38", cell_size = 0.25)+ theme(axis.text = element_text(size = 16)) + theme(legend.text=element_text(size=14))
cd4_monocl_cluster_trajectory
dev.off()



#plot_cells(cd4.cds, color_cells_by = "cluster", show_trajectory_graph = FALSE,label_cell_groups=FALSE)

#plot_cells(cd4.cds, color_cells_by = "seurat_clusters", show_trajectory_graph = FALSE,label_cell_groups=FALSE)

#DimPlot(cd8_islets_cluster, reduction = "umap", label = TRUE, repel = TRUE)


#CD8
cd8.cds <- as.cell_data_set(cd8_islets_cluster)
cd8.cds <- preprocess_cds(cd8.cds, method = "PCA")
cd8.cds <- cluster_cells(cds = cd8.cds, reduction_method = "UMAP")
cd8.cds <- learn_graph(cd8.cds, use_partition = FALSE) # DON'T SEPARATE PARTITIONS


pdf(file = "cd8_monocl_cluster_trajectory.pdf", width = 4, height = 3)
cd8_monocl_cluster_trajectory <- plot_cells(cd8.cds, color_cells_by = "ident", show_trajectory_graph = TRUE,label_cell_groups=TRUE, label_leaves=FALSE,label_branch_points=FALSE, label_roots = FALSE,trajectory_graph_color = "gray38")+ theme(axis.text = element_text(size = 16)) + theme(legend.text=element_text(size=14))
cd8_monocl_cluster_trajectory
dev.off()



pdf("cd8_monocle_partitions.pdf",width=4,height=5)
plot_cells(cd8.cds, color_cells_by = "partition", show_trajectory_graph = FALSE) 
dev.off()
pdf("cd8_monocle_clusters.pdf",width=4,height=5)
plot_cells(cd8.cds, color_cells_by = "cluster", show_trajectory_graph = FALSE) 
dev.off()


#### 5.#Order cells ####
#CD$
#we are asigning cluster 3 as CD4 cells root
cd4.cds <- order_cells(cd4.cds, root_cells = rownames(Integratedt1dislets@meta.data)[Integratedt1dislets@meta.data$seurat_clusters == 3])


#plot trajectories colored by pseudotime
#plot trajectories colored by pseudotime
pdf(file = "cd4.root3_trajectory.pdf", width = 4, height = 2.5)
cd4.root3_trajectory <- plot_cells(cd4.cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           label_roots = FALSE,
           trajectory_graph_color = "gray38") + theme(axis.text = element_text(size = 16)) + theme(legend.text=element_text(size=14))
cd4.root3_trajectory
dev.off()



#CD8
#we are asigning cluster 8 as CD8 cells root
#there was some issues with automatic rooting by the following command

### THERE IS NOT CLUSTER 8?????????
cd8.cds <- order_cells(cd8.cds, root_cells = rownames(Integratedt1dislets@meta.data)[Integratedt1dislets@meta.data$seurat_clusters == 8])


#plot trajectories colored by pseudotime
pdf(file = "cd8.root8_trajectory.pdf", width = 4, height = 2.5)
cd8.root8_trajectory <- plot_cells(cd8.cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           label_roots = FALSE,
           trajectory_graph_color = "gray38") + theme(axis.text = element_text(size = 16)) + theme(legend.text=element_text(size=14))
cd8.root8_trajectory
dev.off()

#more pseudotime plot options
#Extract the pseudotime values and add to the Seurat object
Integratedt1dislets <- AddMetaData(
  object = Integratedt1dislets,
  metadata = cd4.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "CD4"
)

Integratedt1dislets <- AddMetaData(
  object = Integratedt1dislets,
  metadata = cd8.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "CD8"
)
#plot
cd4_clusters = c(0, 1, 3, 4, 6, 7, 9,13,15)
cd8_clusters = c(2, 5, 8, 11)
cd4_cd8_clusters = c(cd4_clusters,cd8_clusters)
cd4_cd8_islets_cluster <- subset(Integratedt1dislets, subset=seurat_clusters %in% cd4_cd8_clusters)

pdf("CD4_CD8_islet_pseudotime_seurat.pdf",width=7,height=3)
FeaturePlot(cd4_cd8_islets_cluster, c("CD4", "CD8"), pt.size = 0.1) & scale_color_viridis_c()
dev.off()

####6. Finding genes that change as a function of pseudotime ####
#Identify top genes 
#Identify genes
#CD4

# remember to do the below before running
#trace('calculateLW', edit = T, where = asNamespace("monocle3"))
# change Matrix::rBind to rbind

cd4.cds_graph_test_res <- graph_test(cd4.cds,
                                       neighbor_graph = "principal_graph",
                                       cores = 6)

head(cd4.cds_graph_test_res)

cd4.cds_graph_test_res_subset <- subset(cd4.cds_graph_test_res, q_value < 0.05)
cd4.cds_graph_test_res_subset_ordered = cd4.cds_graph_test_res_subset[order(cd4.cds_graph_test_res_subset$p_value,(-1*cd4.cds_graph_test_res_subset$morans_I)),]

 contamination_genes = c('Ctrb1', 'Try5', 'Cela2a', 'Gcg', 'Clps', 'Try4', 'Pnlip', 'Cela3b', 'Cpa1', '2210010C04Rik', 'Sycn', 'Ins2', 'Rnase1', 'Reg1', 'Ctrc', 'Sst', 'Cel', 'Cpb1', 'Ins1', 'Pnliprp1', 'Zg16', 'Iapp', 'Ppy', 'Cpa2', 'Pla2g1b', 'Prss2', 'Cela1', 'Ctrl', 'Gm20481', 'Hspa1b', 'Lag3', 'Pdcd1', 'Hspa1a', 'Tnfsf11')

cd4.cds_graph_test_res_subset_ordered=  cd4.cds_graph_test_res_subset_ordered[!rownames(cd4.cds_graph_test_res_subset_ordered)%in%contamination_genes,]
 

cd4.graph_deg_ids <- row.names(cd4.cds_graph_test_res_subset_ordered)

#to write files
write.csv(cd4.cds_graph_test_res_subset_ordered, file ="cd4.cds_graph_test_res_subset_ordered.csv")

#The class issue has to be corrected by 
rowData(cd4.cds)$gene_name <- rownames(cd4.cds)
rowData(cd4.cds)$gene_short_name <- rowData(cd4.cds)$gene_name


pdf("top_pseudo_changing_genes_CD4.pdf",width=5,height=4)
plot_cells(cd4.cds, genes = head(cd4.graph_deg_ids),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE,min_expr=0)+ theme(axis.text = element_text(size = 10))
dev.off()


#CD8
cd8.cds_graph_test_res <- graph_test(cd8.cds,
                                     neighbor_graph = "principal_graph",
                                     cores = 6)
head(cd8.cds_graph_test_res)

cd8.cds_graph_test_res_sig <- subset(cd8.cds_graph_test_res, q_value < 0.05)

cd8.cds_graph_test_res_ordered = cd8.cds_graph_test_res_sig[order(cd8.cds_graph_test_res_sig$q_value,(-1*cd8.cds_graph_test_res_sig$morans_I)),]

 contamination_genes = c('Ctrb1', 'Try5', 'Cela2a', 'Gcg', 'Clps', 'Try4', 'Pnlip', 'Cela3b', 'Cpa1', '2210010C04Rik', 'Sycn', 'Ins2', 'Rnase1', 'Reg1', 'Ctrc', 'Sst', 'Cel', 'Cpb1', 'Ins1', 'Pnliprp1', 'Zg16', 'Iapp', 'Ppy', 'Cpa2', 'Pla2g1b', 'Prss2', 'Cela1', 'Ctrl', 'Gm20481', 'Hspa1b', 'Lag3', 'Pdcd1', 'Hspa1a', 'Tnfsf11')

cd8.cds_graph_test_res_ordered=  cd8.cds_graph_test_res_ordered[!rownames(cd8.cds_graph_test_res_ordered)%in%contamination_genes,]
 
#to write files
write.csv(cd8.cds_graph_test_res_ordered, file ="cd8.cds_graph_test_res_ordered.csv")

#The class issue has to be corrected by 
rowData(cd8.cds)$gene_name <- rownames(cd8.cds)
rowData(cd8.cds)$gene_short_name <- rowData(cd8.cds)$gene_name

# plot gene expression as a function of pseudotime
pdf("top_pseudo_changing_genes_CD8.pdf",width=5,height=4)
plot_cells(cd8.cds, genes = head(rownames(cd8.cds_graph_test_res_ordered)),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE,min_expr=0)+ theme(axis.text = element_text(size = 10))
dev.off()


# use plot_pseudotime_heatmap function

## make function to do heatmap

library(ComplexHeatmap)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(circlize)
library(monocle3)

# taken from https://github.com/cole-trapnell-lab/monocle-release/issues/295
make_heatmap = function(cds_temp,gene_list,normalization) {
  if(normalization=="size_only") {
    pt.matrix <- exprs(cds_temp)[match(gene_list,rownames(rowData(cds_temp))),order(pseudotime(cds_temp))]
  } else {
    pt.matrix = normalized_counts(cds_temp, norm_method = "log")[match(gene_list,rownames(rowData(cds_temp))),order(pseudotime(cds_temp))]
  }
  pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
  pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
  rownames(pt.matrix) <- gene_list

#  htkm <- Heatmap(
#  pt.matrix,
#  name                         = "z-score",
#  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
#  show_row_names               = TRUE,
#  show_column_names            = FALSE,
#  row_names_gp                 = gpar(fontsize = 6),
#  km = 3,
#  row_title_rot                = 0,
#  cluster_rows                 = TRUE,
#  cluster_row_slices           = FALSE,
#  cluster_columns              = FALSE)
  
  hthc <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

  #print(htkm)
  print(hthc)
  return(hthc)
}

pdf("top_diff_genes_heatmap_CD8.pdf",width=3,height=1.2)
top_diff_genes_heatmap_CD8 = make_heatmap(cds_temp=cd8.cds,gene_list=head(rownames(cd8.cds_graph_test_res_ordered)),normalization="log")
dev.off()

pdf("top_diff_genes_heatmap_CD4.pdf",width=3,height=1.2)
top_diff_genes_heatmap_CD4 = make_heatmap(cds_temp=cd4.cds,gene_list=head(cd4.graph_deg_ids),normalization="log")
dev.off()

####7. Identify gene module that change as a function of pseudotime ####

make_heatmap_by_module = function(cds_temp,module_df,normalization,topX_genes=NULL,pval_df=NULL) {
  # first make it so I only keep top X gene in each module
  if(!is.null(topX_genes) & !is.null(pval_df)) {
    
    # sort by pvalue
    pval_df_sorted = pval_df[match(module_df$id,rownames(pval_df)),]
    pval_df_sorted = pval_df_sorted %>% select(p_value,q_value,morans_I)
    module_df = cbind(module_df,pval_df_sorted)
    module_df = module_df[order(module_df$p_value,(-1*module_df$morans_I)),]
    #module_df_sorted = module_df[match(rownames(pval_df_sorted),module_df$id),]
    module_df = module_df %>% group_by(module) %>% slice_head(n=topX_genes)
  }
  if(normalization=="size_only") {
    pt.matrix <- exprs(cds_temp)[match(module_df$id,rownames(rowData(cds_temp))),order(pseudotime(cds_temp))]
  } else {
    pt.matrix = normalized_counts(cds_temp, norm_method = "log")[match(module_df$id,rownames(rowData(cds_temp))),order(pseudotime(cds_temp))]
  }
  pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
  pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
  rownames(pt.matrix) <- module_df$id

  
  hthc <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  row_split                    = module_df$module)

  print(hthc)
  return(hthc)
}

#gene modules by cluster
#CD4
cd4.gene_modules <- find_gene_modules(cd4.cds[cd4.graph_deg_ids,],
                                  resolution=1e-3)

seurat_clusters_cd4 = colData(cd4.cds)$seurat_clusters
seurat_clusters_cd4 = as.character(seurat_clusters_cd4)
seurat_clusters_cd4 = as.numeric(seurat_clusters_cd4)
seurat_clusters_cd4 = as.factor(seurat_clusters_cd4)

names(seurat_clusters_cd4) = rownames(colData(cd4.cds))


cd4.cluster_groups <- tibble::tibble(cell=row.names(colData(cd4.cds)),
                                     cluster_group=seurat_clusters_cd4)
cd4.agg_mat <- aggregate_gene_expression(cd4.cds, cd4.gene_modules, cd4.cluster_groups)
row.names(cd4.agg_mat) <- paste0("Module ", row.names(cd4.agg_mat))

#rename_seurat cluster id
rename.cd4.agg_mat = cd4.agg_mat

colnames(rename.cd4.agg_mat) <- c('Lef1 memory','CD29 CD25','Central memory', 'CD25 CTLA4 Foxp3-', 'CD25 Th1','CD25 CTLA4 Foxp3+', 'CTLA4 memory', 'Th17 like','Th17 memory like')



pdf("cd4_module_vs_cluster_heatmap.pdf",width=3,height=3.5)
pheatmap::pheatmap(rename.cd4.agg_mat,
                   scale="column", clustering_method="ward.D2", fontsize_row = 6, fontsize_col =7)
dev.off()

#to write files
write.csv(cd4.gene_modules, file ="cd4.gene_modules.csv")


#CD8

cd8.gene_modules <- find_gene_modules(cd8.cds[rownames(cd8.cds_graph_test_res_ordered),],
                                      resolution=1e-3)


seurat_clusters_cd8 = colData(cd8.cds)$seurat_clusters
seurat_clusters_cd8 = as.character(seurat_clusters_cd8)
seurat_clusters_cd8 = as.numeric(seurat_clusters_cd8)
seurat_clusters_cd8 = as.factor(seurat_clusters_cd8)

names(seurat_clusters_cd8) = rownames(colData(cd8.cds))
  
cd8.cluster_groups <- tibble::tibble(cell=row.names(colData(cd8.cds)),
                                 cluster_group=seurat_clusters_cd8)
cd8.agg_mat <- aggregate_gene_expression(cd8.cds, cd8.gene_modules, cd8.cluster_groups)
row.names(cd8.agg_mat) <- paste0("Module ", row.names(cd8.agg_mat))

#rename_seurat cluster id
rename.cd8.agg_mat = cd8.agg_mat

colnames(rename.cd8.agg_mat) <- c('Ly6a central memory','Effector memory','Central memory', 'CTLA2a Trm')

pdf("cd8_module_vs_cluster_heatmap.pdf",width=3,height=3.5)
pheatmap::pheatmap(rename.cd8.agg_mat,
                   scale="column", clustering_method="ward.D2", fontsize_row = 8, fontsize_col =8)
dev.off()



moduled_heatmap_cd8 = make_heatmap_by_module(cds_temp=cd8.cds,module_df=cd8.gene_modules,normalization="log",topX_genes=3,pval_df=cd8.cds_graph_test_res_ordered)
pdf("moduled_heatmap_cd8.pdf",width=5,height=2.5)
print(moduled_heatmap_cd8)
dev.off()


# make heatmap for cd4
moduled_heatmap_cd4 = make_heatmap_by_module(cds_temp=cd4.cds,module_df=cd4.gene_modules,normalization="log",topX_genes=3,pval_df=cd4.cds_graph_test_res_subset_ordered)
pdf("moduled_heatmap_cd4.pdf",width=5,height=4.5)
print(moduled_heatmap_cd4)
dev.off()