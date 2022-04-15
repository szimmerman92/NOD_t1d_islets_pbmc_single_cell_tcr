library(Seurat)

# write metadata and gene expression data to separate csv files for integrated PBMCs
T1D_pbmcs = readRDS("/n/data1/joslin/icrb/kostic/zahir/tcell/seurat/DataProcessing/Integrated_rds/Integratedpbmc.rds")

# get metadata
T1D_pbmcs_metadata = T1D_pbmcs@meta.data
T1D_pbmcs_metadata = cbind(barcode=rownames(T1D_pbmcs_metadata),T1D_pbmcs_metadata)

# extract raw counts
T1D_pbmcs_RNA_raw_counts = T1D_pbmcs@assays$RNA@counts
T1D_pbmcs_RNA_raw_counts = as.matrix(T1D_pbmcs_RNA_raw_counts)
T1D_pbmcs_RNA_raw_counts = cbind(geneName=rownames(T1D_pbmcs_RNA_raw_counts),T1D_pbmcs_RNA_raw_counts)
# extract log normalized counts
T1D_pbmcs_logNorm = T1D_pbmcs@assays$integrated@data
T1D_pbmcs_logNorm = as.matrix(T1D_pbmcs_logNorm)
T1D_pbmcs_logNorm = cbind(geneName=rownames(T1D_pbmcs_logNorm),T1D_pbmcs_logNorm)
# extract log normalized and scalled counts
T1D_pbmcs_logNorm_scaled = T1D_pbmcs@assays$integrated@scale.data
T1D_pbmcs_logNorm_scaled = cbind(geneName=rownames(T1D_pbmcs_logNorm_scaled),T1D_pbmcs_logNorm_scaled)

# write files
write.csv(T1D_pbmcs_metadata,file="/n/data1/joslin/icrb/kostic/zahir/tcell/seurat/DataProcessing/Integrated_rds/Integratedpbmc_metadata.csv",row.names=F,quote=T)
write.csv(T1D_pbmcs_logNorm,file="/n/data1/joslin/icrb/kostic/zahir/tcell/seurat/DataProcessing/Integrated_rds/Integratedpbmc_logNorm.csv",row.names=F,quote=T)
write.csv(T1D_pbmcs_logNorm_scaled,file="/n/data1/joslin/icrb/kostic/zahir/tcell/seurat/DataProcessing/Integrated_rds/Integratedpbmc_logNorm_scaled.csv",row.names=F,quote=T)
write.csv(T1D_pbmcs_RNA_raw_counts,file="/n/data1/joslin/icrb/kostic/zahir/tcell/seurat/DataProcessing/Integrated_rds/Integratedpbmc_RNA_rawcounts.csv",row.names=F,quote=T)

# get cluster 5 and cluster 3. cluster 5 are CD8 effector memory CTLs and cluster 3 is CD4 effector memory T cells

# get cluster 5 cells
T1D_pbmcs_cluster5 = subset(x = T1D_pbmcs, idents = "5")
# get cluster 3 cells
T1D_pbmcs_cluster3 = subset(x = T1D_pbmcs, idents = "3")

# get metadata
T1D_pbmcs_cluster5_metadata = T1D_pbmcs_cluster5@meta.data
T1D_pbmcs_cluster5_metadata = cbind(barcode=rownames(T1D_pbmcs_cluster5_metadata),T1D_pbmcs_cluster5_metadata)
# get raw counts
T1D_pbmcs_cluster5_raw_counts = T1D_pbmcs_cluster5@assays$RNA@counts
T1D_pbmcs_cluster5_raw_counts = as.matrix(T1D_pbmcs_cluster5_raw_counts)
T1D_pbmcs_cluster5_raw_counts = cbind(geneName=rownames(T1D_pbmcs_cluster5_raw_counts),T1D_pbmcs_cluster5_raw_counts)
# get metadata
T1D_pbmcs_cluster3_metadata = T1D_pbmcs_cluster3@meta.data
T1D_pbmcs_cluster3_metadata = cbind(barcode=rownames(T1D_pbmcs_cluster3_metadata),T1D_pbmcs_cluster3_metadata)
# get raw counts
T1D_pbmcs_cluster3_raw_counts = T1D_pbmcs_cluster3@assays$RNA@counts
T1D_pbmcs_cluster3_raw_counts = as.matrix(T1D_pbmcs_cluster3_raw_counts)
T1D_pbmcs_cluster3_raw_counts = cbind(geneName=rownames(T1D_pbmcs_cluster3_raw_counts),T1D_pbmcs_cluster3_raw_counts)
# write files
write.csv(T1D_pbmcs_cluster5_metadata,file="/n/data1/joslin/icrb/kostic/zahir/tcell/seurat/DataProcessing/Integrated_rds/Integratedpbmc_cluster5_metadata.csv",row.names=F,quote=T)
write.csv(T1D_pbmcs_cluster5_raw_counts,file="/n/data1/joslin/icrb/kostic/zahir/tcell/seurat/DataProcessing/Integrated_rds/Integratedpbmc_cluster5_RNA_rawcounts.csv",row.names=F,quote=F)

write.csv(T1D_pbmcs_cluster3_metadata,file="/n/data1/joslin/icrb/kostic/zahir/tcell/seurat/DataProcessing/Integrated_rds/Integratedpbmc_cluster3_metadata.csv",row.names=F,quote=T)
write.csv(T1D_pbmcs_cluster3_raw_counts,file="/n/data1/joslin/icrb/kostic/zahir/tcell/seurat/DataProcessing/Integrated_rds/Integratedpbmc_cluster3_RNA_rawcounts.csv",row.names=F,quote=F)
