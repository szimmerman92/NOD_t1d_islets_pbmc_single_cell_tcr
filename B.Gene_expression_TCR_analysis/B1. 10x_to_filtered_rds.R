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

# This script takes the output from the 10X Cell Ranger, filters cells,
# clusters data, and outputs Seurat/rds objects

# Define samples
T1D_islets_samples <- c("A001", "A002", "A003", "A004", "A005", "A006")
T1D_PBMC_samples <- c("A007", "A008", "A009", "A010", "A011", "A012")
NonT1D_PBMC_samples <- c("A013", "A014", "A015", "A016")

# Set working directory

setwd("DataProcessing")


##### step-1: Read data ####
# Read the cellranger output from filtered_feature_bc_matrix directory. 
#The data directory has the following files for each sample
# barcodes.tsv.gz
# features.tsv.gz
# matrix.mtx.gz

data <- Read10X(data.dir = "filtered_feature_bc_matrix/filtered_feature_bc_matrix_A006")

#Set the sample name (example for sample A003)
sample <- "A006"



##### step 2: Create a Seurat object ####
data2 <- CreateSeuratObject(counts = data, project = "sample", min.cells = 3, min.features = 200)
data2

##### step 3: Preprocessing and QC of data ####
#The [[ operator can add columns to object metadata.
#Rename a column
idx = which('orig.ident' == colnames(data2@meta.data))
colnames(data2@meta.data)[idx] <- "Sample"

# Save unfiltered RDS
file_name <- paste0(sample, '_unfiltered.rds')
saveRDS(data2, file = file_name)

##QC step 1: Seurat default- Remove low quality cells based on overall mitochondrial gene expression and aberrant high gene count
data2[["percent.mt"]] <- PercentageFeatureSet(data2, pattern = "^mt-")

#remove low quality cells 
#There are some cells with abnormally high count of nCount_RNA so we remove this. 
data2 <- subset(data2, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 20)
data2

#Keep cells that express either Cd3e or Cd3d or Cd3g
data2 <- subset(x = data2, subset = Cd3e > 0 | Cd3d > 0 | Cd3g > 0)
data2


##QC step 2: Filter cells based on expression of housekeeping genes and mitochondrial genes

#load mouse housekeeping genes
housekeeping_genes <- readLines('GeneLists/HK_Satija_mouse.txt')

# Intersection of genes
s1 <- rownames(data2)
s2 <- housekeeping_genes
housekeeping_genes <- s1[ toupper(s1) %in% toupper(s2) ]
subset_counts_matrix <- data.frame(data2@assays$RNA@counts)[housekeeping_genes,]
subset_counts_matrix[subset_counts_matrix > 0] = 1 
summary(colSums(subset_counts_matrix))

counts <- colSums(subset_counts_matrix)
counts2 <- data.frame(counts)
counts3 <- counts2 %>% tibble::rownames_to_column() %>% filter(counts2 > length(housekeeping_genes)/2)

all.genes <- rownames(x = data2)



##QC step 3: filtering based on list of mitochondrial genes
#load mouse mitochondrial genes
mito_genes <- readLines('GeneLists/mito_genes_mouse.txt', warn = FALSE)
# Intersection of mito genes
# We removed cells if the mito gene expression was higher than 2 standard deviations from the mean
#We took the mean and the SD of the number of mitochondrial genes present in each cell. Then to get the cutoff (mean + (2* SD))
# The number of counts to be filtered is determined on a sample-to-sample basis

s1 <- rownames(data2)
s2 <- mito_genes
mito_genes <- s1[ toupper(s1) %in% toupper(s2) ]
subset_counts_matrix <- data.frame(data2@assays$RNA@counts)[mito_genes,]
subset_counts_matrix[subset_counts_matrix > 0] =1 
summary(colSums(subset_counts_matrix))

counts4 <- colSums(subset_counts_matrix)

mean_prev_mito_genes = mean(counts4)
sd_prev_mito_genes = sd(counts4)
mito_prev_cutoff = mean_prev_mito_genes + (2*sd_prev_mito_genes)

counts5 <- data.frame(counts4)
counts6 <- counts5 %>% tibble::rownames_to_column() %>% filter(counts4 < mito_prev_cutoff) # Change this number

cells <- intersect(counts6$rowname,counts3$rowname)


##### step 4: Normalize and scale data ####
working_so <- NormalizeData(data2)
working_so <- subset(working_so, cells = gsub('.', '-', cells, fixed = TRUE))
working_so <- ScaleData(object = working_so, features = rownames(working_so))

working_so



##### step 5: Find variable genes ####

####Identification of highly variable features (feature selection) for T1D_islets ####
working_so1 <- FindVariableFeatures(object = working_so, mean.function = ExpMean, dispersion.function = LogVMR)
hv.genes <- head(x = VariableFeatures(object = working_so1), 2000)



##### step 6: Run PCA, find neighbors, clustering, and UMAP####
working_so1 <- RunPCA(object = working_so1, pc.genes = hv.genes, pcs.print = 1:5,
                     genes.print = 5, pcs.compute = 50)
working_so1 <- FindNeighbors(object = working_so1, dims = 1:30)
working_so1 <- FindClusters(object = working_so1, resolution = 0.5)
working_so1 <- RunUMAP(object = working_so1, reduction = "pca", dims = 1:15, n_neighbors = 15, min_dist = 0.3)



# Plot the UMAP
DimPlot(object = working_so1, reduction = 'umap')

# Add the UMAP coordinates to the metadata
working_so1@meta.data$UMAP_1 <- working_so1@reductions$umap@cell.embeddings[, 1]
working_so1@meta.data$UMAP_2 <- working_so1@reductions$umap@cell.embeddings[, 2]

# Optional, save rds object after mito/HK filtering
file_name <- paste0(sample, '_postMitoHK.rds')
saveRDS(working_so1, file = file_name)

#optionally read in existing Seurat object
#working_so1 <- readRDS("A003_postMitoHK.rds")





##### step 7: Add TCR metadata from the clone pipeline####
# 
# These include all cells, including unfiltered cells
# Add TCR metadata from the clone pipeline
# These include all cells, including unfiltered cells
TCR_data <- read.csv('sam.clone_v2/alpha_beta_TCR/islet_alpha_beta_TCR_matches_pool_6.csv')

# Remove duplicates (identical rows)
TCR_data <- unique(TCR_data)

# Add rownames to TCR data
rownames(TCR_data) <- TCR_data$Barcode
# remove barcode column
TCR_data = TCR_data[,-1]


# merge the datasets together only keeping cells that are in working_so1@meta.data
# puts any NAs for values in working_so1@meta.data not in TCR_data
meta = merge(working_so1@meta.data,TCR_data,by="row.names",all.x=TRUE)

# reset rownames to be barcodes
rownames(meta) = meta[,1]
meta = meta[,-1]
working_so1@meta.data = meta

# replace NAs with notcr if TCR missing and set frequency to 0
working_so1@meta.data$Matching_pre_filter[is.na(working_so1@meta.data$Matching_pre_filter)] = "notcr"
working_so1@meta.data$Frequency_pre_filter[is.na(working_so1@meta.data$Frequency_pre_filter)] = 0
working_so1@meta.data$TCR[is.na(working_so1@meta.data$TCR)] = "notcr"

# Rename a column
idx = which('seurat_clusters' == colnames(working_so1@meta.data))
colnames(working_so1@meta.data)[idx] <- "Seurat_clusters"


# Save data
mtx_path <- paste0(sample, '_filtered_mtx.csv')
meta_path <- paste0(sample, '_filtered_meta.csv')
rds_path <- paste0(sample, '_filtered.rds')

write.csv(working_so1@meta.data, meta_path, quote = FALSE)
write.csv(working_so1@assays$RNA@counts, mtx_path, quote = FALSE)
saveRDS(working_so1, file = rds_path)


