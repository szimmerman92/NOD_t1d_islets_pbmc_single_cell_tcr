library(Seurat)

# Define samples
samples <- c("A007_pbmc", "A008_pbmc", "A009_pbmc", "A010_pbmc", "A011_pbmc", "A012_pbmc",
             "A001_islets", "A002_islets", "A003_islets", "A004_islets", "A005_islets", "A006_islets",
             "A013_pbmc", "A014_pbmc", "A015_pbmc", "A016_pbmc")

# Use this loop to load rds objects
for(i in 1:length(samples)){
  sample <- samples[i]
  
  # Load file name
  file_name <- paste0('RDSobjects/', sample, '.rds')
  
  # Load rds object
  rds <- readRDS(file_name)
  
  # Assign rds object to sample
  assign(sample, rds)
}



# In each sample, replace 'matching' with 'sample_matching'
# For example, in the sample A001, the level 'matching' in the
# 'Matching' column will be replaced with 'A001_matching' that means this is A001 islets matching
# This is done for convenience, since all the data will be integrated

seurat_objects <- c(A007_pbmc, A008_pbmc, A009_pbmc, A010_pbmc, A011_pbmc, A012_pbmc,
                    A001_islets, A002_islets, A003_islets, A004_islets, A005_islets, A006_islets,
                    A013_pbmc, A014_pbmc, A015_pbmc, A016_pbmc)



library(plyr)
for(i in 1:length(samples)){
  
  sample <- samples[i]
  rds <- seurat_objects[[i]]
  
  # Rename data
  rds$Matching <- revalue(rds$Matching, c("matching" = paste0(sample, "_matching")))
  
  # Assign rds object to sample
  assign(sample, rds)
  
}


#Add sample id as meta data before integrating 
A001_islets$sampleid <- "A001_islets"
A002_islets$sampleid <- "A002_islets"
A003_islets$sampleid <- "A003_islets"
A004_islets$sampleid <- "A004_islets"
A005_islets$sampleid <- "A005_islets"
A006_islets$sampleid <- "A006_islets"
A007_pbmc$sampleid <- "A007_t1Dpbmc"
A008_pbmc$sampleid <- "A008_t1Dpbmc"
A009_pbmc$sampleid <- "A009_t1Dpbmc"
A010_pbmc$sampleid <- "A010_t1Dpbmc"
A011_pbmc$sampleid <- "A011_t1Dpbmc"
A012_pbmc$sampleid <- "A012_t1Dpbmc"
A013_pbmc$sampleid <- "A013_non-t1Dpbmc"
A014_pbmc$sampleid <- "A014_non-t1Dpbmc"
A015_pbmc$sampleid <- "A015_non-t1Dpbmc"
A016_pbmc$sampleid <- "A016_non-t1Dpbmc"

#add sample type before integrating
A001_islets$sampletype <- "t1dislets"
A002_islets$sampletype <- "t1dislets"
A003_islets$sampletype <- "t1dislets"
A004_islets$sampletype <- "t1dislets"
A005_islets$sampletype <- "t1dislets"
A006_islets$sampletype <- "t1dislets"
A007_pbmc$sampletype <- "t1dpbmc"
A008_pbmc$sampletype <- "t1dpbmc"
A009_pbmc$sampletype <- "t1dpbmc"
A010_pbmc$sampletype <- "t1dpbmc"
A011_pbmc$sampletype <- "t1dpbmc"
A012_pbmc$sampletype <- "t1dpbmc"
A013_pbmc$sampletype <- "nont1dpbmc"
A014_pbmc$sampletype <- "nont1dpbmc"
A014_pbmc$sampletype <- "nont1dpbmc"
A014_pbmc$sampletype <- "nont1dpbmc"

#############
#Integrate PBMC data#
##################

#This is very time consuming step. so we run it on server. 
samples.list <- c(A007_pbmc, A008_pbmc, A009_pbmc, A010_pbmc, A011_pbmc, A012_pbmc)

# First run SCTransform
for(i in 1:length(samples.list)){
  samples.list[[i]] <- SCTransform(samples.list[[i]], verbose = T)
}

# To avoid errors,
#Sys.setenv('R_MAX_VSIZE'= 100000000000)
#options(future.globals.maxSize = 3145728000)
# Note: 3000 * 1024^2 bytes to resolve error with 2GB list

# Keep top 3000 genes for integration
samples.features <- SelectIntegrationFeatures(object.list = samples.list, nfeatures = 3000)
samples.list <- PrepSCTIntegration(object.list = samples.list,
                                   anchor.features = samples.features,
                                   verbose = T)

# Find the integration anchors
samples.anchors <- FindIntegrationAnchors(object.list = samples.list,
                                          normalization.method = "SCT", 
                                          anchor.features = samples.features,
                                          verbose = T)

#Optional
#save anchor
# Save RDS
saveRDS(samples.anchors, 'Integrated_rds/samples.anchors_t1dpbmc.rds')

#optional load
samples.anchors <- readRDS("Integrated_rds/samples.anchors_t1dpbmc.rds")

# Integrate the data
samples.integrated <- IntegrateData(anchorset = samples.anchors,
                                    normalization.method = "SCT", 
                                    verbose = T)

# Change the default assay
DefaultAssay(object = samples.integrated) <- "RNA"

# Scale data, find variable genes, run PCA
samples.integrated <- ScaleData(samples.integrated  , verbose = T)
samples.integrated <- FindVariableFeatures(object = samples.integrated, mean.function = ExpMean,                                                 dispersion.function = LogVMR)
#samples.integrated <- FindVariableFeatures(samples.integrated)
samples.integrated <- RunPCA(samples.integrated, npcs = 30, verbose = T)

# Find nearest neighbors, run clustering, run UMAP
# Change resolution to get about 7 clusters
samples.integrated <- FindNeighbors(object = samples.integrated , dims = 1:16)
samples.integrated <- FindClusters(samples.integrated , resolution = 0.5, verbose=T)
samples.integrated <- RunUMAP(samples.integrated , reduction = "pca", dims = 1:16)

# Plot the UMAP
DimPlot(samples.integrated, reduction = "umap")

# Save RDS
saveRDS(samples.integrated, 'Integrated_rds/Integratedt1dpbmc.rds')




##################
# Integrate islets data
##################
samples.list <- c(A001_islets, A002_islets, A003_islets, A004_islets, A005_islets, A006_islets)

# First run SCTransform
for(i in 1:length(samples.list)){
  samples.list[[i]] <- SCTransform(samples.list[[i]], verbose = T)
}

# To avoid errors,
#Sys.setenv('R_MAX_VSIZE'= 100000000000)
#options(future.globals.maxSize = 3145728000)
# Note: 3000 * 1024^2 bytes to resolve error with 2GB list

# Keep top 3000 genes for integration
samples.features <- SelectIntegrationFeatures(object.list = samples.list, nfeatures = 3000)
samples.list <- PrepSCTIntegration(object.list = samples.list,
                                   anchor.features = samples.features,
                                   verbose = T)

# Find the integration anchors
samples.anchors <- FindIntegrationAnchors(object.list = samples.list,
                                          normalization.method = "SCT", 
                                          anchor.features = samples.features,
                                          verbose = T)

#Optional
#save anchor
# Save RDS
saveRDS(samples.anchors, 'Integrated_rds/samples.anchors_t1dislets.rds')

#optional load
samples.anchors <- readRDS("Integrated_rds/samples.anchors_t1dislets.rds")



# Integrate the data
samples.integrated <- IntegrateData(anchorset = samples.anchors,
                                    normalization.method = "SCT", 
                                    verbose = T)

# Change the default assay
DefaultAssay(object = samples.integrated) <- "RNA"

# Scale data, find variable genes, run PCA
samples.integrated <- ScaleData(samples.integrated  , verbose = T)
samples.integrated <- FindVariableFeatures(object = samples.integrated, mean.function = ExpMean,                                                 dispersion.function = LogVMR)
#samples.integrated <- FindVariableFeatures(samples.integrated)
samples.integrated <- RunPCA(samples.integrated, npcs = 30, verbose = T)

# Find nearest neighbors, run clustering, run UMAP
# Change resolution to get about 7 clusters
samples.integrated <- FindNeighbors(object = samples.integrated , dims = 1:16)
samples.integrated <- FindClusters(samples.integrated , resolution = 0.5, verbose=T)
samples.integrated <- RunUMAP(samples.integrated , reduction = "pca", dims = 1:16)

# Plot the UMAP
DimPlot(samples.integrated, reduction = "umap")

# Save RDS
saveRDS(samples.integrated, 'Integrated_rds/Integratedislets.rds')



##################
# Integrate non-t1dpbmc data
##################

# Define samples
samples.list <- c(A013_pbmc, A014_pbmc, A015_pbmc, A016_pbmc)

# First run SCTransform
for(i in 1:length(samples.list)){
  samples.list[[i]] <- SCTransform(samples.list[[i]], verbose = T)
}

# To avoid errors,
#Sys.setenv('R_MAX_VSIZE'= 100000000000)
#options(future.globals.maxSize = 3145728000)
# Note: 3000 * 1024^2 bytes to resolve error with 2GB list

# Keep top 3000 genes for integration
samples.features <- SelectIntegrationFeatures(object.list = samples.list, nfeatures = 3000)
samples.list <- PrepSCTIntegration(object.list = samples.list,
                                   anchor.features = samples.features,
                                   verbose = T)

# Find the integration anchors
samples.anchors <- FindIntegrationAnchors(object.list = samples.list,
                                          normalization.method = "SCT", 
                                          anchor.features = samples.features,
                                          verbose = T)

#Optional
#save anchor
# Save RDS
saveRDS(samples.anchors, 'Integrated_rds/samples.anchors_nont1dpbmc.rds')

#optional load
#samples.anchors <- readRDS("Integrated_rds/samples.anchors_nont1dpbmc.rds")



# Integrate the data
samples.integrated <- IntegrateData(anchorset = samples.anchors,
                                    normalization.method = "SCT", 
                                    verbose = T)

# Change the default assay
DefaultAssay(object = samples.integrated) <- "RNA"

# Scale data, find variable genes, run PCA
samples.integrated <- ScaleData(samples.integrated  , verbose = T)
samples.integrated <- FindVariableFeatures(object = samples.integrated, mean.function = ExpMean,                                                 dispersion.function = LogVMR)
#samples.integrated <- FindVariableFeatures(samples.integrated)
samples.integrated <- RunPCA(samples.integrated, npcs = 30, verbose = T)

# Find nearest neighbors, run clustering, run UMAP
# Change resolution to get about 7 clusters
samples.integrated <- FindNeighbors(object = samples.integrated , dims = 1:16)
samples.integrated <- FindClusters(samples.integrated , resolution = 0.5, verbose=T)
samples.integrated <- RunUMAP(samples.integrated , reduction = "pca", dims = 1:16)

# Plot the UMAP
DimPlot(samples.integrated, reduction = "umap")

# Save RDS
saveRDS(samples.integrated, 'Integrated_rds/Integratednont1dpbmc.rds')



##################
# Integrate all PBMC (T1D non-t1D) samples together = integratedt1dnont1dpbmc
#Integrate all sample (n=16) = integratedallsamples
##################

#we first will merge samples 

#merging all 6 islets samples 

#we followed merge nornalized data as our data is already normalized 
#https://satijalab.org/seurat/articles/merge_vignette.html
#https://github.com/satijalab/seurat/issues/2676

####Merging T1D islets object (A001+A002+A003+A004A005+A006) ####
t1disletsMerged <- merge(A001_islets, y = c(A002_islets, A003_islets, A004_islets, A005_islets, A006_islets), add.cell.ids = c("A001", "A002", "A003", "A004", "A005", "A006"), project = "t1disletsMerged", merge.data = TRUE)

t1disletsMerged

t1disletsMerged$mergedsampletype <- "T1D_Islets"

#Save unfiltered_merged seurat object
# Save RDS
saveRDS(t1disletsMerged, 'Integrated_rds/t1disletsMerged.rds')

#optional load
t1disletsMerged <- readRDS("Integrated_rds/t1disletsMerged.rds")



####Merging T1D pbmc object (A007+A008+A009+A010+A011+A012) ####
t1dpbmcMerged <- merge(A007_pbmc, y = c(A008_pbmc, A009_pbmc, A010_pbmc, A011_pbmc, A012_pbmc), add.cell.ids = c("A007", "A008", "A009", "A010", "A011", "A012"), project = "t1dpbmcMerged", merge.data = TRUE)

t1dpbmcMerged

t1dpbmcMerged$mergedsampletype <- "T1D_pbmc"

#Save unfiltered_merged seurat object
# Save RDS
saveRDS(t1dpbmcMerged, 'Integrated_rds/t1dpbmcMerged.rds')

#optional load
t1dpbmcMerged <- readRDS("Integrated_rds/t1dpbmcMerged.rds")



####Merging npnT1D pbmc object (A013+A014+A015+A016) ####
nont1dpbmcMerged <- merge(A013_pbmc, y = c(A014_pbmc, A015_pbmc, A016_pbmc), add.cell.ids = c("A013", "A014", "A015", "A016"), project = "nont1dpbmcMerged", merge.data = TRUE)

nont1dpbmcMerged

nont1dpbmcMerged$mergedsampletype <- "nonT1D_pbmc"

#Save unfiltered_merged seurat object
# Save RDS
saveRDS(nont1dpbmcMerged, 'Integrated_rds/nont1dpbmcMerged.rds')

#optional load
nont1dpbmcMerged <- readRDS("Integrated_rds/nont1dpbmcMerged.rds")



####Now start integreting all these three large merged object

# Define samples
samples.list <- c(t1disletsMerged, t1dpbmcMerged, nont1dpbmcMerged)

# First run SCTransform
for(i in 1:length(samples.list)){
  samples.list[[i]] <- SCTransform(samples.list[[i]], verbose = T)
}


# Keep top 3000 genes for integration
samples.features <- SelectIntegrationFeatures(object.list = samples.list, nfeatures = 3000)
samples.list <- PrepSCTIntegration(object.list = samples.list,
                                   anchor.features = samples.features,
                                   verbose = T)

# Find the integration anchors
samples.anchors <- FindIntegrationAnchors(object.list = samples.list,
                                          normalization.method = "SCT", 
                                          anchor.features = samples.features,
                                          verbose = T)
#Optional
#save anchor
# Save RDS
saveRDS(samples.anchors, 'Integrated_rds/samples.anchors_allmergedsamples.rds')

#optional load
samples.anchors <- readRDS("Integrated_rds/samples.anchors_allmergedsamples.rds")



# Integrate the data
samples.integrated <- IntegrateData(anchorset = samples.anchors,
                                    normalization.method = "SCT", 
                                    verbose = T)

# Change the default assay
DefaultAssay(object = samples.integrated) <- "RNA"

# Scale data, find variable genes, run PCA
samples.integrated <- ScaleData(samples.integrated  , verbose = T)
samples.integrated <- FindVariableFeatures(object = samples.integrated, mean.function = ExpMean,                                                 dispersion.function = LogVMR)
#samples.integrated <- FindVariableFeatures(samples.integrated)
samples.integrated <- RunPCA(samples.integrated, npcs = 30, verbose = T)

# Find nearest neighbors, run clustering, run UMAP
# Change resolution to get about 7 clusters
samples.integrated <- FindNeighbors(object = samples.integrated , dims = 1:18)
samples.integrated <- FindClusters(samples.integrated , resolution = 0.5, verbose=T)
samples.integrated <- RunUMAP(samples.integrated , reduction = "pca", dims = 1:18)

# Plot the UMAP
#DimPlot(samples.integrated, reduction = "umap")

# Save RDS
saveRDS(samples.integrated, 'Integrated_rds/Integrated_allmergedsamples.rds')

#end




##################
# Integrate all 16 samples together
##################

# Define samples
samples.list <- c(A007_pbmc, A008_pbmc, A009_pbmc, A010_pbmc, A011_pbmc, A012_pbmc,
                    A001_islets, A002_islets, A003_islets, A004_islets, A005_islets, A006_islets,
                    A013_pbmc, A014_pbmc, A015_pbmc, A016_pbmc)

# First run SCTransform
for(i in 1:length(samples.list)){
  samples.list[[i]] <- SCTransform(samples.list[[i]], verbose = T)
}


# Keep top 3000 genes for integration
samples.features <- SelectIntegrationFeatures(object.list = samples.list, nfeatures = 3000)
samples.list <- PrepSCTIntegration(object.list = samples.list,
                                   anchor.features = samples.features,
                                   verbose = T)

# Find the integration anchors
samples.anchors <- FindIntegrationAnchors(object.list = samples.list,
                                          normalization.method = "SCT", 
                                          anchor.features = samples.features,
                                          verbose = T)

#Optional
#save anchor
# Save RDS
saveRDS(samples.anchors, 'Integrated_rds/samples.anchors_allsamples.rds')

#optional load
samples.anchors <- readRDS("Integrated_rds/samples.anchors_allsamples.rds")



# Integrate the data
samples.integrated <- IntegrateData(anchorset = samples.anchors,
                                    normalization.method = "SCT", 
                                    verbose = T)

# Change the default assay
DefaultAssay(object = samples.integrated) <- "RNA"

# Scale data, find variable genes, run PCA
samples.integrated <- ScaleData(samples.integrated  , verbose = T)
samples.integrated <- FindVariableFeatures(object = samples.integrated, mean.function = ExpMean,                                                 dispersion.function = LogVMR)
#samples.integrated <- FindVariableFeatures(samples.integrated)
samples.integrated <- RunPCA(samples.integrated, npcs = 30, verbose = T)

# Find nearest neighbors, run clustering, run UMAP
# Change resolution to get about 7 clusters
samples.integrated <- FindNeighbors(object = samples.integrated , dims = 1:16)
samples.integrated <- FindClusters(samples.integrated , resolution = 0.5, verbose=T)
samples.integrated <- RunUMAP(samples.integrated , reduction = "pca", dims = 1:16)

# Plot the UMAP
#DimPlot(samples.integrated, reduction = "umap")

# Save RDS
saveRDS(samples.integrated, 'Integrated_rds/Integrated_allmergedsamples.rds')
#optional load
samples.integrated <- readRDS("Integrated_rds/Integrated_allmergedsamples.rds")
#end


