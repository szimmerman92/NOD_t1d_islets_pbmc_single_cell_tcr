#!/usr/bin/env Rscript
library(Seurat)
library(dplyr)
library(statmod)
T1D_pbmcs = readRDS("/n/data1/joslin/icrb/kostic/zahir/tcell/seurat/DataProcessing/Integrated_rds/Integratedpbmc.rds")
T1D_islets = readRDS("/n/data1/joslin/icrb/kostic/zahir/tcell/seurat/DataProcessing/Integrated_rds/Integratedislets.rds")
# subset cells to only include matching ones
T1D_pbmcs_matching <- subset(T1D_pbmcs, subset = Matching == "A007_pbmc_matching" | Matching == "A008_pbmc_matching" | Matching == "A009_pbmc_matching" | Matching == "A010_pbmc_matching" | Matching == "A011_pbmc_matching" | Matching == "A012_pbmc_matching")
T1D_islets_matching <- subset(T1D_islets, subset = Matching == "A001_islets_matching" | Matching == "A002_islets_matching" | Matching == "A003_islets_matching" | Matching == "A004_islets_matching" | Matching == "A005_islets_matching" | Matching == "A006_islets_matching")
# combine cells based on their shared clonotype
T1D_pbmcs_clonotypes = unique(T1D_pbmcs_matching@meta.data$TCR)
T1D_pbmcs_raw_counts_integrated_df = T1D_pbmcs_matching@assays$RNA@counts
T1D_pbmcs_raw_counts_integrated_df = as.data.frame(T1D_pbmcs_raw_counts_integrated_df)
# make clonotype by gene matrix for PBMCs 
T1D_pbmcs_matching_bulk_expression = sapply(T1D_pbmcs_clonotypes, function(clonotype) {
  cells_same_clonotype = rownames(T1D_pbmcs_matching@meta.data[T1D_pbmcs_matching@meta.data$TCR == clonotype,])
  expression_with_clonotype = T1D_pbmcs_raw_counts_integrated_df[,cells_same_clonotype,drop=FALSE]
  bulk_expression = rowSums(expression_with_clonotype)
  return(bulk_expression)
})

T1D_islet_clonotypes = unique(T1D_islets_matching@meta.data$TCR)
T1D_islets_raw_counts_integrated_df = T1D_islets_matching@assays$RNA@counts
T1D_islets_raw_counts_integrated_df = as.data.frame(T1D_islets_raw_counts_integrated_df)

# make clonotype by gene matrix for islets 
T1D_iselts_matching_bulk_expression = sapply(T1D_islet_clonotypes, function(clonotype) {
  cells_same_clonotype = rownames(T1D_islets_matching@meta.data[T1D_islets_matching@meta.data$TCR == clonotype,])
  expression_with_clonotype = T1D_islets_raw_counts_integrated_df[,cells_same_clonotype,drop=FALSE]
  bulk_expression = rowSums(expression_with_clonotype)
  return(bulk_expression)
})

# create metadata data frame
bulk_metadata = data.frame(clonotype=c(colnames(T1D_pbmcs_matching_bulk_expression),colnames(T1D_iselts_matching_bulk_expression)),tissue=c(rep("PBMC",ncol(T1D_pbmcs_matching_bulk_expression)),rep("islet",ncol(T1D_iselts_matching_bulk_expression))))
# now combine the 2 dataframe
all_bulk_expression_data = merge(T1D_pbmcs_matching_bulk_expression,T1D_iselts_matching_bulk_expression,by="row.names",all=TRUE)
# replace NA with 0. 
all_bulk_expression_data[is.na(all_bulk_expression_data)] <- 0
rownames(all_bulk_expression_data) = all_bulk_expression_data[,1]
all_bulk_expression_data = all_bulk_expression_data[,-1]
all_bulk_expression_data = as.matrix(all_bulk_expression_data)
# check to make sure colnames are correct order
all.equal(sapply(strsplit(colnames(all_bulk_expression_data),split="[.]"), function(x) x[1]),bulk_metadata[,1])

# now time to run edgeR
library(edgeR)
# Create a DGEList object
dgList <- DGEList(counts=all_bulk_expression_data, genes = row.names(all_bulk_expression_data))
# Filtering genes
countsPerMillion <- cpm(dgList)
countCheck <- countsPerMillion > 1
keep <- which(rowSums(countCheck) >= 2)
dgList <- dgList[keep,]
#summary(cpm(dgList))
# Normalization
dgList <- calcNormFactors(dgList, method="TMM")
# Set up design matrix
Clone <- factor(bulk_metadata[,1])
Tissue <- factor(bulk_metadata[,2],levels=c("PBMC","islet"))
designMat <- model.matrix(~Clone+Tissue)
row.names(designMat) <- colnames(dgList)
# Estimating Dispersion
dgList <- estimateDisp(dgList, designMat, robust = TRUE)
saveRDS(dgList,file="matching_bulk_edgeR_object.rds")

# DE genes
fit <- glmFit(dgList, designMat)
lrt <- glmLRT(fit)
topTags(lrt)
# Checks
o <- order(lrt$table$PValue)
cpm(dgList)[o[1:10],]

summary(decideTests(lrt))

pdf("MD_plot.pdf")
plotMD(lrt)
abline(h=c(-1, 1), col="blue")
dev.off()
# Save results
edgeR_result <- topTags(lrt, n=length(row.names(lrt$table)))
DEgenes <- edgeR_result$table[edgeR_result$table$FDR < 0.05,]
all_DEgenes <- edgeR_result$table
DEgenes

write.csv(x=DEgenes, file='MousePseudoBulk_DEgenes_unsorted.csv', row.names = F)
write.csv(x=all_DEgenes, file='MousePseudoBulk_all_DEgenes_unsorted.csv', row.names = F)