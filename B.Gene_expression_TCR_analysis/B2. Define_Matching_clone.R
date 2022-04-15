# This file defines matching status and counts clone size between
# blood and islets samples using the filtered rds objects.
# The additional restriction is that matching status can only be
# applied to AB cells (or 'alpha/beta cells'), which are cells
# that have at least one alpha AND one beta chain in the TCR.
# All the TCRs are written in the form
#
#  alpha_chain|beta_chain
#
# For example, the TCR sequence CALSDRGGSNAKLTF|CASTNSAETLYF belongs
# to an AB cell, but the TCR sequence |CASSDGGAYAEQFF is missing an
# alpha chain, and would therefore not be considered an AB cell.
# Additionally, the clone size of all non-AB cells are defined to be 0.

# Load libraries
#install.packages("hash")
library(hash)
library(plyr)
library(dplyr)


############################################################
#                                                          #
#    Create function to enforce AB matching criteria       #
#                                                          #
############################################################


# Create function to enforce AB matching criteria
AB_Matching <- function(pbmc_rds, islets_rds){
  
  # pbmc = pbmc Seurat/rds object
  # islets = islets Seurat/rds object
  
  # Get just the metadata
  pbmc_df <- pbmc_rds@meta.data
  islets_df <- islets_rds@meta.data
  
  # Turn TCR into character values instead of factors
  pbmc_df$TCR <- as.character(pbmc_df$TCR)
  islets_df$TCR <- as.character(islets_df$TCR)
  
  # Replace NA with 'notcr' in the 'TCR' column
  pbmc_df$TCR[ is.na(pbmc_df$TCR) ] <- 'notcr'
  islets_df$TCR[ is.na(islets_df$TCR) ] <- 'notcr'
  
  # Get a vector of all the unique TCRs for PBMC and islets
  Unique_pbmc_TCR <- unique(pbmc_df$TCR)
  Unique_islets_TCR <- unique(islets_df$TCR)
  
  # From the vector, remove 'notcr'
  Unique_pbmc_TCR <- Unique_pbmc_TCR[ Unique_pbmc_TCR != 'notcr' ]
  Unique_islets_TCR <- Unique_islets_TCR[ Unique_islets_TCR != 'notcr' ]

  
  # Count clone sizes
  # Non-AB cells will be dealt with later
  pbmc_sizes <- dplyr::count(pbmc_df, TCR, sort = TRUE)
  islets_sizes <- dplyr::count(islets_df, TCR, sort = TRUE)
  
  # Now get a vector of all the TCRs in PBMC and islets
  all_pbmc_TCRs <- pbmc_df$TCR
  all_islets_TCRs <- islets_df$TCR
  
  # Run matching analysis on PBMC
  # Initialize empty vectors
  matching <- character(length = nrow(pbmc_df))
  clone_size <- rep(0, nrow(pbmc_df))
  
  for(i in 1:length(all_pbmc_TCRs)){
    # For each TCR
    TCR <- all_pbmc_TCRs[i]
    
    # If it has no sequenced TCR, call it 'notcr' and give it a clone size of 0
    if(TCR == 'notcr'){
      matching[i] <- 'notcr'
      clone_size[i] <- 0
    }
    
    # If it's missing an alpha chain, call it 'beta_chain_only', give it a
    # clone size of 0
    else if(substr(TCR, 1, 1) == "|"){                   # Check first character of TCR
      matching[i] <- "beta_chain_only"
      clone_size[i] <- 0
    }
    
    # If it's missing a beta chain, call it 'alpha_chain_only', give it a
    # clone size of 0
    else if(substr(TCR, nchar(TCR), nchar(TCR)) == "|"){ # Check last character of TCR
      matching[i] <- "alpha_chain_only"
      clone_size[i] <- 0
    }
    
    # If the TCR is in Unique_islets_TCR, call it 'matching'
    else if(TCR %in% Unique_islets_TCR){
      matching[i] <- "matching"
      clone_size[i] <- pbmc_sizes[pbmc_sizes$TCR == TCR, 'n']
    }
    
    # If the TCR is not in Unique_islets_TCR, call it non-matching'
    else {
      matching[i] <- "not_matching"
      clone_size[i] <- pbmc_sizes[pbmc_sizes$TCR == TCR, 'n']
    }
    
  }
  
  # Add the results back to the pbmc dataframe
  pbmc_df$Matching <- matching
  pbmc_df$Clone.size <- clone_size
  
  
  
  # Run matching analysis on islets
  # Initialize empty vectors
  matching <- character(length = nrow(islets_df))
  clone_size <- rep(0, nrow(islets_df))
  
  for(i in 1:length(all_islets_TCRs)){
    # For each TCR
    TCR <- all_islets_TCRs[i]
    
    # If it has no sequenced TCR, call it 'notcr' and give it a clone size of 0
    if(TCR == 'notcr'){
      matching[i] <- 'notcr'
      clone_size[i] <- 0
    }
    
    # If it's missing an alpha chain, call it 'beta_chain_only', give it a
    # clone size of 0
    else if(substr(TCR, 1, 1) == "|"){                   # Check first character of TCR
      matching[i] <- "beta_chain_only"
      clone_size[i] <- 0
    }
    
    # If it's missing a beta chain, call it 'alpha_chain_only', give it a
    # clone size of 0
    else if(substr(TCR, nchar(TCR), nchar(TCR)) == "|"){ # Check last character of TCR
      matching[i] <- "alpha_chain_only"
      clone_size[i] <- 0
    }
    
    # If the TCR is in Unique_blood_TCR, call it 'matching'
    else if(TCR %in% Unique_pbmc_TCR){
      matching[i] <- "matching"
      clone_size[i] <- islets_sizes[islets_sizes$TCR == TCR, 'n']
    }
    
    # If the TCR is not in Unique_pbmc_TCR, call it non-matching'
    else {
      matching[i] <- "not_matching"
      clone_size[i] <- islets_sizes[islets_sizes$TCR == TCR, 'n']
    }
    
  }
  
  # Add the results back to the islets dataframe
  islets_df$Matching <- matching
  islets_df$Clone.size <- clone_size
  
  
  
  # Now that the matching is done, add the dataframes back to the Seurat objects
  pbmc_rds@meta.data <- pbmc_df
  islets_rds@meta.data <- islets_df
  
  return(list(pbmc_rds, islets_rds))
  
}


#########################################################
#                                                       #
#     Run matching on filtered T1D PBMC and t1D islets rds objects        #
#                                                       #
#########################################################

# Run matching on filtered mouse rds objects
pbmc_samples <- c("A007", "A008", "A009", "A010", "A011", "A012")
islets_samples <- c("A001", "A002", "A003", "A004", "A005", "A006")

# Loop over each paired sample, run matching
for(i in 1:length(pbmc_samples)){
  
  # Load Seurat objects
  pbmc_rds <- readRDS(paste0(pbmc_samples[i], '_filtered.rds'))
  islets_rds <- readRDS(paste0(islets_samples[i], '_filtered.rds'))
  
  # Run matching
  output <- AB_Matching(pbmc_rds, islets_rds)
  
  pbmc_rds <- output[[1]]
  islets_rds <- output[[2]]
  
  # Save updated rds objects
  saveRDS(pbmc_rds, paste0('RDSobjects/', pbmc_samples[i], '.rds'))
  saveRDS(islets_rds, paste0('RDSobjects/', islets_samples[i], '.rds'))
  
}

