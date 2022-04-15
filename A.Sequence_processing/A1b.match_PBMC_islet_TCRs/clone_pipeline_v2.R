#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

filtered_contig_annotations_islet_file = args[1]
filtered_contig_annotations_pbmc_file = args[2]
output_folder = args[3]
pool_num = args[4]

filtered_contig_annotations_islet = read.csv(filtered_contig_annotations_islet_file)
filtered_contig_annotations_pbmc = read.csv(filtered_contig_annotations_pbmc_file)

# just include info I want in filtered_contig_annotations
filtered_contig_annotations_islet = filtered_contig_annotations_islet[,c("barcode","contig_id","chain","cdr3","umis","raw_clonotype_id")]
filtered_contig_annotations_pbmc = filtered_contig_annotations_pbmc[,c("barcode","contig_id","chain","cdr3","umis","raw_clonotype_id")]

# next step is to have a row for each barcode with all chains collapsed together. currently alpha and beta chains are in diff rows
# if there are more than one alpha and/or beta chains, pick the chain with the highest UMI
# if there are more than one chain with the max UMI keep both
filtered_contig_annotations_islet_split = split(filtered_contig_annotations_islet,filtered_contig_annotations_islet$barcode)
filtered_contig_annotations_pbmc_split = split(filtered_contig_annotations_pbmc,filtered_contig_annotations_pbmc$barcode)

get_max_umi_barcodes = function(my_list) {
  barcode_chain_info = lapply(my_list, function(x) {
    # first get TRA chain that has highest umi
    barcode_temp = unique(x$barcode)
    max_umi_tra_chain = ""
    max_umi_trb_chain = ""
    if("TRA"%in%x$chain) {      
      tra_df = x[x$chain=="TRA",]
      max_tra_umi_num = max(tra_df$umis)
      tra_max_umis = tra_df[tra_df$umis==max_tra_umi_num,]
      max_umi_tra_chain = paste(sort(tra_max_umis$cdr3),collapse=",")
    }
    if ("TRB"%in%x$chain) {
      trb_df = x[x$chain=="TRB",]
      max_trb_umi_num =	max(trb_df$umis)
      trb_max_umis = trb_df[trb_df$umis==max_trb_umi_num,]
      max_umi_trb_chain = paste(sort(trb_max_umis$cdr3),collapse=",")
    }
    if(max_umi_tra_chain == "" & max_umi_trb_chain == "") {
      trb_tra = "notcr"
    } else {
      trb_tra = paste(max_umi_tra_chain,max_umi_trb_chain,sep="|")
    }
    if(max_umi_tra_chain == "") {
      max_umi_tra_chain = "notcr"
    }
    if(max_umi_trb_chain == "") {
      max_umi_trb_chain = "notcr"
    }
    barcode_with_umi_resolved_chains = c(barcode_temp,max_umi_tra_chain,max_umi_trb_chain,trb_tra)
    return(barcode_with_umi_resolved_chains)
  })
barcode_chain_info = do.call("rbind",barcode_chain_info)
return(barcode_chain_info)
}

islet_barcode_chain_info = get_max_umi_barcodes(filtered_contig_annotations_islet_split)
pbmc_barcode_chain_info = get_max_umi_barcodes(filtered_contig_annotations_pbmc_split)

# get frequency of chains
islet_beta_freq = table(islet_barcode_chain_info[,3])
islet_alpha_beta_freq = table(islet_barcode_chain_info[,4])
pbmc_beta_freq = table(pbmc_barcode_chain_info[,3])
pbmc_alpha_beta_freq = table(pbmc_barcode_chain_info[,4])
# if notcr is found set frequency to 0. otherwise frequency will be the number of cells with no TCR
islet_beta_freq[names(islet_beta_freq) == "notcr"] = 0
islet_alpha_beta_freq[names(islet_alpha_beta_freq) == "notcr"] = 0
pbmc_beta_freq[names(pbmc_beta_freq) == "notcr"] = 0
pbmc_alpha_beta_freq[names(pbmc_alpha_beta_freq) == "notcr"] = 0

## put frequencies in tables
islet_barcode_chain_info = cbind(islet_barcode_chain_info,islet_beta_freq[match(islet_barcode_chain_info[,3],names(islet_beta_freq))])
islet_barcode_chain_info = cbind(islet_barcode_chain_info,islet_alpha_beta_freq[match(islet_barcode_chain_info[,4],names(islet_alpha_beta_freq))])
pbmc_barcode_chain_info = cbind(pbmc_barcode_chain_info,pbmc_beta_freq[match(pbmc_barcode_chain_info[,3],names(pbmc_beta_freq))])
pbmc_barcode_chain_info	= cbind(pbmc_barcode_chain_info,pbmc_alpha_beta_freq[match(pbmc_barcode_chain_info[,4],names(pbmc_alpha_beta_freq))])


# get matches based on alpha and beta chain
get_alpha_plus_beta_matches = function(mydf,other_df) {
  matches_df = apply(mydf, 1, function(myrow) {
    alpha_beta_chain = myrow[4]
    # only compare if we have either an alpha and beta chain
    if(alpha_beta_chain != "notcr") {
      other_df_matches = other_df[other_df[,4] == alpha_beta_chain,,drop=FALSE]
      other_df_matches_names = rownames(other_df_matches)
      other_df_matches_names_pasted = paste(other_df_matches_names,collapse="_")
      if(length(other_df_matches_names) >0) {
        return(c(myrow[1],"matching",other_df_matches_names_pasted))
      } else {
        return(c(myrow[1],"not_matching",other_df_matches_names_pasted))
      }
    } else {
      return(c(myrow[1],"notcr",""))
    }
  })
}

pbmc_alpha_beta_matches = t(get_alpha_plus_beta_matches(pbmc_barcode_chain_info,islet_barcode_chain_info))
islet_alpha_beta_matches = t(get_alpha_plus_beta_matches(islet_barcode_chain_info,pbmc_barcode_chain_info))

## Next get matches by beta chain only

get_beta_matches = function(mydf,other_df) {
  matches_df = apply(mydf, 1, function(myrow) {
    beta_chain = myrow[3]
    # only compare if we have a beta chain
    if(beta_chain != "notcr") {
      other_df_matches = other_df[other_df[,3] == beta_chain,,drop=FALSE]
      other_df_matches_names = rownames(other_df_matches)
      other_df_matches_names_pasted = paste(other_df_matches_names,collapse="_")
      if(length(other_df_matches_names) >0) {
        return(c(myrow[1],"beta_matching",other_df_matches_names_pasted))
      } else {
        return(c(myrow[1],"not_matching",other_df_matches_names_pasted))
      }
    } else {
      return(c(myrow[1],"notcr",""))
    }
  })
}

pbmc_beta_matches = t(get_beta_matches(pbmc_barcode_chain_info,islet_barcode_chain_info))
islet_beta_matches = t(get_beta_matches(islet_barcode_chain_info,pbmc_barcode_chain_info))

# create output dataframes
pbmc_barcode_chain_info_match_info_beta_alpha = cbind(pbmc_barcode_chain_info[,1],pbmc_alpha_beta_matches[,2],pbmc_barcode_chain_info[,6],pbmc_barcode_chain_info[,4])
pbmc_barcode_chain_info_match_info_beta = cbind(pbmc_barcode_chain_info[,1],pbmc_beta_matches[,2],pbmc_barcode_chain_info[,5],pbmc_barcode_chain_info[,3])
islet_barcode_chain_info_match_info_beta_alpha = cbind(islet_barcode_chain_info[,1],islet_alpha_beta_matches[,2],islet_barcode_chain_info[,6],islet_barcode_chain_info[,4])
islet_barcode_chain_info_match_info_beta = cbind(islet_barcode_chain_info[,1],islet_beta_matches[,2],islet_barcode_chain_info[,5],islet_barcode_chain_info[,3])

colnames(pbmc_barcode_chain_info_match_info_beta_alpha) = c("Barcode","Matching_pre_filter","Frequency_pre_filter","TCR")
colnames(pbmc_barcode_chain_info_match_info_beta) = c("Barcode","Matching_pre_filter","Frequency_pre_filter","TCR")
colnames(islet_barcode_chain_info_match_info_beta_alpha) = c("Barcode","Matching_pre_filter","Frequency_pre_filter","TCR")
colnames(islet_barcode_chain_info_match_info_beta) = c("Barcode","Matching_pre_filter","Frequency_pre_filter","TCR")

dir.create(paste(output_folder,"/alpha_beta_TCR",sep=""))
dir.create(paste(output_folder,"/beta_TCR",sep=""))


write.csv(pbmc_barcode_chain_info_match_info_beta_alpha,file=paste(output_folder,"/alpha_beta_TCR/PBMC_alpha_beta_TCR_matches_pool_",pool_num,".csv",sep=""),row.names=F)
write.csv(islet_barcode_chain_info_match_info_beta_alpha,file=paste(output_folder,"/alpha_beta_TCR/islet_alpha_beta_TCR_matches_pool_",pool_num,".csv",sep=""),row.names=F)
write.csv(pbmc_barcode_chain_info_match_info_beta,file=paste(output_folder,"/beta_TCR/PBMC_beta_TCR_matches_pool_",pool_num,".csv",sep=""),row.names=F)
write.csv(islet_barcode_chain_info_match_info_beta,file=paste(output_folder,"/beta_TCR/islet_beta_TCR_matches_pool_",pool_num,".csv",sep=""),row.names=F)



