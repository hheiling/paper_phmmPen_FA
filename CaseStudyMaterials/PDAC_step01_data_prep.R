# Merge individual-study pancreatic cancer datasets into single dataset 
# with gene expression information

# load libraries
library(stringr)

# Specify path to working directory
path = "~/GitHub/cloud/proj3/PDAC"
# Specify path to data of interest
path_data = sprintf("%s/data_asis_refseq_star_rsem", path)
# Save step01 results
path_save_results = sprintf("%s/PDAC_TSP",path)

# Study labels
data_names = c("Aguirre",
               "CPTAC",
               "Dijk",
               "Moffitt_GEO_array",
               "PACA_AU_array",
               "Puleo_array",
               "TCGA_PAAD")
# gene expression data for each study
datasets = str_c(data_names,".rds")
# survival data for each study
data_surv = str_c(data_names,".survival_data.rds")

# Relevant genes to focus on in analyses.
genelist = read.table(sprintf("%s/v3_20140714_all_filter-14-200-genelist-TumorType_spar250.txt", path),
                      header = TRUE)

# record the 500 tumor-related genes of interest
symbols_500 = genelist$symbol

#################################################################################################
# Clean up expression datasets 
#################################################################################################

# vector of common symbols across all datasets, and contained within the genelist specified above
## Initialize with the 500 gene list
symbols_common = symbols_500

# Save original datasets
data_lst = list()
# Save cleaned datasets
data_clean = list()

for(d in 1:length(datasets)){ 
  cat(sprintf("\n Start dataset %i \n",d))
  
  # Extract relevant dataset and save
  data = readRDS(sprintf("%s/%s",path_data,datasets[d]))
  data_lst[[data_names[d]]] = data
  
  # Extract relevant output from data list object
  ex = data$ex
  featInfo = data$featInfo
  SYMBOL_all = featInfo$SYMBOL
  
  # Extract sample IDs from dataset
  sampID = colnames(ex)
  
  cat(sprintf("Subjects total: %i \n",length(sampID)))
  
  # unique gene symbols present in dataset
  cat(sprintf("Total gene symbols in dataset: %i \n", length(SYMBOL_all)))
  symbols = unique(intersect(SYMBOL_all, symbols_500))
  cat(sprintf("Total symbols in dataset to use (unique, within 500 tumor gene list): %i \n", length(symbols)))
  # vector to save names of rows in cleaned ex matrix
  ex_clean_rows = character(length(symbols))
  
  for(s in 1:length(symbols)){
    # select gene symbol of interest
    sym = symbols[s]
    # rownames of the ex matrix - gene symbol names
    ex_rows = SYMBOL_all
    # select row or rows corresponding to specified gene symbol
    ex_tmp = ex[which(ex_rows == sym),1:ncol(ex),drop=FALSE]
    
    if(nrow(ex_tmp) == 1){
      # If single single row corresponding to gene symbol, not additional processing necessary
      ex_use = ex_tmp
    }else if(nrow(ex_tmp) >= 2){
      # If 2+ rows corresponding to gene symbol, do the following:
      #   Average across values if microarray dataset
      #   Sum across values if RNA-seq dataset
      # Note: subjects on columns, genes on rows. Therefore, within a particular column (subject),
      # add/average values
      if(str_detect(datasets[d],"array")){
        ex_use = colMeans(ex_tmp)
      }else{
        ex_use = colSums(ex_tmp)
      }
    }
    
    if(s == 1){
      ex_clean = ex_use
    }else{
      ex_clean = rbind(ex_clean, ex_use)
    }
    
    ex_clean_rows[s] = sym
    
  }
  # assign appropriate symbol names to rows in cleaned ex matrix
  rownames(ex_clean) = ex_clean_rows
  
  
  # Extract relevant subject-level information
  sampInfo = data.frame(sampID = sampID, study = rep(data_names[d],times=length(sampID)))
  rownames(sampInfo) = sampID
  
  # combine relevant output into list object: ex_clean plus relevant sample-level information
  dataSet = list(ex = ex_clean, sampInfo = sampInfo)
  
  data_clean[[data_names[d]]] = dataSet 
  
  # Update list of common genes across processed datasets
  symbols_common = intersect(symbols_common, ex_clean_rows)
  
}

length(symbols_common) # 420

#################################################################################################
# Merge expression data with survival data 
#################################################################################################

data_merge = list()

for(d in 1:length(datasets)){
  # Extract gene expression information
  ex = data_clean[[data_names[d]]]$ex
  # Extract sampID and study information
  sampInfo = data_clean[[data_names[d]]]$sampInfo
  # Transpose gene expression information so that both gene expression info and 
  #   sample information is arraged with rows = subjects, columns = genes/covariates
  ex_trans = data.frame(sampInfo, t(ex))
  # Load survival information
  surv_info = readRDS(sprintf("%s/%s",path_data,data_surv[d]))
  # Merge sample information, gene expression information, and survival information together
  ex_merge = merge(surv_info, ex_trans, by = "sampID")
  
  data_merge[[data_names[d]]] = ex_merge
  
}

#################################################################################################
# Remove subjects with problematic data, create single dataset
# Check sufficient events within groups
## Issue: whitelist = FALSE (results not from a primary tumor)
## Issue: missing survival data or times == 0
# Combine the datasets, only using the symbols that are (a) contained in the above genelist
# and (b) present in all datasets (i.e. in the symbols_common vector)
#################################################################################################

# ex matrix with all observations, using a common set of genes, and containing 
#   the relevant sample-level information as well
ex_full = NULL

# non-gene covariates to keep
sampInfo_keep = c("sampID","study","time","event")
# total covariates to keep
cols_use = c(sampInfo_keep, symbols_common)

for(d in 1:length(datasets)){
  tmp0 = data_merge[[data_names[d]]]
  tmp1 = tmp0[which(tmp0$whitelist == TRUE),]
  if(any(is.na(tmp1$time)) | any(tmp1$time == 0)){
    tmp2 = tmp1[which(!is.na(tmp1$time) & (tmp1$time > 0)),]
  }else{
    tmp2 = tmp1
  }
  
  cat(sprintf("%s sample size: %i \n", data_names[d], nrow(tmp2)))
  cat(sprintf("%s events: %i \n", data_names[d], sum(as.numeric(tmp2$event))))
  
  df_use = tmp2[,cols_use]
  if(d == 1){
    ex_full = df_use
  }else{
    ex_full = rbind(ex_full, df_use) # Columns are all in the same order
  }
  
}

# Check important covariates are of correct format
ex_full$event = as.numeric(ex_full$event)
ex_full$study = factor(ex_full$study)
ex_full$time = as.numeric(ex_full$time)

dim(ex_full) # ncol: length(symbols_common) + 4; nrow: sum of all observations, 879

colnames(ex_full)[which(!(colnames(ex_full) %in% symbols_common))]
head(colnames(ex_full))

table(ex_full$study, ex_full$event) # 0 = censored, 1 = event

mean(ex_full$event)
sum(ex_full$event) # 539

save(ex_full, file = sprintf("%s/ex_full.RData",path_save_results))

#################################################################################################
# Filter genes based on expression levels
## Rank transform the expression within each sample
## Compute mean rank for each gene
## Drop bottom 20% of genes based on these mean ranks
## Check: plot log mean expression (x) vs mean rank (y) 
#################################################################################################

## Note for ex_full data.frame: rows = subjects, columns = genes and other covariates

# Select gene expression columns
ex_tmp = ex_full[,which(colnames(ex_full) %in% symbols_common)]
# Calculate mean expression for each gene across all subjects
ex_avg = colMeans(ex_tmp)
# Set-up rank values for each subject
ex_rank = matrix(0, nrow = nrow(ex_full), ncol = ncol(ex_tmp))
colnames(ex_rank) = colnames(ex_tmp)

# Compute rank information
for(i in 1:nrow(ex_tmp)){
  ex_rank[i,] = rank(ex_tmp[i,])
}

rank_avg = colMeans(ex_rank)
names(rank_avg) = colnames(ex_rank)
plot(x = log(ex_avg), y = rank_avg) # Approx linear relationship, removing bottom 20% should be fine


# Order average ranks (lowest average rank to highest average rank)
rank_avg_order = rank_avg[order(rank_avg)]
# Remove bottom 20% of genes based on mean rank
upper80 = names(rank_avg_order[-c(1:round(0.2*length(rank_avg_order)))])
length(upper80)

# Create updated dataset with these genes removed
ex_final = ex_full[,c(sampInfo_keep,upper80)]
dim(ex_final)

# Set ex_tmp to NULL
ex_tmp = NULL

#################################################################################################
# Create TSP covariates
#################################################################################################

# train_sub is the n x p expression matrix 
# mat is a two column matrix, where each row is a tsp, column 1 is gene A in the TSP and 
#   column 2 is gene B. This can be enumerated from the list of candidate genes
ind_fun = function(train_sub, mat){
  indmat = matrix(-1, nrow = nrow(train_sub), ncol = nrow(mat))
  for(i in 1:nrow(mat)){
    p1 = which(colnames(train_sub) == mat[i,1])
    p2 = which(colnames(train_sub) == mat[i,2])
    indmat[,i] = (train_sub[,p1] > train_sub[,p2])^2
  }
  rownames(indmat) = rownames(train_sub)
  colnames(indmat) = str_c(mat[,1],"_",mat[,2])
  return(indmat)
}

TSP = t(combn(upper80, 2))
dim(TSP)

ex_TSP0 = ex_final[,which(colnames(ex_final) %in% upper80)]
ex_TSP1 = ind_fun(train_sub = ex_TSP0, mat = TSP)

TSP_names = colnames(ex_TSP1)

ex_TSP = data.frame(ex_final[,c(sampInfo_keep)], ex_TSP1)
dim(ex_TSP)

#################################################################################################
# Remove TSPs that do not have enough variation in the data
## If the mean TSP < 0.1 or > 0.9, then remove
#################################################################################################

TSP_keep = numeric(length(TSP_names))
names(TSP_keep) = TSP_names
for(t in 1:length(TSP_names)){
  tsp = TSP_names[t]
  tsp_mean = mean(ex_TSP[,tsp])
  if((tsp_mean >= 0.10) & (tsp_mean <= 0.90)){
    TSP_keep[t] = 1
  }
}

sum(TSP_keep)

symbols_use_TSP = names(TSP_keep[which(TSP_keep == 1)])
length(symbols_use_TSP)

ex_TSP = ex_TSP[,c(sampInfo_keep,symbols_use_TSP)]
dim(ex_TSP)

#################################################################################################
# Save step01 results
#################################################################################################

# Save step01 result
save(ex_TSP, file = sprintf("%s/step01_TSP.RData",path_save_results))

#################################################################################################

#################################################################################################