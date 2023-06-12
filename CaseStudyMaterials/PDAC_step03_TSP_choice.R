# Merge individual-study pancreatic cancer datasets into single dataset 
# with gene expression information

# load libraries
library(stringr)
library(reshape2)

# Path to saved step01 and step02 results
path_data = "~/GitHub/cloud/proj3/PDAC/PDAC_TSP"

files = list.files(path = sprintf("%s/TSP_loglik_v2",path_data), pattern = ".RData", full.names = TRUE)

# load ex_TSP data.frame from step01 output
load(sprintf("%s/step01_TSP.RData", path_data))

#################################################################################################
# Load and merge TSP log-likelihood results into single vector
## If in step02 calculated all log-likelihoods at once instead of performing parallel
## calculations, this step can be skipped/modified
#################################################################################################

TSP_loglik = NULL
for(i in 1:length(files)){
  # load loglik_vals object
  load(file = files[i])
  TSP_loglik = c(TSP_loglik, loglik_vals)
}

length(TSP_loglik)
summary(TSP_loglik)

dim(ex_TSP)

# Remove NA values, if any
## Note: Found no NA values
# TSP_loglik = TSP_loglik[which(!is.na(TSP_loglik))]
# length(TSP_loglik)

#################################################################################################
# Filter out TSPs with correlation due to shared genes
## Order step02 log-likelihoods
## If a gene is in a top TSP, remove all other TSPs with less-optimal
##    log-likelihoods that share one of the genes in this top TSP
#################################################################################################

# Goals:
## Choose top TSPs that maximize the log-likelihood.
## Remove TSPs with high correlation: if gene is in a top TSP, remove all other TSPs with
##    less optimal log-likelihoods that share this gene
## Maximum unique genes: 336
## Goal: Pick TSPs with greatest log-likelihoods (i.e. maximize the log-likelihoods) 
loglik_order = TSP_loglik[order(TSP_loglik, decreasing = TRUE)]

TSP_keep_ll = loglik_order
TSP_gene_names = str_split(names(TSP_keep_ll),pattern="_",n=2,simplify=TRUE)
length(unique(c(TSP_gene_names))) # 336
for(i in 1:length(unique(c(TSP_gene_names)))){
  tsp = names(TSP_keep_ll)[i]
  gene1 = TSP_gene_names[i,1]
  gene2 = TSP_gene_names[i,2]
  tsp_shared_genes = which((TSP_gene_names[,1] == gene1) | (TSP_gene_names[,2] == gene2) | (TSP_gene_names[,1] == gene2) | (TSP_gene_names[,2] == gene1))
  tsp_remove = tsp_shared_genes[which(tsp_shared_genes > i)]
  if(length(tsp_remove) >= 1){
    TSP_keep_ll = TSP_keep_ll[-tsp_remove]
    TSP_gene_names = TSP_gene_names[-tsp_remove,]
  }
  if(length(TSP_keep_ll) == i){
    break
  }
}

length(TSP_keep_ll) # 168
dim(TSP_gene_names)

# Check if duplicate gene names included
## Should have 2 * length(TSP_keep_ll) unique gene names
length(unique(c(TSP_gene_names)))

TSP_gene_names_use = names(TSP_keep_ll)

#################################################################################################
# Check correlations of TSPs kept in model
## Possibly: Remove TSPs with too high of correlation
#################################################################################################

sampInfo_keep = c("sampID","study","time","event")
# Remove sample information variables, keep gene expression only
ex_tmp = ex_TSP[,TSP_gene_names_use]
# Calculate spearman correlations
## Examine absolute values
cor_mat = abs(cor(ex_tmp, method = "spearman")) 
cor_mat[lower.tri(cor_mat, diag = TRUE)] = NA
summary(c(cor_mat))
hist(c(cor_mat))

df = melt(cor_mat, na.rm = TRUE)
df[which(df$value > 0.9),]
df[which((df$value > 0.8) & (df$value < 0.9)),]
cor_limit = 0.8
length(unique(c(df$Var1[which(df$value > cor_limit)], df$Var2[which(df$value > cor_limit)])))

# No additional removals of TSP covariates needed based on excessively high correlation (cor > 0.8)
TSP_gene_names_final = TSP_gene_names_use

length(TSP_gene_names_use)
length(TSP_gene_names_final)

#################################################################################################
# Save final dataset: PDAC_surv
## Save as data.frame with: sample IDs, study, time, event, and relevant TSPs
#################################################################################################

PDAC_surv = ex_TSP[,c("sampID","study","time","event",TSP_gene_names_final)]
class(PDAC_surv)
# Check dimensions
dim(PDAC_surv)
# Save
save(PDAC_surv, file = sprintf("%s/PDAC_surv_TSP.RData", path_data))

#################################################################################################
# Descriptive statistics per study - sample size, number events
#################################################################################################

# Sample size
table(PDAC_surv$study)
# Number events
table(PDAC_surv$study, PDAC_surv$event)

#################################################################################################
# 
#################################################################################################