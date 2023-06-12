# Merge individual-study pancreatic cancer datasets into single dataset 
# with gene expression information

# load libraries
library(stringr)
library(coxme)
library(survival)

# 1-44
array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Specify path to saved step01 output
path_data = "~/Project3/PDAC_phmmPen_FA_TSP"
# Specify path to output data
path_output = str_c(path_data,"/TSP_loglik_v2")
if(!dir.exists(path_output)){dir.create(path_output, recursive = TRUE)}

# load ex_TSP data.frame object
load(file = sprintf("%s/step01_TSP.RData",path_data))
dim(ex_TSP)

# Non-TSP variable names
sampInfo_keep = c("sampID","study","time","event")
# TSP variable names
TSP_names = colnames(ex_TSP)[which(!(colnames(ex_TSP) %in% sampInfo_keep))]
length(TSP_names)

#################################################################################################
# Compute log-likelihood for all TSPs
## Will use these log-likelhoods in Step03
#################################################################################################

batch = array_val

if(batch <= 43){
  idx = 1:1000 + (batch-1)*1000
}else if(batch == 44){
  idx = 43001:length(TSP_names)
}

if(batch <= 9){
  batch_label = str_c("00",batch)
}else if(batch <= 99){
  batch_label = str_c("0",batch)
}else{
  batch_label = as.character(batch)
}

# data.frame to use in fitting the following logistic mixed effects model:
## y ~ tsp + (tsp | study)
ex_fit = ex_TSP

# Calculate log-likelihood for all TSPs
loglik_vals = numeric(length(idx))
names(loglik_vals) = TSP_names[idx]
for(j in 1:length(idx)){
  t = idx[j]
  print(sprintf("t %i",t))
  tsp = TSP_names[t]
  y = Surv(ex_TSP$time, ex_TSP$event)
  x = ex_fit[,tsp]
  study = factor(ex_fit$study)
  # Fit model with tsp covariate, random intercept, and random effect for tsp covariate
  fit_coxme = try(coxme(y ~ x + (1 + x | study)))
  if(!inherits(fit_coxme,"try-error")){
    loglik_vals[j] = logLik(fit_coxme)
  }else{
    loglik_vals[j] = NA
  }
  
}
save(loglik_vals, file = sprintf("%s/step02_TSP_loglik_Batch%s.RData", path_output, batch_label))



################################################################################################

print(gc(full = TRUE))

q(save="no")

################################################################################################