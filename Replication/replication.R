# Code to replicate all tables in paper

# Load libraries
library(stringr)
library(xtable)
library(ggplot2)

# set working directory
setwd(path)

# True number of predictors with non-zero fixed effects
## Number of non-zero random effects: p_true + 1, which also includes a random intercept
p_true = 5 

######################################################################################################################
# Variable selection for piecewise exponential mixed effects simulations, p=100
######################################################################################################################

print("Variable selection for piecewise exponential mixed effects simulations, p=100")

load("Paper_Results/phmmPen_FA_p100.RData")
# This will load the 'res' list object with the relevant output

# True/False Positives and Timing
print(xtable(res$out_mat[,-c(1:(p_true))]))

# r estimate results
idx = seq(from = 2, to = length(res$rest_avg_pseudo), by = 2)
r_out = data.frame(Avg_r = res$rest_avg_pseudo, res$rest_pseudo)
print(xtable(r_out[idx,], digits = c(0,2,0,0,0)))

######################################################################################################################
# Variable selection for piecewise exponential mixed effects simulations, p=500
######################################################################################################################

print("Variable selection for piecewise exponential mixed effects simulations, p=500")

load("Paper_Results/phmmPen_FA_p500.RData")
# This will load the 'res' list object with the relevant output

# True/False Positives and Timing
print(xtable(res$out_mat[,-c(1:(p_true))]))

# r estimate results
idx = 1:8
r_out = data.frame(Avg_r = res$rest_avg_pseudo, res$rest_pseudo)
print(xtable(r_out[idx,], digits = c(0,2,0,0,0)))

######################################################################################################################
# Variable selection using ncvreg::cv.ncvsurv() (fixed-effects only)
######################################################################################################################

print("Variable selection using ncvreg fixed effects only selection, p=100")

load("Paper_Results/ncvsurv.RData")
# This will load the 'res' list object with the relevant output

# True/False Positives and Timing
print(xtable(res$out_mat[,-c(1:(p_true))]))


######################################################################################################################
# Variable selection applied to Weibull-simulated mixed effects survival data, p=100
######################################################################################################################

print("Weibull-simulated mixed effects survival data variable selection results, p=100")

load("Paper_Results/Weibull.RData")
# This will load the 'res' list object with the relevant output

# True/False Positives and Timing
print(xtable(res$out_mat[,-c(1:(p_true))]))

# r estimate results
r_out = data.frame(Avg_r = res$rest_avg_pseudo, res$rest_pseudo)
print(xtable(r_out, digits = c(0,2,0,0,0)))


######################################################################################################################
# Variable selection results when the number of fixed effects does not equal the number of random effects, p=100
######################################################################################################################

print("Variable selection results when the number of fixed effects does not equal the number of random effects, p=100")

p_true = 10

load("Paper_Results/alt_num_ranef.RData")
# This will load the 'res' list object with the relevant output

# True/False Positives and Timing
print(xtable(res$out_mat[,-c(1:(p_true))]))

# r estimate results
r_out = data.frame(Avg_r = res$rest_avg_pseudo, res$rest_pseudo)
print(xtable(r_out, digits = c(0,2,0,0,0)))



######################################################################################################################
# Variable selection results when purposefully underestimating the number of latent factors r, p=100
######################################################################################################################

print("Variable selection results when purposefully underestimating the number of latent factors r, p=100")

p_true = 5

load("Paper_Results/alt_rval.RData")
# This will load the 'res' list object with the relevant output

# True/False Positives and Timing
print(xtable(res$out_mat[,-c(1:(p_true))]))

######################################################################################################################
# Variable selection results using different numbers of time intervals J=5,6, p=100
######################################################################################################################

print("Variable selection results using different numbers of time intervals J=5,6, p=100")

load("Paper_Results/alt_J_5and6.RData")
# This will load the 'res' list object with the relevant output

# True/False Positives and Timing
print(xtable(res$out_mat[,-c(1:(p_true))]))

# r estimate results
r_out = data.frame(Avg_r = res$rest_avg_pseudo, res$rest_pseudo)
print(xtable(r_out, digits = c(0,2,0,0,0)))

######################################################################################################################
# Variable selection results using different numbers of time intervals J=9,10, p=100
######################################################################################################################

print("Variable selection results using different numbers of time intervals J=9,10, p=100")

load("Paper_Results/alt_J_9and10.RData")
# This will load the 'res' list object with the relevant output

# True/False Positives and Timing
print(xtable(res$out_mat[,-c(1:(p_true))]))

# r estimate results
r_out = data.frame(Avg_r = res$rest_avg_pseudo, res$rest_pseudo)
print(xtable(r_out, digits = c(0,2,0,0,0)))



######################################################################################################################
# Case Study phmmPen_FA
######################################################################################################################
print("Case Study: phmmPen_FA results")
# Estimated $r$ for the Growth Ratio procedure

load("Paper_Results/PDAC_Selection_Results_revision.RData")
# This will load the 'res' list object with the relevant output

# C-index values for each combination of elastic net parameter and latent factor r values
cidx_mat = matrix(NA, nrow = length(res), ncol = 1)
colnames(cidx_mat) = c("phmmPen_FA")
rownames(cidx_mat) = names(res)
for(i in 1:length(res)){
  cidx_mat[i,1] = round(res[[i]]$c_index[1],4)
}
print(cidx_mat)

# Estimated r from the Growth Ratio procedure
idx = which(str_detect(names(res),"GR_est"))
for(i in idx){
  print(sprintf("%s: r = %i", names(res)[i], res[[i]]$r_est))
}

# Time to complete the algorithm in hours
for(i in 1:length(res)){
  print(sprintf("%s: hours = %.1f", names(res)[i], res[[i]]$time_mat[,3] / 3600))
}

# Number of non-zero fixed effects estimates in the best model
for(i in 1:length(res)){
  print(sprintf("%s: number non-zero fixed effects = %i", names(res)[i], sum(res[[i]]$coef_vals != 0)))
}

# Select output reported in main paper
## See "replication_case_study_sensitivity.R" for additional output
i = which(names(res) == "Alpha_0.9_r_3")

# Bar graph summaries of fixed effects
  fixef_all = res[[i]]$coef_vals
  fixef_non0 = fixef_all[which((fixef_all != 0))]
  df = data.frame(fixef = fixef_non0, TSP = names(fixef_non0))
  p = ggplot(data = df) + geom_col(mapping = aes(y = fixef, x = TSP)) +
    theme(axis.text.x = element_text(angle = 270)) + # , vjust = 0.5, hjust=1
    ylab("Log Hazard Ratio") 
  # ggsave(file = sprintf("Figures/PDAC_Fixef_Coef_%s.pdf",names(res)[i]),
  #        plot = p, units = "in", width = 6, height = 4)
  p = p + ggtitle(sprintf("%s fixed effects", names(res)[i]))
  print(p)


# Numerical summaries of selected fixed effects
  fixef_all = res[[i]]$coef_vals
  fixef_non0 = fixef_all[which((fixef_all != 0))]
  df = data.frame(fit_type = names(res)[i], 
                  fixef = fixef_non0, TSP = names(fixef_non0))
  rownames(df) = NULL
  print(df)



# Number non-zero random effects in best model:
  
for(i in 1:length(res)){
  print(sprintf("%s: number non-zero random effects = %i", names(res)[i],
                sum(as.numeric(res[[i]]$vars_vals[-1]) != 0)))
}

# Summary of random effect variance estimates:
  
for(i in 1:length(res)){
  var_all = res[[i]]$vars_vals
  var_non0 = var_all[which((var_all != 0))]
  df = data.frame(fit = names(res)[i], var_val = var_non0,
                  TSP = c("(Intercept)",str_sub(names(var_non0)[-1], start = 8))) #  names(var_non0)
  rownames(df) = NULL
  print(df)
}


######################################################################################################################
# 
######################################################################################################################