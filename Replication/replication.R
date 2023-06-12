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
# Case Study phmmPen_FA
######################################################################################################################
print("Case Study: phmmPen_FA results")
# Estimated $r$ for the Growth Ratio procedure

load("Paper_Results/PDAC_Selection_Results.RData")
# This will load the 'res' list object with the relevant output

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

# Bar graph summaries of fixed effects
for(i in 1:length(res)){
  fixef_all = res[[i]]$coef_vals
  fixef_non0 = fixef_all[which((fixef_all != 0))]
  df = data.frame(fixef = fixef_non0, TSP = str_sub(names(fixef_non0), start = 2))
  # str_sub(names(fixef_non0), start = 2) # names(fixef_non0)
  p = ggplot(data = df) + geom_col(mapping = aes(y = fixef, x = TSP)) +
    theme(axis.text.x = element_text(angle = 270)) + # , vjust = 0.5, hjust=1
    ylab("Log Hazard Ratio") 
  # ggsave(file = sprintf("Figures/PDAC_Fixef_Coef_%s.pdf",names(res)[i]),
  #        plot = p, units = "in", width = 6, height = 4)
  p = p + ggtitle(sprintf("%s fixed effects", names(res)[i]))
  print(p)
}

# Numerical summaries of selected fixed effects
for(i in 1:length(res)){
  fixef_all = res[[i]]$coef_vals
  fixef_non0 = fixef_all[which((fixef_all != 0))]
  df = data.frame(fit_type = names(res)[i], 
                  fixef = fixef_non0, TSP = str_sub(names(fixef_non0), start = 2))
  rownames(df) = NULL
  print(df)
}



# Overlap of fixed effects when comparing alpha values \{0.8,0.9\} for r = GR estimate
idx = seq(from = 3, to = 6, by = 2)
comp_idx = idx
fixef_overlap = NULL
fixef_any = NULL
df = NULL
for(j in 1:length(comp_idx)){
  fixef_all = res[[comp_idx[j]]]$coef_vals
  fixef_non0 = str_sub(names(fixef_all[which((fixef_all != 0))]),start=2)
  if(is.null(fixef_overlap)){
    fixef_overlap = fixef_non0
    fixef_any = fixef_non0
  }else{
    fixef_overlap = intersect(fixef_overlap, fixef_non0)
    fixef_any = union(fixef_any, fixef_non0)
  }
  
  df_tmp = data.frame(TSP = str_sub(names(fixef_all),start=2),
                      coef = fixef_all, Scenario = names(res)[comp_idx[j]])
  if(is.null(df)){
    df = df_tmp
  }else{
    
    df = rbind(df,df_tmp)
  }
}
# print(fixef_overlap)
# print(length(fixef_overlap))

df = df[which(df$TSP %in% fixef_any),]
p = ggplot(data = df) + geom_col(mapping = aes(y = coef, x = TSP, fill = Scenario),
                                 position = "dodge") +
  theme(axis.text.x = element_text(angle = 270)) + # , vjust = 0.5, hjust=1
  ylab("Log Hazard Ratio") +
  ggtitle("Overlapping Fixef: alpha = (0.8,0.9), r = GR")
print(p)


# Overlap of fixed effects when comparing GR vs r=3 for a particular set of conditions
idx = seq(from = 1, to = length(res), by = 2)
for(i in idx){
  # print(sprintf("Overlapping Fixef: %s",names(res)[i]))
  comp_idx = c(i,i+1)
  fixef_overlap = NULL
  fixef_any = NULL
  df = NULL
  for(j in 1:length(comp_idx)){
    fixef_all = res[[comp_idx[j]]]$coef_vals
    fixef_non0 = str_sub(names(fixef_all[which((fixef_all != 0))]),start=2)
    if(is.null(fixef_overlap)){
      fixef_overlap = fixef_non0
      fixef_any = fixef_non0
    }else{
      fixef_overlap = intersect(fixef_overlap, fixef_non0)
      fixef_any = union(fixef_any, fixef_non0)
    }
    
    df_tmp = data.frame(TSP = str_sub(names(fixef_all),start=2),
                        coef = fixef_all, Scenario = names(res)[comp_idx[j]])
    if(is.null(df)){
      df = df_tmp
    }else{
      
      df = rbind(df,df_tmp)
    }
  }
  # print(fixef_overlap)
  # print(length(fixef_overlap))
  
  df = df[which(df$TSP %in% fixef_any),]
  p = ggplot(data = df) + geom_col(mapping = aes(y = coef, x = TSP, fill = Scenario),
                                   position = "dodge") +
    theme(axis.text.x = element_text(angle = 270)) + # , vjust = 0.5, hjust=1
    ylab("Log Hazard Ratio") +
    ggtitle(sprintf("%s fixed effects", names(res)[i]))
  print(p)
}

# Overlap of fixed effects when comparing alpha values \{0.7,0.8,0.9\} for a particular value of r (GR estimate)

idx = seq(from = 1, to = 6, by = 2)
# print("Overlapping Fixef: alpha = (0.7,0.8,0.9), r = GR")
comp_idx = idx
fixef_overlap = NULL
fixef_any = NULL
df = NULL
for(j in 1:length(comp_idx)){
  fixef_all = res[[comp_idx[j]]]$coef_vals
  fixef_non0 = str_sub(names(fixef_all[which((fixef_all != 0))]),start=2)
  if(is.null(fixef_overlap)){
    fixef_overlap = fixef_non0
    fixef_any = fixef_non0
  }else{
    fixef_overlap = intersect(fixef_overlap, fixef_non0)
    fixef_any = union(fixef_any, fixef_non0)
  }
  
  df_tmp = data.frame(TSP = str_sub(names(fixef_all),start=2),
                      coef = fixef_all, Scenario = names(res)[comp_idx[j]])
  if(is.null(df)){
    df = df_tmp
  }else{
    
    df = rbind(df,df_tmp)
  }
}
# print(fixef_overlap)
# print(length(fixef_overlap))

df = df[which(df$TSP %in% fixef_any),]
p = ggplot(data = df) + geom_col(mapping = aes(y = coef, x = TSP, fill = Scenario),
                                 position = "dodge") +
  theme(axis.text.x = element_text(angle = 270)) + # , vjust = 0.5, hjust=1
  ylab("Log Hazard Ratio") +
  ggtitle("Overlapping Fixef: alpha = (0.7,0.8,0.9), r = GR")
print(p)



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
                  TSP = str_sub(names(var_non0), start = 2)) #  names(var_non0)
  rownames(df) = NULL
  print(df)
}


######################################################################################################################
# 
######################################################################################################################