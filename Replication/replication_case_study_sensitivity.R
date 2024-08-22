# Code to replicate all tables in paper

# Load libraries
library(stringr)
library(xtable)
library(ggplot2)

# set working directory
setwd(path)


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

# Bar graph summaries of fixed effects
for(i in 1:length(res)){
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
}

# Numerical summaries of selected fixed effects
for(i in 1:length(res)){
  fixef_all = res[[i]]$coef_vals
  fixef_non0 = fixef_all[which((fixef_all != 0))]
  df = data.frame(fit_type = names(res)[i], 
                  fixef = fixef_non0, TSP = names(fixef_non0))
  rownames(df) = NULL
  print(df)
}

# Overlap of fixed effects when comparing alpha values \{0.7,0.8,0.9\} for manually-set r=3
idx = seq(from = 2, to = 7, by = 2)
# print("Overlapping Fixef: alpha = (0.7,0.8,0.9), r = GR")
comp_idx = idx
fixef_overlap = NULL
fixef_any = NULL
df = NULL
for(j in 1:length(comp_idx)){
  fixef_all = res[[comp_idx[j]]]$coef_vals
  fixef_non0 = str_sub(names(fixef_all[which((fixef_all != 0))]),start=1)
  if(is.null(fixef_overlap)){
    fixef_overlap = fixef_non0
    fixef_any = fixef_non0
  }else{
    fixef_overlap = intersect(fixef_overlap, fixef_non0)
    fixef_any = union(fixef_any, fixef_non0)
  }
  
  df_tmp = data.frame(TSP = str_sub(names(fixef_all),start=1),
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
  ggtitle("Overlapping Fixef: alpha = (0.7,0.8,0.9), r = 3")
print(p)




# Overlap of fixed effects when comparing alpha values \{0.7,0.8,0.9\} for r = GR estimate
idx = seq(from = 1, to = 6, by = 2)
# print("Overlapping Fixef: alpha = (0.7,0.8,0.9), r = GR")
comp_idx = idx
fixef_overlap = NULL
fixef_any = NULL
df = NULL
for(j in 1:length(comp_idx)){
  fixef_all = res[[comp_idx[j]]]$coef_vals
  fixef_non0 = str_sub(names(fixef_all[which((fixef_all != 0))]),start=1)
  if(is.null(fixef_overlap)){
    fixef_overlap = fixef_non0
    fixef_any = fixef_non0
  }else{
    fixef_overlap = intersect(fixef_overlap, fixef_non0)
    fixef_any = union(fixef_any, fixef_non0)
  }
  
  df_tmp = data.frame(TSP = str_sub(names(fixef_all),start=1),
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



# Overlap of fixed effects when comparing GR vs r=3 for a particular value of elastic net parameter (0.7 to 1.0)
idx = seq(from = 1, to = length(res), by = 2)
for(i in idx){
  # print(sprintf("Overlapping Fixef: %s",names(res)[i]))
  comp_idx = c(i,i+1)
  fixef_overlap = NULL
  fixef_any = NULL
  df = NULL
  for(j in 1:length(comp_idx)){
    fixef_all = res[[comp_idx[j]]]$coef_vals
    fixef_non0 = str_sub(names(fixef_all[which((fixef_all != 0))]),start=1)
    if(is.null(fixef_overlap)){
      fixef_overlap = fixef_non0
      fixef_any = fixef_non0
    }else{
      fixef_overlap = intersect(fixef_overlap, fixef_non0)
      fixef_any = union(fixef_any, fixef_non0)
    }
    
    df_tmp = data.frame(TSP = str_sub(names(fixef_all),start=1),
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