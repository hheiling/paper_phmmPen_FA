# Compile data from the simulation results

# Define directory where all simulation results are kept
path_sims = "~/"
# Define directory where to store simulation results
path_output = "~/Paper_Results"

# extract_r: Function to extract relevant information from simulation results
extract_r = function(path, sims, sim_type,
                     p_tot = 100, p_true = 5, q = 6,
                     r_true_vec = 3, trace = 0,
                     beta_true = rep(1,length(sims))){
  
  # Results of r estimation procedure
  ## Percent of times too small, too large, or correct
  rest_pseudo = matrix(0, nrow = length(sims), ncol = 3)
  rownames(rest_pseudo) = sim_type
  colnames(rest_pseudo) = c("Too_small_Pct","Correct_Pct","Too_Large_Pct")
  ## Average estimate of r
  rest_avg_pseudo = numeric(length(sims))
  names(rest_avg_pseudo) = sim_type
  
  rest_lst = list()
  
  # Save summary output in out_mat
  ## labs = labels of columns of out_mat
  labs = c(str_c("Beta",1:(p_true)),"TP % Fixef","FP % Fixef",
           "TP % Ranef","FP % Ranef",
           "Median Time (hrs)","Abs. Dev. (Mean)")
  out_mat = matrix(0, nrow = length(sims), ncol = length(labs))
  colnames(out_mat) = labs
  rownames(out_mat) = sim_type
  
  # Save additional variance and fixed effect information
  vars_lst = list()
  beta_lst = list()
  
  r_input_lst = list()
  y_lst = list()
  
  # Save eigenvalues and true variances of random effect covariance matrix used in simulations
  eigen_values = matrix(0, nrow = length(sims), ncol = max(r_true_vec))
  rownames(eigen_values) = sim_type
  vars_true = matrix(0, nrow = length(sims), ncol = q)
  rownames(vars_true) = sim_type
  
  # Average values for all fixed effects selected in best model as well as
  ## average variances selected in best models
  beta_all_avg = matrix(0, nrow = length(sims), ncol = p_tot)
  rownames(beta_all_avg) = sim_type
  colnames(beta_all_avg) = str_c("Beta",1:p_tot)
  
  vars_all_avg =  matrix(0, nrow = length(sims), ncol = p_tot+1)
  rownames(vars_all_avg) = sim_type
  colnames(vars_all_avg) = str_c("Var",0:p_tot)
  
  # Save bias information
  abs_dev_lst = list()
  # Save time information
  time_lst = list()
  
  # Save death rate and censoring rate information
  event_rate_lst = list()
  # Save follow-up time information
  # fwp_time = list()
  
  # Save convergence information
  EM_iter_lst = list()
  
  for(m in 1:length(sims)){
    
    files = list.files(path = str_c(path, sims[m]), full.names = T)
    if(trace == 1){
      cat("completed replicates: ", length(files), "\n")
    }
    
    if(length(r_true_vec) == 1){
      r_true = r_true_vec
    }else if(length(r_true_vec) > 1){
      r_true = r_true_vec[m]
    }
    
    time_mat = numeric(length(files))
    r_est = numeric(length(files))
    
    beta_mat = matrix(0, nrow = length(files), ncol = p_tot)
    abs_dev = numeric(length(files)) # Mean absolute deviation
    vars_mat = matrix(0, nrow = length(files), ncol = p_tot+1)
    # presc_mat = matrix(0, nrow = length(files), ncol = q)
    EM_iter = NULL
    opt_res = NULL
    select_res = list()
    r_input = numeric(length(files))
    y_summary = matrix(NA, nrow = length(files), ncol = 6)
    event_rate_mat = matrix(NA, nrow = length(files), ncol = 2)
    colnames(event_rate_mat) = c("death_rate","censor_rate")
    # fwp_summary =  matrix(NA, nrow = length(files, ncol = 6))
    
    for(f in 1:length(files)){
      # load output list object
      load(files[f])
      beta_mat[f,] = output$coef_mat
      beta_vec = output$coef_mat[,c(1:p_true)]
      beta_non0 = beta_vec[which(beta_vec != 0)]
      if(length(beta_non0) >= 1){
        abs_dev[f] = mean(abs(beta_non0 - rep(beta_true[m],times=length(beta_non0))))
      }else{
        abs_dev[f] = NA
      }
      vars_mat[f,] = output$vars_mat
      # presc_mat[f,] = output$presc_mat
      opt_res = rbind(opt_res, output$results_optim[[1]])
      EM_iter = rbind(EM_iter,output$results_all[,"EM_iter"])
      time_mat[f] = output$time_mat[1,3]
      select_res[[f]] = output$results_all
      if(!is.null(output$r_est)){
        r_est[f] = output$r_est
      }
      y_tmp = output$y_times
      if(!is.null(y_tmp)){
        y_summary[f,] = summary(y_tmp)[1:6]
        if(f == 1){
          colnames(y_summary) = names(summary(y_tmp))[1:6]
        }
      }
      if(!is.null(output$death_rate)){
        event_rate_mat[f,1] = output$death_rate
        event_rate_mat[f,2] = 1 - output$death_rate
      }
      
      if(!is.null(output$r_input)){
        r_input[f] = output$r_input
      }else{
        r_input[f] = NA
      }
      
    } # End f for loop
    
    abs_dev_lst[[m]] = abs_dev
    
    EM_iter_lst[[m]] = EM_iter
    
    time_lst[[m]] = time_mat
    
    event_rate_lst[[m]] = event_rate_mat
    
    # True and false positives - fixed effects, random effects,
    # "TP PreSc","FP Presc","Med Presc"
    out = numeric(length(labs))
    names(out) = labs
    
    
    for(i in 1:p_true){
      out[i] = mean(beta_mat[which(beta_mat[,i] != 0),i])
    }
    for(i in 1:q){
      # Note: will be NA if never != 0
      beta_all_avg[m,i] = mean(beta_mat[which(beta_mat[,i] != 0),i])
      vars_all_avg[m,i] = mean(vars_mat[which(vars_mat[,i] != 0),i])
    }
    out[p_true+1] = sum(beta_mat[,c(1:p_true)] != 0) / length(files) * (1 / (p_true) * 100)
    out[p_true+2] = sum(beta_mat[,-c(1:p_true)] != 0) / length(files) * (1 / (p_tot - p_true) * 100)
    out[p_true+3] = sum(vars_mat[,c(2:q)] != 0) / length(files) * (1 / (q - 1) * 100)
    out[p_true+4] = sum(vars_mat[,-c(1:q)] != 0) / length(files)*(1 / (p_tot-q+1) * 100)
    out[p_true+5] = median(time_mat / 3600)
    out[p_true+6] = mean(abs_dev, na.rm = TRUE)
    
    out_mat[m,] = out
    
    # summary of r estimation using pseudo random effects estimates
    rest_pseudo[m,1] = sum(r_est < r_true) / length(files) * 100
    rest_pseudo[m,2] = sum(r_est == r_true) / length(files) * 100
    rest_pseudo[m,3] = sum(r_est > r_true) / length(files) * 100
    
    rest_avg_pseudo[m] = mean(r_est)
    
    rest_lst[[m]] = r_est
    vars_lst[[m]] = vars_mat
    beta_lst[[m]] = beta_mat
    r_input_lst[[m]] = r_input
    
    y_lst[[m]] = y_summary
    
    eigen_values[m,] = output$eigen_vals[1:max(r_true_vec)]
    B_true = output$B_true
    vars_true[m,] = diag(B_true %*% t(B_true))
  }
  
  
  
  out_mat = round(out_mat, digits = 2)
  rest_pseudo = round(rest_pseudo, digits = 2)
  rest_avg_pseudo = round(rest_avg_pseudo, 2)
  
  vars_val = round(diag(output$B_true %*% t(output$B_true)), digits = 2)
  
  return(list(out_mat = out_mat,
              EM_iter = EM_iter_lst,
              rest_pseudo = rest_pseudo,
              rest_avg_pseudo = rest_avg_pseudo,
              vars_val = vars_val, rest_lst = rest_lst,
              vars_lst = vars_lst, beta_lst = beta_lst,
              beta_all_avg = beta_all_avg, vars_all_avg = vars_all_avg,
              r_input = r_input_lst,
              y_lst = y_lst,
              eigen_values = eigen_values,
              vars_true = vars_true,
              time_lst = time_lst,
              abs_dev_lst = abs_dev_lst,
              event_rate_lst = event_rate_lst))
  
}

####################################################################################################################
# Variable selection for piecewise exponential mixed models, p=100
####################################################################################################################


# Path to full simulation results
path = sprintf("%s/phmmPen_FA_p100_02/",path_sims)
# Sub-folders with individual simulation results
sims = str_c("sim",1:16)
# Description of individual simulation results
sim_type = str_c("Beta_",rep(c("0.5","1.0"),each=8),
                 "_K_",rep(rep(c(5,10),each=4),times=2),
                 "_B_",rep(rep(c("Small","Moder"),each=2),times=4),
                 "_r_",rep(c("Truth","GR"),times=8))
res = extract_r(path = path, sims = sims, sim_type = sim_type,
                p_tot = 100, q = 6, p_true = 5, r_true_vec = rep(3,length(sims)),
                trace = 0, beta_true = rep(c(0.5,1.0),each=8))
save(res, file = sprintf("%s/phmmPen_FA_p100.RData",path_output))


####################################################################################################################
# Variable selection for piecewise exponential mixed models, p=500
####################################################################################################################

# Path to full simulation results
path = sprintf("%s/phmmPen_FA_p500/",path_sims)
# Sub-folders with individual simulation results
sims = str_c("sim",1:8)
# Description of individual simulation results
sim_type = str_c("Beta_",rep(c("0.5","1.0"),each=4),
                 "_K_",rep(rep(c(5,10),each=2),times=2),
                 "_B_",rep(c("Small","Moder"),times=4),
                 "_r_", "GR")
res = extract_r(path = path, sims = sims, sim_type = sim_type,
                p_tot = 500, q = 6, p_true = 5, r_true_vec = rep(3,length(sims)),
                trace = 0, beta_true = rep(c(0.5,1.0),each=4))
save(res, file = sprintf("%s/phmmPen_FA_p500.RData",path_output))

####################################################################################################################
# 
####################################################################################################################