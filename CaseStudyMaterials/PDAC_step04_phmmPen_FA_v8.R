# Perform variable selection on Piecewise Exponential data using the glmmPen_FA method 

library(glmmPen) # Version 1.5.4.2
library(stringr)
library(survival)

# Arrays 1-8
array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Path to home folder
header_nas = "~/"
# Path to additional scratch space
header_work = "/work/users/h/h/hheiling/"
# Name of subfolders to save results
sub_folder = "Step04_Results_v8/"
# prefix: folder to save results
## prefix0: location to save main output files
prefix = str_c(header_nas, "Project3/PDAC_phmmPen_FA_TSP/",sub_folder)
if(!dir.exists(prefix)){dir.create(prefix, recursive = TRUE)}
## prefix0_post: location to save files of posterior draws needed for BIC-ICQ calculations
prefix_BICq = str_c(header_work, "Project3/",sub_folder)
if(!dir.exists(prefix_BICq)){dir.create(prefix_BICq, recursive = TRUE)}

# Load PDAC survival data - PDAC_surv object (data.frame)
load(sprintf("%sProject3/PDAC_phmmPen_FA_TSP/PDAC_surv_TSP.RData",header_nas))
## extract relevant values
y_times = PDAC_surv$time * 30 # Convert from months to days
y_status = PDAC_surv$event
study = PDAC_surv$study
# data checks
class(y_times) # Should be numeric
class(y_status) # Should be numeric
class(study) # Should be a factor
# TSP covariate matrix
X = as.matrix(PDAC_surv[,-c(1:4)]) # Remove sampID, study, time, event; leave all TSP covariates
print(colnames(X)) # Check: only TSP covariates
print(dim(X)) # Should have 168 columns (TSP covariates)

cut_num_val = c(8)
var_start_val = c(0.02)
conv_CD_val = c(10^-3)
alpha_val = seq(from = 0.7, to = 1.0, by = 0.1)
est_type = c("GR","truth_3") # Compare GR estimate vs guesses of 3 or 4 to see if obvious differences in results
seq_min_val = c(0.10)
# All simulation combinations
combos = expand.grid(cut_num_val, var_start_val, conv_CD_val, est_type, alpha_val,seq_min_val)
colnames(combos) = c("cut_num_val","var_start_val","conv_CD_val","est_type","alpha_val","seq_min_val")
combos

batch = array_val

print(batch)
print(combos[batch,])

if(batch <= 9){
  batch_label = str_c("00",batch)
}else if(batch <= 99){
  batch_label = str_c("0",batch)
}else{
  batch_label = as.character(batch)
}

# Extract relevant values from combinations
cut_num = combos[batch,"cut_num_val"]
var_start = combos[batch,"var_start_val"]
conv_CD = combos[batch,"conv_CD_val"]
est_type = combos[batch,"est_type"]
alpha = combos[batch,"alpha_val"]
seq_min = combos[batch,"seq_min_val"]


if(est_type == "GR"){
  r_input = NULL
  r_max = 5
  sample = FALSE
}else if(est_type == "truth_3"){
  r_input = 3
  r_max = NULL
  sample = FALSE
}


start = proc.time()

set.seed(2023)
fit = phmmPen_FA(formula = Surv(y_times, y_status) ~ X + (X | study), 
                 alpha = alpha,
                optim_options = optimControl(nMC_start = 100, # Starting number of posterior draws per E-step
                                             nMC_burnin = 100, # Burn-in number of posterior draws per E-step
                                             nMC_max = 500, # Maximum number of posterior draws per E-step
                                             convEM_type = "AvgEuclid1", # Type of convergence to use for EM algorithm
                                             conv_EM = 0.0015, # Value EM convergence criteria needs to meet
                                             maxitEM = 25, # Maximum number of EM iterations per model fit
                                             conv_CD = conv_CD, # Value of convergence criteria for M-step
                                             t = 2, mcc = 2, # To evaluate convergence, compare most recent coefficient vector with vector t=2 iterations back, and meet this criteria mcc = 2 times
                                             B_init_type = "deterministic",
                                             var_start = var_start,
                                             step_size = 1.0,
                                             var_restrictions = "fixef",
                                             standardization = TRUE),
                tuning_options = selectControl(nlambda = 10, lambda.min = seq_min, lambda.min.presc = seq_min,
                                               search = "abbrev",
                                               BIC_option = "BICq"),
                survival_options = survivalControl(cut_num = cut_num,
                                                   interval_type = "equal"),
                r_estimation = rControl(r = r_input, r_max = r_max, size = 25, sample = sample),
                BICq_posterior = sprintf("%s/BICq_PDAC_Batch_%s", # Location to save posterior draws needed to calculate BIC-ICQ posterior
                                         prefix_BICq, batch_label),
                trace = 0, progress = TRUE)
              

end = proc.time()


# Save relevant output

methods = c("phmmPen")

lhaz_tot = cut_num

# lhaz estimates
lhaz_vals = fixef(fit)[1:lhaz_tot]

# Fixed effects
coef_vals = fixef(fit)[-c(1:lhaz_tot)]

# Random effects variance
vars_vals = diag(fit$sigma)
names(vars_vals) = rownames(fit$sigma)

# Timing
time_mat = matrix(0, nrow = 1, ncol = 3)
time_mat[1,] = (end - start)[1:3]
rownames(time_mat) = methods

output = list(coef_vals = coef_vals, 
              lhaz_vals = lhaz_vals,
              vars_vals = vars_vals,
              time_mat = time_mat,
              lhaz_tot = lhaz_tot,
              results_all = fit$results_all,
              r_est = fit$r_estimation$r, r_input = r_input,
              fit = fit) 


save(output, file = sprintf("%s/Output_%s.RData", prefix, batch_label))

file.remove(sprintf("%s/BICq_PDAC_Batch_%s.bin", prefix_BICq, batch_label))
file.remove(sprintf("%s/BICq_PDAC_Batch_%s.desc", prefix_BICq, batch_label))

################################################################################################

print(gc(full = TRUE))

q(save="no")

################################################################################################
