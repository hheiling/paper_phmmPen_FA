# Perform variable selection on Piecewise Exponential data using the glmmPen_FA method 

library(glmmPen) # Version 1.5.4.1
library(stringr)
library(survival)
library(ncvreg)

# Arrays 1-800 - 8 simulation types, 100 replicates per simulation type
array_val <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# Path to home folder
header_nas = "~/"
# Path to additional scratch space
header_work = "/work/users/h/h/hheiling/"
# Name of subfolders to save results
sub_folder = "ncvsurv_01/"
# prefix: folder to save results
## prefix0: location to save main output files
prefix0 = str_c(header_nas, "Project3/Simulation_Results/",sub_folder)
if(!dir.exists(prefix0)){dir.create(prefix0, recursive = TRUE)}

# 100 simulation replicates needed
sim_total = 100

p_tot = 100

# p = total number of predictors - no intercept in coxph model
p_val = c(5)
beta_val = c(0.5,1.0)
cut_num_val = c(8)
B_type = c("small","moderate")
var_start_val = c(0.020)
conv_CD_val = c(10^-3)
K_val = c(5,10)
# All simulation combinations
combos = expand.grid(p_val, var_start_val, conv_CD_val, cut_num_val, 
                     B_type, K_val, beta_val)
colnames(combos) = c("p_val", "var_start_val", "conv_CD_val", "cut_num_val",
                    "B_type", "K_val", "beta_val")
combos

r = 3

# Place each simulation set-up into separate sub-folder
if((array_val %% sim_total) != 0){
  batch = array_val %% sim_total
  sim_num = array_val %/% sim_total + 1
}else{
  batch = sim_total # 100
  sim_num = array_val %/% sim_total 
}

print(batch)
print(sim_num)
print(combos[sim_num,])

# prefix and prefix_BICq: folder to save results and posterior needed for BIC-ICQ, respectively
prefix = str_c(prefix0, "sim", sim_num) 
if(!dir.exists(prefix)){dir.create(prefix, recursive = TRUE)}

if(batch <= 9){
  batch_label = str_c("00",batch)
}else if(batch <= 99){
  batch_label = str_c("0",batch)
}else{
  batch_label = as.character(batch)
}

# Extract relevant values from combinations
p = combos[sim_num,"p_val"]
beta = combos[sim_num,"beta_val"]
cut_num = combos[sim_num,"cut_num_val"]
B_size = combos[sim_num,"B_type"]
var_start = combos[sim_num,"var_start_val"]
conv_CD = combos[sim_num,"conv_CD_val"]
K = combos[sim_num,"K_val"]

q = p + 1
B0 = cbind(rep(1,q),
           c(rep(c(-1, 1), each=3)),
           c(rep(c(-1, 0, 1), times=2)) )
if(B_size == "large"){
  B = B0
  cov_mat = B %*% t(B)
  print(round(diag(cov_mat), 2))
  SVD = svd(cov_mat)
  round(SVD$d,2)
  eigen_vals = SVD$d
}else if(B_size == "moderate"){
  B = B0 * 0.75
  cov_mat = B %*% t(B)
  print(round(diag(cov_mat), 2))
  SVD = svd(cov_mat)
  round(SVD$d,2)
  eigen_vals = SVD$d
}else if(B_size == "small"){
  B = B0 * 0.50
  cov_mat = B %*% t(B)
  print(round(diag(cov_mat), 2))
  SVD = svd(cov_mat)
  round(SVD$d,2)
  eigen_vals = SVD$d
}

# Simulation of data
# Assume true model only has q-1 non-zero fixed and random effect covariates (ignoring intercept)
# Number total covariates p = 100 (ignore intercept)

N = 1000

set.seed(2022) 
seeds = sample(1000:9999, size = sim_total, replace = FALSE)


cut_points_input = seq(from = 0, to = 2, by = 0.5)
lhaz_vals = c(-1.5,1.0,2.7,3.7,6.8)
dat = sim.data.piecewise.exp(n = N, ptot = 100, pnonzero = p, nstudies = K,
                             sd_raneff = 0, B = B, r = r,
                             cens_type = "unif", cens_max = 5,
                             lhaz_vals = lhaz_vals, cut_points = cut_points_input,
                             seed = seeds[batch], imbalance = 0, beta = rep(beta,p), 
                             pnonzerovar = 0, sd_x = 1.0)

y_times = dat$y_times * 365 # times - convert from years to days
y_status = dat$y_status # event indicator
X = dat$X
Z = dat$Z
group = dat$group 

start = proc.time()

set.seed(seeds[batch])
# All other arguments not specified: default cv.ncvsurv arguments (nfolds = 10, ...)
fit = cv.ncvsurv(X = X, y = Surv(y_times, y_status), penalty = "MCP",
                 lambda.min = 0.05, # Same as phmmPen_FA procedure
                 nlambda = 100) # Default ncvreg length of sequences

end = proc.time()

# Save relevant output

methods = c("ncvsurv")

# Fixed effects
coef_mat = matrix(coef(fit), nrow = 1)
rownames(coef_mat) = methods

# Timing
time_mat = matrix(0, nrow = 1, ncol = 3)
time_mat[1,] = (end - start)[1:3]
rownames(time_mat) = methods[1]

output = list(coef_mat = coef_mat, 
              time_mat = time_mat,
              death_rate = dat$death_rate,
              y_times = dat$y_times,
              y_status = dat$y_status) 


save(output, file = sprintf("%s/Output_%s.RData", prefix, batch_label))

file.remove(sprintf("%s/BICq_proj3_Batch_%s.bin", prefix_BICq, batch_label))
file.remove(sprintf("%s/BICq_proj3_Batch_%s.desc", prefix_BICq, batch_label))

################################################################################################

print(gc(full = TRUE))

q(save="no")

################################################################################################
