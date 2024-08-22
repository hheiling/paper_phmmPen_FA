# Introduction and Installation
This GitHub repository contains the simulation and case study materials used within the paper "Efficient computation of high-dimensional penalized proportional hazard mixed models using the piecewise exponential approximation" by Heiling et al., which describes the phmmPen_FA method. 

The code to run the phmmPen_FA method is currently available in the GitHub respository https://github.com/hheiling/glmmPen and in the glmmPen R package available on CRAN at https://cran.r-project.org/web/packages/glmmPen/index.html. 

Install glmmPen R package from the GitHub repository using the following code:

```
# GitHub repository download:
library(devtools)
install_github("hheiling/glmmPen")

# CRAN package installation: see documentation for glmmPen::phmmPen_FA() function
install.packages("glmmPen")
library(glmmPen)
```

This respository contains the code to load the contents of the simulations and case study results and create the tables in the paper as well as the code to run the simulations and the code to run the case study analyses.

# Replication of Paper Tables and Content: Simulation and Case Study Output

The folder Replication/ contains RData files in the Replication/Paper_Results/ directory which contain summary results for each set of simulations outlined in the paper as well as the summary results for the case study analyses. The code to replicate the tables and content outlined in the paper is given in the "replication.R" code file. 

Note on wording: "Piecewise exponential" model interchangeable with "piecewise constant hazard" model

In order to replicate the tables and content in the paper, run the following R code:

```
# Define path to the Replication/ folder contents
path = "~/paper_phmmPen_FA/Replication"
# Run code
source(sprintf("%s/replication.R",path))
```

# Case Study Sensitivity Analyses

This repository includes additional investigations into the case study output for alternative combinations of the elastic net parameter (values {0.7,0.8,0.9,1.0}) and number of latent factors r (GR estimate of 2 or manually-specified value of 3).

The folder Replication/Paper_Results/ contains the case study summary results for all elastic net and latent factor r parameters listed above in the output object "PDAC_Comparison_Results_revision.RData".

In order to examine output such as the c-index values, times, and graphical and numerical summaries of the selected predictors, please run the following code:

```
# Define path to the Replication/ folder contents
path = "~/paper_phmmPen_FA/Replication"
# Run code
source(sprintf("%s/replication_case_study_sensitivity.R",path))
```

# Running Simulations

The Simulations/ folder contains the code used to run all of the simulations. There are separate files for each set of simulations described in the paper:

* p=100 penalized piecewise exponential mixed effects simulations: "select100_phmmPen_FA_02.R" (described in main paper)
* p=500 penalized piecewise exponential mixed effects simulations: "select500_phmmPen_FA.R" (described in main paper)
* Comparison to 'naive' variable selection using fixed effects only method from ncvreg R package: "select100_ncvsurv_01.R" (described in supplementary materials document)
* phmmPen\_FA variable selection applied to Weibull-simulated survival mixed effects data: "select100_Weibull_01.R" (described in supplementary materials docment)
* Penalized piecewise exponential mixed effects simulations where the number of true fixed effects does not equal the number of true random effects: "select100_alt_num_ranef_01.R" (described in supplementary materials document)
* Penalized piecewise exponential mixed effects simulations where the number of latent factors r is purposefully underestimated: "select100_alt_rval_01.R" (described in supplementary materials document)
* Penalized piecewise exponential mixed effects simulations where alternative values for the number of time intervals J are considered: "select100_alt_J_01.R" (values 5-6) and "select100_alt_J_01B.R" (values 9-10) (described in supplementary materials document)

To run the simulations, each code file should be manually modified to adjust the following arguments:

* prefix0 - path to location where simulation results should be stored
* prefix0_post - path to location where the MCMC posterior samples from the 'full model' (model fit with minimum penalty values) are temporarily stored for use in calculating the BIC-ICQ model selection criteria

Each file specifies the code needed to simulate a single dataset, perform variable selection on that dataset, and save the relevant results. In order to get all of the simulation results specified in the paper, which cover several simulation conditions and 100 replicates per condition, the code needs to be submitted to a computing cluster to run many replicates in parallel. 

Linux code to submit the jobs to a computing cluster:

```
sbatch --array=1-1600 -N 1 -t 36:00:00 --mem=2g -n 1 --output=PHFA_100_%a.out --wrap="R CMD BATCH select100_phmmPen_FA_02.R PHFA_100_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-800 -N 1 -t 80:00:00 --mem=2g -n 1 --output=PHFA_500_%a.out --wrap="R CMD BATCH select500_phmmPen_FA.R PHFA_500_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-800 -N 1 -t 2:00:00 --mem=2g -n 1 --output=ncvsurv_01_%a.out --wrap="R CMD BATCH select100_ncvsurv_01.R ncvsurv_01_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-800 -n 1 -t 36:00:00 --mem=2g -n 1 --output=WEIBULL_01_%a.out --wrap="R CMD BATCH select100_Weibull_01.R WEIBULL_01_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-400 -n 1 -t 36:00:00 --mem=2g -n 1 --output=altnumref_01_%a.out --wrap="R CMD BATCH select100_alt_num_ranef_01.R numranef_01_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-600 -n 1 -t 36:00:00 --mem=2g -n 1 --output=altrval_01_%a.out --wrap="R CMD BATCH select100_alt_rval_01.R rval_01_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-800 -n 1 -t 36:00:00 --mem=2g -n 1 --output=alt_J_01_%a.out --wrap="R CMD BATCH select100_alt_J_01.R alt_J_01_$SLURM_ARRAY_TASK_ID.Rout"

sbatch --array=1-800 -n 1 -t 36:00:00 --mem=2g -n 1 --output=alt_J_01B_%a.out --wrap="R CMD BATCH select100_alt_J_01B.R alt_J_01B_$SLURM_ARRAY_TASK_ID.Rout"
```

Once all of the simulations are run, the code "combine_sim_results.R" can be used to create the RData output files given in Replication/Paper_Results folder. The "path_sim" and "path_output" arguments may need to be manually adjusted in this file.

# Case Study Materials

The CaseStudyMaterials/ folder contains the following items:

* "PDAC_step01_data_prep.R" - this code cleans the data, merges the individual datasets together, and creates all possible top scoring pair (TSP) covariates. The individual data files are not included in this repository at this time, so this code cannot be run directly.
* "step01_TSP.RData" - this RData object contains output from step 01. This output is a data.frame that contains the study id, survival information (time and event/censoring), study membership, and all TSP values for each subject. This object is used in "step02" described below.
* "PDAC_step02_TSP_loglik_parallel.R" - this code loads the "step01_TSP.RData" object and calculates the log-likelihood for each TSP covariate. These log-likelihoods are calculated for a model with the TSP, a random intercept, and a random slope using the coxme::coxme() function.
* "PDAC_step03_TSP_choice.R" - this code reads in the log-likelihood values from step 02 and uses this information to select the best TSP covariates to use in the variable selection procedure. This also outputs a data.frame object titled "PDAC_surv_TSP.RData" that is used in the variable selection procedure.
* "PDAC_surv_TSP.RData" - this RData file contains the output dataset created in the variable selection procedure given in "PDAC_step04_phmmPen_Fa_v8.R" and is used in the case study analyses.
* "PDAC_step04_phmmPen_FA_comparison_02.R" - this code performs variable selection on the PDAC_surv_TSP.RData dataset using the phmmPen_FA algorithm with elastic net penalization.
* "compile_casestudy_results.R" - this code takes the results from the "PDAC_step04_phmmPen_FA_comparison_02.R" procedures and creates the "PDAC_Selection_Results.RData" output file given in Replication/Paper_Results/ folder. The "path_sim" and "path_output" arguments may need to be manually adjusted in this file.

To run the "PDAC_step04_phmmPen_FA_comparison_02.R" procedures, the code files must first have the following arguments manually adjusted:

* prefix - path to location where case study results should be stored
* prefix_BICq - path to location where the MCMC posterior samples from the 'full model' (model fit with minimum penalty values) are temporarily stored for use in calculating the BIC-ICQ model selection criteria
* load(sprintf("%sProject3/PDAC_phmmPen_FA_TSP/PDAC_surv_TSP.RData",header_nas)) - this line may need to be updated for the location where the PDAC_basal.RData dataset is stored 

The code can be submitted to the computing cluster using the commands outlined below:

```
sbatch --array=1-8 -N 1 -t 12:00:00 --mem=2g -n 1 --output=pdac04_%a.out --wrap="R CMD BATCH PDAC_step04_phmmPen_FA_comparison_02.R pdac04_$SLURM_ARRAY_TASK_ID.Rout"
```
