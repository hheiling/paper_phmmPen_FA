# compile case study results

# Define directory where all case study results are kept
path_sims = "~/PDAC_phmmPen_FA_TSP"
# Define directory where to store simulation results
path_output = "~/Paper_Results"

######################################################################################################################
# phmmPen_FA case study results
######################################################################################################################

# Load results
path = sprintf("%s/Step04_Results_v8/",path_sims)
files = list.files(path = path, pattern = ".RData", full.names = TRUE)
labels = str_c("Alpha_",rep(seq(from = 0.7, to = 1.0, by = 0.1), each = 2),
               "_r_",rep(c("GR_est","3"),times=4))


res = list()
for(i in 1:length(files)){
  # load output list object
  load(file = files[i])
  output_tmp = output
  output_tmp$coef_vals = output$coef_mat[1,]
  output_tmp$coef_mat = NULL
  output_tmp$time_mat = output_tmp$time_mat[1,,drop=FALSE]
  output_tmp$c_index = output_tmp$cindex_lst[[1]]
  res[[labels[i]]] = output_tmp
}
save(res, file = "Paper_Results/PDAC_Comparison_Results_revision.RData")

######################################################################################################################

######################################################################################################################