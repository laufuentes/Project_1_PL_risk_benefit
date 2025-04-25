# load_results.R
folder <- readRDS(file.path("run_MC_simu","current_scenario.rds"))
B <- readRDS(file.path(folder,"results","data","MC_iter.rds"))
Lambda <- readRDS(file.path(folder,"results","data","Lambda.rds"))
L <- length(Lambda)

source("run_MC_simu/R_files/utils.R")

#################################################
# Process multiple folds for estimation results #
#################################################
res_dir_mfolds <- file.path(folder,"results", "estimated","individual_results_per_fold")
or_res_dir_mfolds <- file.path(folder,"results", "estimated", "oracular","or_individual_results_per_fold")

# Loop through datasets and parameter combinations
for (b in 1:B) {
  for (l in 1:L) {
    # Generate the file pattern for each combination of b and l
    pattern <- sprintf("^res_%d_%d_.*\\.rds$", b, l)
    
    # List the files in the respective directories
    res_files <- list.files(path = res_dir_mfolds, pattern = pattern, full.names = TRUE)
    res_files_or <- list.files(path = or_res_dir_mfolds, pattern = pattern, full.names = TRUE)
    
    # Process the results from both the original and the "or" directories
    results <- process_results(res_files)
    results_or <- process_results(res_files_or)
    
    # Summarize the results
    result_T <- summarize_results(results)
    result_T_or <- summarize_results(results_or)
    
    # Save the summarized results
    write.csv(result_T, file = file.path(folder, "results", "estimated", "individual_results", paste0("result_", b, "_", l, ".csv")))
    write.csv(result_T_or, file = file.path(folder, "results", "estimated", "oracular", "or_individual_results", paste0("result_", b, "_", l, ".csv")))
  }
}

est_res_dir_ind <- file.path(folder,"results", "estimated","individual_results")
or_est_res_dir_ind <- file.path(folder,"results", "estimated", "oracular","or_individual_results")

for (b in 1:B){
   pattern <- sprintf("result_%d_\\d+\\.csv$", b)
   or_res_files_ind <- list.files(path = or_est_res_dir_ind, pattern = pattern, full.names = TRUE)
   or_res_files_ind <- or_res_files_ind[order(sapply(or_res_files_ind, extract_j_csv))]
   or_all_results <- lapply(or_res_files_ind, function(x)read.csv(x, stringsAsFactors=FALSE))
   or_results_i<- do.call(rbind, or_all_results)
   or_idx_opt <- which(or_results_i$lambda==min(or_results_i$lambda[or_results_i$constraint<=0]))
   saveRDS(or_results_i$lambda[[or_idx_opt]], file = file.path(folder,"results","estimated", "oracular", "lambda_opt",paste0("opt_lambda_",b,".rds")))
   write.csv(or_results_i, file = file.path(folder,"results","estimated", "oracular", "group_by_ind",paste0("individual_",b,".csv")))

   res_files_ind <- list.files(path = est_res_dir_ind, pattern = pattern, full.names = TRUE)
   res_files_ind <- res_files_ind[order(sapply(res_files_ind, extract_j_csv))]
   all_results <- lapply(res_files_ind, function(x)read.csv(x, stringsAsFactors=FALSE))
   results_i<- do.call(rbind, all_results)
   idx_opt <- which(results_i$lambda==min(results_i$lambda[results_i$constraint<=0]))
   saveRDS(results_i$lambda[[idx_opt]], file = file.path(folder,"results","estimated", "lambda_opt",paste0("opt_lambda_",b,".rds")))
   write.csv(results_i, file = file.path(folder,"results","estimated", "group_by_ind",paste0("individual_",b,".csv")))
}

for (l in 1:L){
  # List all result files
  pattern <- sprintf("^result_\\d+_%d\\.csv$", l)

  res_files_l <- list.files(path = est_res_dir_ind, pattern = pattern, full.names = TRUE)
  all_results <- lapply(res_files_l, function(x)read.csv(x, stringsAsFactors=FALSE))
  results_l<- do.call(rbind,all_results)
  write.csv(results_l, file = file.path(folder,"results","estimated", "group_by_lambda",paste0("lambda_",l,".csv")))
  
  or_res_files_l <- list.files(path = or_est_res_dir_ind, pattern = pattern, full.names = TRUE)
  or_all_results <- lapply(or_res_files_l, function(x)read.csv(x, stringsAsFactors=FALSE))
  or_results_l<- do.call(rbind,or_all_results)
  write.csv(or_results_l, file = file.path(folder,"results","estimated", "oracular", "group_by_lambda",paste0("lambda_",l,".csv")))

}


##################################
# Process oracular #
##################################

# Directory where individual result files are saved
res_dir_or <- file.path(folder,"results", "oracular","individual_results")

for (b in 1:B){
   pattern <- sprintf("^res_%d_\\d+\\.rds$", b)
   res_files_or <- list.files(path = res_dir_or, pattern = pattern, full.names = TRUE)
   res_files_or <- res_files_or[order(sapply(res_files_or, extract_j))]
   results_or_i <- process_results(res_files_or)
    idx_opt <- which(results_or_i$lambda==min(results_or_i$lambda[results_or_i$constraint<=0]))
    saveRDS(results_or_i$lambda[[idx_opt]], file = file.path(folder,"results","oracular", "lambda_opt",paste0("opt_lambda_",b,".rds")))
    file.copy(from = file.path(folder,"results","oracular","theta_opt",paste0("theta_",b,"_",idx_opt,".rds")), to = file.path(folder,"results","oracular","theta_opt","opt",paste0("theta_opt_",b,".rds")))
    write.csv(results_or_i, file = file.path(folder,"results","oracular","group_by_ind",paste0("individual_",b,".rds")), row.names = FALSE)
}

for (l in 1:L){
  # List all result files
  pattern <- sprintf("^res_\\d+_%d\\.rds$", l)
  res_files_or_l <- list.files(path = res_dir_or, pattern = pattern, full.names = TRUE)
  res_files_or_l <- res_files_or_l[order(sapply(res_files_or_l, extract_i))]
  results_or_l <- process_results(res_files_or_l)
  write.csv(results_or_l, file = file.path(folder,"results","oracular", "group_by_lambda",paste0("lambda_",l,".csv")), row.names = FALSE)
}
