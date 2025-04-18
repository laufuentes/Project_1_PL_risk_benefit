# load_results.R

B <- readRDS(file.path("MC","results","data","MC_iter.rds"))
Lambda <- readRDS(file.path("MC","results","data","Lambda.rds"))
L <- length(Lambda)

 # Directory where individual result files are saved
 res_dir <- file.path("MC","results", "estimated","individual_results_per_fold")

extract_j <- function(filename) {
  as.integer(sub(".*res_\\d+_(\\d+)\\.rds$", "\\1", filename))
}

for (b in 1:3){
  for (l in 1:L){
    pattern <- sprintf("^res_%d_%d_.*\\.rds$", b, l)
    res_files <- list.files(path = res_dir, pattern = pattern, full.names = TRUE)
    all_results <- lapply(res_files, readRDS)
    results<- do.call(rbind, lapply(all_results, function(x) {
    data.frame(
        lambda = x[[1]],
        beta = x[[2]],
        risk = x[[4]],
        constraint = x[[5]],
        obj = x[[6]],
        policy_value = x[[7]])}))

    result_T <- results %>%
        summarise(
        obj = mean(obj),
        risk = mean(risk),
        constraint = mean(constraint),
        policy_value = mean(policy_value),
        .groups = "drop")%>% mutate(lambda=Lambda[l])
    saveRDS(result_T, file = file.path("MC","results","estimated","individual_results",paste0("result_",b,"_",l,".rds")))      
  }
}

res_dir_ind <- file.path("MC","results", "estimated","individual_results")
extract_j <- function(filename) {
  as.integer(sub(".*result_\\d+_(\\d+)\\.rds$", "\\1", filename))
}


for (b in 1:B){
   pattern <- sprintf("result_%d_\\d+\\.rds$", b)
   res_files_ind <- list.files(path = res_dir_ind, pattern = pattern, full.names = TRUE)
    res_files_ind <- res_files_ind[order(sapply(res_files_ind, extract_j))]
    all_results <- lapply(res_files_ind, readRDS)
    results<- do.call(rbind, all_results)
    idx_opt <- which(results$lambda==min(results$lambda[results$constraint<=0]))
    saveRDS(results$lambda[[idx_opt]], file = file.path("MC","results","estimated", "lambda_opt",paste0("opt_lambda_",b,".rds")))
    saveRDS(results, file = file.path("MC","results","estimated", "group_by_ind",paste0("individual_",b,".rds")))
}

for (l in 1:L){
  # List all result files
  pattern <- sprintf("^result_\\d+_%d\\.rds$", l)
  res_files_l <- list.files(path = res_dir_ind, pattern = pattern, full.names = TRUE)
  all_results <- lapply(res_files_l, readRDS)
  results<- do.call(rbind,all_results)
  saveRDS(results, file = file.path("MC","results","estimated", "group_by_lambda",paste0("lambda_",l,".rds")))
}


 # Load all results into a list
 
 # Optionally, combine into a data.frame if they are compatible (e.g., all named lists or rows)
 # This step depends on what `res` actually is in your process
 # If each `res` is a named list or a data.frame row:

 
 # Save combined results
