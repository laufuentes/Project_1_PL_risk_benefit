# load_results.R
B <- readRDS(file.path("MC","results","data","MC_iter.rds"))
L <- length(readRDS(file.path("MC","results","data","Lambda.rds")))

 # Directory where individual result files are saved
 res_dir <- file.path("MC","results", "oracular","individual_results")

extract_j <- function(filename) {
  as.integer(sub(".*res_\\d+_(\\d+)\\.rds$", "\\1", filename))
}

for (b in 1:B){
   pattern <- sprintf("^res_%d_\\d+\\.rds$", b)
   res_files <- list.files(path = res_dir, pattern = pattern, full.names = TRUE)
    res_files <- res_files[order(sapply(res_files, extract_j))]
    all_results <- lapply(res_files, readRDS)
    results<- do.call(rbind, lapply(all_results, function(x) {
    data.frame(
        lambda = x[[1]],
        beta = x[[2]],
        optimal_x = I(x[[3]]),
        risk = x[[4]],
        constraint = x[[5]],
        obj = x[[6]],
        policy_value = x[[7]])}))
    idx_opt <- which(results$lambda==min(results$lambda[results$constraint<=0]))
    saveRDS(results$lambda[[idx_opt]], file = file.path("MC","results","oracular", "lambda_opt",paste0("opt_lambda_",b,".rds")))
    file.copy(from = file.path("MC","results","oracular","theta_opt",paste0("res_",b,"_",idx_opt,".rds")), to = file.path("MC","results","oracular","theta_opt","opt",paste0("theta_opt_",b,".rds")))
    write.csv(results, file = file.path("MC","results","oracular","group_by_ind",paste0("individual_",b,".rds")), row.names = FALSE)
}

for (l in 1:L){
  # List all result files
  pattern <- sprintf("^res_\\d+_%d\\.rds$", l)
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
    write.csv(results, file = file.path("MC","results","oracular", "group_by_lambda",paste0("lambda_",l,".csv")), row.names = FALSE)
}


 # Load all results into a list
 
 # Optionally, combine into a data.frame if they are compatible (e.g., all named lists or rows)
 # This step depends on what `res` actually is in your process
 # If each `res` is a named list or a data.frame row:

 
 # Save combined results
