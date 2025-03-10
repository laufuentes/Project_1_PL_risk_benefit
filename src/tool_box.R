set.seed(2025)
setwd("~/Project_1_PL_risk_benefit")


option_det <- function(string, split_char){
  res<- strsplit(string, split = split_char)
  return(res[[1]])
}