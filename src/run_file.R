library(tidyverse)
library(optimx)
library(lbfgs)
library(boot)
library(graphics)
library(ggplot2)
library(grid)
library(rgenoud)
library(locfit)
library(madness)

source("src/objective_functions.R")
source("src/optim_functions.R")
source("src/synthetic_data.R")

n <- 1e3
setting <- "Other_2"
option <- option_det(setting, "_")
centered <- TRUE
epsilon <- 0.03
beta <- 0.5 
lambda <- 5
alpha <- 0.1

exp <- data_gen(n, option)
df <- exp[[2]]

X <- df %>% select(starts_with("X"))





