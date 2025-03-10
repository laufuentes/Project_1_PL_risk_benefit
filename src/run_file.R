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

source("src/estimators.R")
source("src/objective_functions.R")
source("src/optim_functions.R")
source("src/synthetic_data.R")

n <- 1e3
setting <- "Other_2"
option <- option_det(setting, "_")
centered <- TRUE
epsilon <- 0.03
Jfold<- 5
technique <- "cv.glmnet"
B <- seq(0.05, 2, 0.01)
Lambda <- seq(0,15, 0.1)
alpha <- 0.1

exp <- data_gen(n, option)
df <- exp[[2]]

X <- df %>% select(starts_with("X")) %>% as.matrix()

c_data <- expand.grid(m=c(1,2,3),a=c(0,1))
c_data$c <- c(0,-1,-1, 1,1,1)

gdot1 <- function(X){
    rep(0, dim(as.matrix(X))[1])
}

gdot2 <-function(X){
    rep(1, dim(as.matrix(X))[1])
}

gdot3 <- function(X){
    as.matrix(rep(1,dim(as.matrix(X))[1]))
}

gdot_functs <- list(gdot1, gdot2, gdot3)

Hm<- function(a,x,m,gdot_functs,c_data, mu_nj){
    gm_dot <- gdot_functs[[m]]
    return(c_data$c[which(c_data$m==m & c_data$a==a)]*gm_dot(mu_nj(a,x)))
}

h1 <- function(a,x){
  Hm(a,x,m=1,gdot_functs,c_data, mu_nj)
}
h2 <- function(a,x){
  Hm(a,x,m=2,gdot_functs,c_data, mu_nj)
}

h_functs <- list(h1, h2)



