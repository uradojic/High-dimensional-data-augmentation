# Implementation of the HDPA-method from
# "Dimension estimation in PCA model using high-dimensional data augmentation"
# By U. Radojicic and J. Virta

# In the below code:
# sigma = hardcoded to 1
# lambda = hardcoded to c(5, 4, 3)
# maximum estimated d we accept hardcoded to 20



# For M-P quantile function
library(RMTstat)

# For PA
library(ICtest)

# For data manipulation
library(tidyverse)




# Eigenvalue debiasing function
#
# tau = the eigenvalue
# sigma2 = noise variance
# pdivn = the ratio of p to n
f_lambda <- function(tau, sigma2, pdivn){
  temp <- sigma2*(1 + pdivn) - tau
  (-1*temp + sqrt(max(c(temp^2 - 4*pdivn*sigma2^2, 0))))/2
}



# The HDPA estimator
#
# x = the data matrix
# r = the amount of augmentations
# sigma2 = the noise variance (true value or an estimate)
our_estimator <- function(x, r, sigma2){
  n <- nrow(x)
  p <- ncol(x)
  
  S0 <- cov(x)
  eval_S0 <- eigen(S0)$values
  
  est_lambda <- sapply(eval_S0, FUN = function(s) f_lambda(s, sigma2, p/n))
  
  
  # Augment and get norms
  
  xa <- sigma2*matrix(rnorm(n*r), n, r)
  S <- cov(cbind(x, xa))
  eig_S <- eigen(S)
  U <- eig_S$vectors[(p + 1):(p + r), 1:p]
  sqnorm_vec <- colSums(U^2)
  
  
  # Final criterion
  crit <- (sqnorm_vec*est_lambda)*(est_lambda + (p/n + r/n)*sigma2)/(est_lambda + sigma2)
  
  which.min(diff(crit[1:min(p, 20)]))
}



# Schott's high-dimensional subsphericity test
# 
# x = the data matrix
# k = the dimension for which the test is carried out
PCAschott_own <- function(x, k){
  n <- nrow(x)
  p <- ncol(x)
  S <- cov(x)
  eval_S <- eigen(S, symmetric = TRUE)$values[(k + 1):p]
  test_stat <- ((n - k)*(mean(eval_S^2)/(mean(eval_S))^2 - 1) - (p - k) - 1)/2
  2*(1 - pnorm(abs(test_stat)))
}



# Forward dimension estimation based on Schott's subsphericity test
# 
# x = the data matrix
schott_test <- function(x){
  p <- ncol(x)
  cond <- TRUE
  k <- -1
  while(cond){
    k <- k + 1
    if(PCAschott_own(x, k) >= 0.05){
      cond <- FALSE
    }
    if(k == min(20, p - 2)){
      cond <- FALSE
    }
  }
  k
}



# Single simulation run
#
# n = sample size
# gp = the ration of p to n
# gr = the ration of r to n
single_simu <- function(n, gp, gr){
  p <- floor(n*gp)
  r <- floor(n*gr)
  
  sigma2 <- 1
  lambda_1 <- 5
  lambda_2 <- 4
  lambda_3 <- 3
  
  # Generate data
  x <- matrix(rnorm(n*p, 0, sqrt(sigma2)), n, p)
  x[, 1] <- x[, 1]*sqrt(1 + lambda_1/sigma2)
  x[, 2] <- x[, 2]*sqrt(1 + lambda_2/sigma2)
  x[, 3] <- x[, 3]*sqrt(1 + lambda_3/sigma2)
  
  # Estimate sigma
  val <- eigen(cov(x))$values[1:min(n - 1, p)]
  if(n - 1 >= p){
    sigma2_est <- median(val)/qmp(1/2, ndf = n - 1, pdim = p, var = sigma2)
  }
  if(n - 1 < p){
    sigma2_est <- median(val)/(qmp(1/2, ndf = p, pdim = n - 1, var = sigma2)*p/(n - 1))
  }
  
  
  # PA by Luo & Li
  est_1 <- PCAaug(x, noise = "known", naug = r, nrep = 1, sigma2 = sigma2)$k
  est_2 <- PCAaug(x, noise = "known", naug = r, nrep = 1, sigma2 = sigma2_est)$k
  
  # HDPA
  est_3 <- our_estimator(x, r, sigma2)
  est_4 <- our_estimator(x, r, sigma2_est)
  
  # Schott
  est_5 <- tryCatch(schott_test(x), error = function(e) NA)
  
  
  data.frame(n = rep(n, 5),
             gp = rep(gp, 5),
             gr = rep(gr, 5),
             est = c(est_1, est_2, est_3, est_4, est_5),
             method = c(1, 2, 3, 4, 5))
  
}







# 10 replicates of results for all combinations of parameters

reps <- 10

n_list <- c(100, 200)
gp_list <- c(0.05, 0.2, 0.5, 1, 1.5)
gr_list <- c(0.05, 0.2, 0.5, 1, 1.5)

par_grid <- expand.grid(n = n_list, gp = gp_list, gr = gr_list)

par_n <- nrow(par_grid)

res <- NULL

for(i in 1:reps){
  temp <- lapply(1:par_n, function(i) single_simu(par_grid$n[i], par_grid$gp[i], par_grid$gr[i]))
  temp_2 <- bind_rows(temp)
  res <- bind_rows(res, temp_2)
  print(i)
}


