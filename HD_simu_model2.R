# Simulation function for HD-Aug

# sigma = hardcoded to 1
# lambda = hardcoded to c(5, 4.8, 4.6, ..., 3)
# d = 11
# maximum d we accept hardcoded to 30

#setwd("C:/Users/uradojicic/Desktop/HD_simus")

# For M-P quantile function
library(RMTstat)


# For PCAaug
library(ICtest)


# For data manipulation
library(tidyverse)




# Eigenvalue debiasing function
f_lambda <- function(tau, sigma2, pdivn){
  temp <- sigma2*(1 + pdivn) - tau
  (-1*temp + sqrt(max(c(temp^2 - 4*pdivn*sigma2^2, 0))))/2
}



# The HD-aug estimator
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
  
  which.min(diff(crit[1:min(p, 30)]))
}



# Schott's high-dimensional test
PCAschott_own <- function(x, k){
  n <- nrow(x)
  p <- ncol(x)
  S <- cov(x)
  eval_S <- eigen(S, symmetric = TRUE)$values[(k + 1):p]
  test_stat <- ((n - k)*(mean(eval_S^2)/(mean(eval_S))^2 - 1) - (p - k) - 1)/2
  2*(1 - pnorm(abs(test_stat)))
}



# Forward dimension estimation based on Schott's subsphericity test
schott_test <- function(x){
  p <- ncol(x)
  cond <- TRUE
  k <- -1
  while(cond){
    k <- k + 1
    if(PCAschott_own(x, k) >= 0.05){
      cond <- FALSE
    }
    if(k == min(30, p - 2)){
      cond <- FALSE
    }
  }
  k
}



# Single simulation run
single_simu <- function(n, gp, gr,seed=sample(1:10000,1)){
  p <- floor(n*gp)
  r <- floor(n*gr)
  
  sigma2 <- 1; 
  d <- 11;
  p<-max(p,d);
  
  lambdas <- seq(5,3,length.out=d)
  # Generate data

  set.seed(seed)
  x <- sqrt(sigma2)*2*(matrix(rbinom(n*p, 1, 0.5), n, p) - 0.5)
  for (i in 1:d) {
    x[, i] <- x[, i]*sqrt(1 + lambdas[i]/sigma2)
  }
  
  # Estimate sigma
  val <- eigen(cov(x))$values[1:min(n - 1, p)]
  if(n - 1 >= p){
    sigma2_est <- median(val)/qmp(1/2, ndf = n - 1, pdim = p, var = sigma2)
  }
  if(n - 1 < p){
    sigma2_est <- median(val)/(qmp(1/2, ndf = p, pdim = n - 1, var = sigma2)*p/(n - 1))
  }
  
  
  # Luo & Li
  est_1 <- PCAaug(x, noise = "known", naug = r, nrep = 1, sigma2 = sigma2)$k
  est_2 <- PCAaug(x, noise = "known", naug = r, nrep = 1, sigma2 = sigma2_est)$k
  
  # Us
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


# Results of one simulation
# 
# n <- 200
# gp <- 0.8
# gr <- 3
# 
# single_simu(n, gp, gr)
# 


# 1000 replicates of results for all combinations of parameters
# Paralelize by the chunks of 10 repetitions

n_list <- c(100,200,500,1000)
gp_list <- c(0.05, 0.2, 0.5, 1, 1.5)
gr_list <- c(0.05, 0.2, 0.5, 1, 1.5, 2.5,5)


# Define the number of cores to use
num_cores <- 100  # Leave 1 core free


library(doParallel)
library(foreach)

# Register parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Run the simulations in parallel
results <- foreach(m = 1:100, .combine = rbind, .packages = c("RMTstat", "ICtest", "tidyverse")) %dopar% {
  
  # Set seed for reproducibility
  set.seed(m)
  
  # Generate seeds for the current set of simulations
  seeds <- sample(1:10000, 10)
  
  # Generate parameter grid
  par_grid <- expand.grid(n = n_list, gp = gp_list, gr = gr_list, seeds = seeds)
  par_n <- nrow(par_grid)
  
  # Run the simulation for each parameter combination in this chunk
  temp <- lapply(1:par_n, function(i) {
    single_simu(par_grid$n[i], par_grid$gp[i], par_grid$gr[i], par_grid$seeds[i])
  })
  
  # Combine results from this chunk
  bind_rows(temp)
}

# Write final results to a single CSV
write.csv2(results, file = "HDSimu_model2.csv")