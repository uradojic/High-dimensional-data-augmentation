# High-dimensional-data-augmentation
This project contains supporting R scripts for reproducing results in ADD REF to ARXIV

The goal of the project is to develop a high-dimensional version of a recent dimension estimation procedure proposed in  
Luo, W. and Li, B. (2021) that determines the dimension via the introduction of augmented noise variables into the data.

All code here is written in R and requires the packages RMTstat, ICtest, tidyverse, ggplot2, and mvtnorm.

The main functions in the source file are:

  - f_lambda: eigenvalue debiasing function.

  - our_estimator: HDPA estimator (ADD ref to ARXIV file).

  - PCAschott_own: high-dimensional subsphericity test by Schott, J. R. (2006). 

  - schott_test: forward dimension estimation based on subsphericity test by Schott, J. R. (2006).

In addition, the project contains two scripts for reproducing the simulation results in Section 6 of the paper, designed as "stand-alone scripts" to ease the reproduction of the simulation results. 

  - HDsimu_model1.r - script for reproducing simulation results presented in Figure 2 for Model 1
  - HDsimu_model2.r - script for reproducing simulation results presented in Figure 2 for Model 2 

Authors

Radojicic, U. and Virta, J.

License

GNU GPLv3

References

Luo, W. and Li, B. (2021): On order determination by predictor augmentation. Biometrika, 108(3), pp.557-574.

Radojicic, U. and Virta J. (2025): Dimension estimation in PCA model using high-dimensional data augmentation. arXiv preprint arXiv:2502.04220.  

Schott, J. R. (2006): A high-dimensional test for the equality of the smallest eigenvalues of a covariance matrix. 
Journal of Multivariate Analysis 97, 827–843.
