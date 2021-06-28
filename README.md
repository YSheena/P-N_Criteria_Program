# P-N_Criteria_Program
The programs used in "The convergence speed of MLE to the information projection of an exponential family --a criteria for the model dimension and the sample size-- "
# Five main programs
1. expo_est_Db=1_wine_beta_ver4.R: The risk cacluation for the full model (p=66)
2. expo_est_Db=1_wine_beta_ver4_sub.R: The risk cacluation for the submodel (xi's with high-correlaitions are omitted)
3. expo_est_Db=1_wine_beta_ver4_classifying.R: The crossvalidation of the classifyer using the full model.
4. expo_est_Db=1_wine_beta_ver4_classifying_submodel.R: The crossvalidation of the classifyer using the submodel (xi's with high-correlaitions are omitted).
5. expo_est_Db=1_wine_beta_ver4_2_sub.R: The risk cacluation for the submodel (xi's with high-correlaitions are omitted) using different datasets for model formulation and parameter estimation.
# Two Stan programs: All the main programs need these programs.
1. expo_stan_Db=1_ver2_beta.stan
2. expo_stan_Db=1_ver3_beta.stan
# One subroutine program
1. grid_making.R: used for the main programs 1 & 2.
2. grid_making_ver2.R : used for the main program 5.
