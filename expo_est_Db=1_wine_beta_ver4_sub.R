# This program is aimed to calculate the estimation risk 
# for an exponential family distribution specialized for wine data.
#
# This program is for the submodel with reduced xi's 
# (the base distribution remains the same)


# The main variables of the model
# x1:"p1"-dimensional continuous variables
# x2:1-dimensional discrete variables from 1 to "Rq"
# Rq:the number of possible outcomes of x2
# Rqn:the Rq-dimensional vector of log probability mass function of x2
# n_base: the number of the xi's in the full model
# n_base_s: the number of the xi's in the submodel
# i.e. "p", the dimension of the exponential family distribution

### Model Description ###
#
# 1. the base p.d.f, i.e. dmu/dx ;
# x1_i(i=1,...,p1), x2 are all independent and
# x1_i~beta(a_i,b_i), x2 ~ multinomial dist. 
# with the parameters given empirically from the sample.
# log of the base p.d.f.
# (a_1-1)log(x1_1)+(b_1-1)log(1-x1_1)
#+(a_2-1)log(x1_2)+(b_2-1)log(1-x1_2)+...
#+log(I(x2=1)*p(x2=1)+...+I(x2=Rq)*p(x2=Rq))
#
# 2. xi's
# all the possible cross between the elements of x=(x1,x2)
# In this submodel program, one will removed from the pair of xi's 
# with high correlation. It will be carried out later.

xi <- function(x){
  count <- 0;
  n <- length(x)
  y <- rep(0,(n-1)*n/2)
  for (i in 1:(n-2)){
    for (j in (i+1):(n-1)){
      count <- count + 1
      y[count] <- x[i]*x[j]
    }
  }
  for (i in 1:(n-1)){
    count <- count + 1
    y[count] <- x[n]*x[i]
  }
  return(y)
}
### End of Model Description ###

.libPaths("./library")
.libPaths()
install.packages("rstan",dependencies = TRUE)
install.packages("nleqslv")
install.packages("gtools")
install.packages("tidyverse")
install.packages("tensorA")
install.packages("rlang")
library("rlang")
library("rstan")
library("MASS") #We use "ginv" function
library("nleqslv") #Nonlinear Equation Solver
#library("tidyr") #For data polishing
library("parallel") # Parallel calculation.
library("tensorA") # tensor, Einstein summation
rstan_options(auto_write = TRUE)
options(mc.cores=parallel::detectCores()) 
source("grid_making.R")

winedata <- read.csv("winequality-red.csv",header=TRUE)

# Scaling the continuous variables x1 so that
# they are distributed within [0,1]
max_x1 <- apply(winedata[,-12],2,max)
winedata_con <- as.matrix(winedata[,-12])%*%diag(1/max_x1)

# Change "quality" so that the range 1:Rq
winedata_dis <- winedata[12]-2

# Normalize data set as data.frame 
winedata_norm <- as.data.frame(cbind(winedata_con,winedata_dis))
names(winedata_norm) <- names(winedata)


p1 <- 11; # the dimension x1 "chemical substances"
Rq <- 6; # "quality" min 3, max 8
n_base <- p1*(p1+1)/2 # the number of xi's for the full model
N <- nrow(winedata_norm) # The number of the individuals of the data

# Base distribution for x2
# Let prob_i be the estimated probability P(x2=i) 
# prob = c(p_1,...,p_Rq)
prob <- table(winedata_norm$quality)
minp <- min(log(prob))

# the log of probability mass function of x2 up to a constant
Rqn <- log(prob)-rep(minp,Rq)

# log(p.d.f) for the base distribution of x1
# log of p.d.f. of beta(a,b) : (a-1)log(x)+(b-1)log(1-x)
# mean m=a/(a+b), variance v=(ab)/((a+b)^2*(a+b+1))
# ==> a = {m^2(1-m)\pm m(sqrt(m^2(1-m)^2-4v^2)}/(2v)
#     b = a(1-m)/m

# the mean and variance of x1 
m <- apply(winedata_norm[,-12],2,mean)
v <- apply(winedata_norm[,-12],2,var)

#  empirically given parameters of the base distribution of x1
beta_1 <- m^2/v*(1-m)-m
beta_2 <- beta_1*(1-m)/m

# Empirical Data of xi's
xi_emp <- as.data.frame(t(apply(winedata_norm,1,xi)))

#Find high correlated xi's
Corr_s <- cor(xi_emp)
Corr_s[lower.tri(Corr_s,diag = TRUE)] <- 0
# the index of the xi's that will be removed
# one of the pair whose correlation is over 0.98
zsi <- (which(abs(Corr_s) > 0.95, arr.ind=T))[,"col"]
zsi <- unique(zsi)
n_zsi <- length(zsi)
# the number of xi's in the submodel
n_base_s <- n_base - n_zsi


# Sample mean of xi's in the submodel
eta0 <- apply(xi_emp,2,mean)
eta0_s <- eta0[-zsi]
# Sample covariance of xi's
Sigma <- cov(xi_emp)
Sigma_s <- Sigma[-zsi,-zsi]

# Function "rv_model"
# Generataor of rv's from g(x;theta)
# Input: "theta", the parameter of g(x;theta)
#        "seedchange", the seed for function "stan"
# Output: all random objects in the output of "stan"
rv_model <- function(theta,seedchange) {
  data <- list(p1=p1,n_base=n_base,Rq=Rq,Rqn=Rqn,theta=theta,
               beta_1=beta_1,beta_2=beta_2)
  fit <- stan(file = "expo_stan_Db=1_ver3_beta.stan", data=data, 
              control = list(adapt_delta = 0.99,max_treedepth = 15),
              chains=4, seed= seedchange, iter=10000)
  return(rstan::extract(fit) )
}

# Function "eta_DdPsi"
# Input : "theta", the parameter of g(x;theta)
#        "seedchange", the seed for function "stan"
# Ouput : List of two objects
# 1. "eta":simulated mean of xi's under g(x;theta)
# 2. "DdPsi":simulated covariance matrix of xi's under g(x;theta)
# Note that the "eta", "DdPsi" is calculated from the "theoretical" expectation
# over the distribution of x2, hence more accurate than those calculated 
# "empirically" from the samples from "rv_model"
eta_DdPsi <- function(theta,seedchange) {
  data <- list(p1=p1,n_base=n_base,Rq=Rq,Rqn=Rqn,theta=theta,
               beta_1=beta_1,beta_2=beta_2)
  fit <- stan(file = "expo_stan_Db=1_ver2_beta.stan", data=data, 
              control = list(adapt_delta = 0.99,max_treedepth = 15),
              chains=4, seed= seedchange, iter=10000)
  rx <- rstan::extract(fit) 
  eta <- apply(rx$xs,2,mean) 
  Ddpsi <- apply(rx$xsq,c(2,3),mean) - eta%*%t(eta)
  return(list(eta = eta, Ddpsi = Ddpsi))
}

# grids for the calculation of moments or cumulants 
# This takes several minutes.
grids <- grid_making(n_base)
M2 <- grids$M2 ; M3 <- grids$M3 ; M4 <- grids$M4

# gridvec3s, gridvec4s is the new grids for the submodel
index_s <- (1:n_base)[-zsi]
gridvec3s <- expand.grid(index_s,index_s,index_s)
gridvec4s <- expand.grid(index_s,index_s,index_s,index_s)

# The sample moment calculation for the chosen variables
prodfunc<- function(data,x){
  mean(apply(data[,x],1,prod))
}

# New function "theta_f" 
# Input: theta_s, the theta for the submodel
# Output: the theta for the full model, i.e. the zero-interporated vector 
theta_f <- function(theta_s,zsi){
  for (i in 1:n_zsi){
    if (zsi[i] < n_base){
      theta_s <- c(theta_s[1:(zsi[i]-1)],0,theta_s[zsi[i]:length(theta_s)])
    }else{
      theta_s <- c(theta_s[1:(zsi[i]-1)],0)
    }
  }
  return(theta_s)
}

##### Repeat the risk calucualtion with different seeds in stan #####
result_m <- c()
for (seedchange in 1234:1244) {

  ## The solution of theta_* using "nleqslv"
  # Procedure 
  # Step 1: Under the given theta(t), calculate
  # eta(t), the mean of xi's, and DdPsi(t), the covariance matrix of xi's
  # from the samples generated by MCMC.
  # "eta_DdPsi" function with "expo_stan_Db=12_beta.stan" is used.
  # Step 2: By Newton-Raphson method (together with another searching method), 
  # search the value of theta_* iteratively.
  # theta(t+1) = theta(t)-(eta(t)-eta)%*%DdPsi(t)^(-1)
  
  theta_z1 <- rep(0,n_base_s) # initial value of theta
  result <- nleqslv(theta_z1,fn = function(x) {
    eta_DdPsi(theta_f(x,zsi),seedchange)$eta[-zsi] - eta0_s},
    jac = function(x) {eta_DdPsi(theta_f(x,zsi),seedchange)$Ddpsi[-zsi,-zsi]} ,
    method = "Newton",
    global = "cline",
    xscalm = "auto",jacobian=TRUE,control=list(trace=1))
  
  theta_star <- result$x # the solution of theta_*
  saveRDS(result,"result.rds")
  
  # N^{-1} order term coeff. of the estimation risk
  Psin <- ginv(result$jac) #just for the theta's in the submodel !
  term0 <- sum(diag(Psin%*%Sigma_s))
  
  # Sample from g(x;theta_*)
  x_gstar <- rv_model(theta_f(theta_star,zsi),seedchange)
  x_model <- cbind(x_gstar$x1,x_gstar$x2)
  
  # Data of xi's from the model
  xi_gstar <- as.data.frame(t(apply(x_model,1,xi)))
  # Sample mean of xi's
  eta1 <- apply(xi_gstar,2,mean)
  # Sample covariance of xi's
  DdPsi <- cov(xi_gstar)
 
  
  # N^{-1} order term coeff. calculated in another way
  term0a <- sum(diag(ginv(DdPsi[-zsi,-zsi])%*%Sigma_s))
  
  ## First-order term of ED ##
  FirN <- term0/(2*N)
  
  
  ### N^{-2}-order term calculation ###
  
  # Higher order moments from the empirical data and the model g(x;eta_*)
  detectCores() 
  cl1 <- makeCluster(detectCores())
  clusterExport(cl1, varlist = 
                  list('xi_gstar','xi_emp','prodfunc'))
  moment2e <- parSapply(cl1,M2,function(x) {prodfunc(xi_emp,x)})
  moment3e <- parSapply(cl1,M3,function(x) {prodfunc(xi_emp,x)})
  moment2m <- parSapply(cl1,M2,function(x) {prodfunc(xi_gstar,x)})
  moment3m <- parSapply(cl1,M3,function(x) {prodfunc(xi_gstar,x)})
  #the most time-consuming, about 1 hr
  moment4m <- parSapply(cl1,M4,function(x) {prodfunc(xi_gstar,x)}) 
  stopCluster(cl1)
  
  saveRDS(moment4m,"moment4m.rds")
  
  # Higher order cumulants from the empirical data and the model g(x;eta_*)
  cl2 <- makeCluster(24)
  clusterExport(cl2, varlist=
                  list('a_v','M3','M4','eta0','eta1', 'moment2e', 'moment2m',
                       'moment3e','moment3m','moment4m'))
  
  
  cumu3e <- parSapply(cl2, M3, 
                      function(x){
                        x <- unlist(x)
                        moment3e[a_v(x)]-
                          eta0[x[1]]*moment2e[a_v(x[c(2,3)])]-
                          eta0[x[2]]*moment2e[a_v(x[c(1,3)])]-
                          eta0[x[3]]*moment2e[a_v(x[c(1,2)])]+
                          2*eta0[x[1]]*eta0[x[2]]*eta0[x[3]]
                      })
  
  cumu3m <- parSapply(cl2, M3, 
                      function(x){
                        x <- unlist(x)
                        moment3m[a_v(x)]-
                          eta1[x[1]]*moment2m[a_v(x[c(2,3)])]-
                          eta1[x[2]]*moment2m[a_v(x[c(1,3)])]-
                          eta1[x[3]]*moment2m[a_v(x[c(1,2)])]+
                          2*eta1[x[1]]*eta1[x[2]]*eta1[x[3]]
                      })
  
  cumu4m <- parSapply(cl2, M4, 
                      function(x){
                        x <- unlist(x)
                        moment4m[a_v(x)]-
                          moment3m[a_v(x[c(1,2,3)])]*eta1[x[4]]-
                          moment3m[a_v(x[c(1,2,4)])]*eta1[x[3]]-
                          moment3m[a_v(x[c(1,3,4)])]*eta1[x[2]]-
                          moment3m[a_v(x[c(2,3,4)])]*eta1[x[1]]-
                          moment2m[a_v(x[c(1,2)])]*moment2m[a_v(x[c(3,4)])]-
                          moment2m[a_v(x[c(1,3)])]*moment2m[a_v(x[c(2,4)])]-
                          moment2m[a_v(x[c(1,4)])]*moment2m[a_v(x[c(2,3)])]+
                          2*moment2m[a_v(x[c(1,2)])]*eta1[x[3]]*eta1[x[4]]+
                          2*moment2m[a_v(x[c(1,3)])]*eta1[x[2]]*eta1[x[4]]+
                          2*moment2m[a_v(x[c(1,4)])]*eta1[x[2]]*eta1[x[3]]+
                          2*moment2m[a_v(x[c(2,3)])]*eta1[x[1]]*eta1[x[4]]+
                          2*moment2m[a_v(x[c(2,4)])]*eta1[x[1]]*eta1[x[3]]+
                          2*moment2m[a_v(x[c(3,4)])]*eta1[x[1]]*eta1[x[2]]-
                          6*eta1[x[1]]*eta1[x[2]]*eta1[x[3]]*eta1[x[4]]
                      })
  
  stopCluster(cl2)
  
  
  # Change them to cover all the possible indexes,  1 <= i, j, k, (l)<= n_base_s 
  cl3 <- makeCluster(24)
  
  clusterExport(cl3,varlist = list('cumu3e','cumu3m','cumu4m','a_v'))
  
  # these cumulants are reduced for the subodel
  cumu3ef <- parApply(cl3,gridvec3s,1,FUN = function (x) {
    cumu3e[a_v(sort(x))]})
  
  cumu3mf <- parApply(cl3,gridvec3s,1,FUN = function (x) {
    cumu3m[a_v(sort(x))]})
  
  # takes several minutes
  cumu4mf <- parApply(cl3,gridvec4s,1,FUN = function (x) {
    cumu4m[a_v(sort(x))]})
  
  stopCluster(cl3)
 
  # First term of the N^{-2}-order term
  
  # 1. Change all objects into tensors
  
  Psin_t1 <- as.tensor(Psin, dims=c(I=n_base_s,J=n_base_s))
  Psin_t2 <- as.tensor(Psin, dims=c(K=n_base_s,L=n_base_s))
  Psin_t3 <- as.tensor(Psin, dims=c(M=n_base_s,S=n_base_s))
  cumu3ef_t <- to.tensor(cumu3ef,dims=c(I=n_base_s,K=n_base_s,M=n_base_s))
  cumu3mf_t <- to.tensor(cumu3mf,dims=c(J=n_base_s,L=n_base_s,S=n_base_s))
  
  
  # 2. Multiplication between tensors

  # It is important to avoid outer product ! 
  # Otherwise the size of the tensors easilly get over the memory.
  
  term1 <- cumu3mf_t %e% Psin_t3 %e% Psin_t2 %e% Psin_t1 %e% cumu3ef_t
  
  
  # Second term of the N^{-2}-order term

  Psin_t1 <- as.tensor(Psin, dims=c(I=n_base_s,J=n_base_s))
  Psin_t2 <- as.tensor(Psin, dims=c(K=n_base_s,L=n_base_s))
  Psin_t3 <- as.tensor(Psin, dims=c(M=n_base_s,S=n_base_s))
  Psin_t4 <- as.tensor(Psin, dims=c(O=n_base_s,P=n_base_s))
  Psin_t5 <- as.tensor(Psin, dims=c(U=n_base_s,V=n_base_s))
  Sigma_t1 <- as.tensor(Sigma_s, dims=c(J=n_base_s,L=n_base_s))
  Sigma_t2 <- as.tensor(Sigma_s, dims=c(S=n_base_s,P=n_base_s))
  Sigma_t3 <- as.tensor(Sigma_s, dims=c(J=n_base_s,S=n_base_s))
  Sigma_t4 <- as.tensor(Sigma_s, dims=c(L=n_base_s,P=n_base_s))
  Sigma_t5 <- as.tensor(Sigma_s, dims=c(J=n_base_s,P=n_base_s))
  Sigma_t6 <- as.tensor(Sigma_s, dims=c(L=n_base_s,S=n_base_s))
  cumu3mf_t1 <- to.tensor(cumu3mf,dims=c(I=n_base_s,K=n_base_s,U=n_base_s))
  cumu3mf_t2 <- to.tensor(cumu3mf,dims=c(M=n_base_s,O=n_base_s,V=n_base_s))
  
  term21 <- Sigma_t2 %e% Psin_t4 %e% Psin_t3 %e% cumu3mf_t2 %e%
    Sigma_t1 %e% Psin_t5 %e% Psin_t2 %e% Psin_t1 %e% cumu3mf_t1
  
  term22 <- Sigma_t4 %e% Psin_t4 %e% Psin_t3 %e% cumu3mf_t2 %e%
    Sigma_t3 %e% Psin_t5 %e% Psin_t2 %e% Psin_t1 %e% cumu3mf_t1
  
  term23 <- Sigma_t6 %e% Psin_t4 %e% Psin_t3 %e% cumu3mf_t2 %e%
    Sigma_t5 %e% Psin_t5 %e% Psin_t2 %e% Psin_t1 %e% cumu3mf_t1
  
  term2 <- term21+term22+term23
  
  # Third term of the N^{-2}-order term
  
  Psin_t1 <- as.tensor(Psin, dims=c(I=n_base_s,J=n_base_s))
  Psin_t2 <- as.tensor(Psin, dims=c(K=n_base_s,L=n_base_s))
  Psin_t3 <- as.tensor(Psin, dims=c(M=n_base_s,S=n_base_s))
  Psin_t4 <- as.tensor(Psin, dims=c(O=n_base_s,P=n_base_s))
  Sigma_t1 <- as.tensor(Sigma_s, dims=c(J=n_base_s,L=n_base_s))
  Sigma_t2 <- as.tensor(Sigma_s, dims=c(S=n_base_s,P=n_base_s))
  Sigma_t3 <- as.tensor(Sigma_s, dims=c(J=n_base_s,S=n_base_s))
  Sigma_t4 <- as.tensor(Sigma_s, dims=c(L=n_base_s,P=n_base_s))
  Sigma_t5 <- as.tensor(Sigma_s, dims=c(J=n_base_s,P=n_base_s))
  Sigma_t6 <- as.tensor(Sigma_s, dims=c(L=n_base_s,S=n_base_s))
  cumu4mf_t <- to.tensor(cumu4mf, dims=c(I=n_base_s,K=n_base_s,M=n_base_s,O=n_base_s))
  
  term31 <- cumu4mf_t %e% Psin_t1 %e% Psin_t2 %e%
    Psin_t3 %e% Psin_t4 %e% Sigma_t1  %e% Sigma_t2
  
  term32 <- cumu4mf_t %e% Psin_t1 %e% Psin_t2 %e% 
    Psin_t3 %e% Psin_t4 %e% Sigma_t3  %e% Sigma_t4
  
  term33 <- cumu4mf_t %e% Psin_t1 %e% Psin_t2 %e% 
    Psin_t3 %e% Psin_t4 %e% Sigma_t5  %e% Sigma_t6
  
  term3 <- term31+term32+term33
  
  # The N^{-2}-order term
  SecN <- (-8*term1 + 9*term2 -3*term3)/(24*N^2)
  
  # The estimation risk in total
  EstR <- FirN + SecN
  
  cat("FirN=",FirN,"\n SecN=",SecN, "\n EstR=",EstR, "\n")
  
  # Saving the results in "result_m"
  result_m <- rbind(result_m, 
                    c(theta_star,eta1[-zsi],term0,term1,term2,term3,FirN,SecN,EstR))
} 
##### End of the repetition #####
saveRDS(result_m,"result_m.rds")

# Summary for the repetition 
risk_mean <- apply(result_m,2,mean)
risk_sd <- apply(result_m,2,sd)
result_m <- rbind(rbind(result_m,risk_mean),risk_sd)

# the colum names of the data set of the repetition result
theta_star_name <- c()
eta1_name <- c()
for (i in 1:n_base_s) {
  theta_star_name <- append(theta_star_name,sprintf("theta_*_%2d",i))
  eta1_name <- append(eta1_name, sprintf("eta1_%2d",i))
  i <- i + 1
}

# Making of the data.frame for the repetition result
colnames(result_m) <- c(theta_star_name, eta1_name, 
                        "term0","term1","term2","term3","FirN","SecN","EstR")
result_d <- as.data.frame(result_m)

write.csv(result_d, "result_d.csv")
