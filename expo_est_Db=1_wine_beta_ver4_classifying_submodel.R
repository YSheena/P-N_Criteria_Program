# This program is aimed to carry out the cross validation of 
# the classifier made from the exponential family model g(x;theta)
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

# Function "log_exp_pdf"
# Input: "x1": the continuous variable
#        "x2": the discrete variable
#        "Rqn": log p.m.f.of x2 up to a constant
#        "theta": the parameter vector of the exponential distribution model
# Output: the non-constant part w.r.t. x2 of the log p.d.f. of x=(x1,x2)
log_exp_pdf <- function(x1, x2, theta, Rqn){
  return (Rqn[x2] + theta%*%(xi(c(x1,x2))))
}

# Function "prpb_x2_given_x1"
# Input: "x1": the continuous variable
#        "Rqn":log p.m.f.of x2
#        "theta": the parameter vector of the exponential distribution model
# Output: the vector of  the log p.d.f. for each value of x2 given x1 
# up to a constant
prpb_x2_given_x1 <- function(x1, theta, Rqn){
  n = length(Rqn)
  y <- c()
  for (i in 1:n){
    y[i] = log_exp_pdf(x1,i,theta, Rqn);
  }
  return(y)
}


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


# Function "theta_f" 
# Input: theta_s, the theta for the submodel
# Output: the theta for the full model, i.e. the zero-interporated vector 
theta_f <- function(theta_s,zsi){
  for (i in 1:length(zsi)){
    if (zsi[i]< n_base){
      theta_s <- c(theta_s[1:(zsi[i]-1)],0,theta_s[zsi[i]:length(theta_s)])
    }else{
      theta_s <- c(theta_s[1:(zsi[i]-1)],0)
    }
  }
  return(theta_s)
}

# Bayes discriminant function
bayes_dis <- function(x1, theta, Rqn){
  pvec <- prpb_x2_given_x1(x1, theta, Rqn)
  return(which(pvec == max(pvec))+2)
}

.libPaths("./library")
.libPaths()
install.packages("rstan",dependencies = TRUE)
install.packages("nleqslv")
#install.packages("gtools")
#install.packages("tidyverse")
#install.packages("tensorA")
install.packages("rlang")
library("rlang")
library("rstan")
library("MASS") #We use "ginv" function
library("nleqslv") #Nonlinear Equation Solver
#library("tidyr") #For data polishing
#library("parallel") # Parallel calculation.
#library("tensorA") # tensor, Einstein summation
rstan_options(auto_write = TRUE)
options(mc.cores=parallel::detectCores()) 


winedata <- read.csv("winequality-red.csv",header=TRUE)

# the record of the classification result for the training data
all_result_train<- data.frame(class_pre_train=factor(c(0)),
                              class_true_train=factor(c(0)))
# the record of the classification result for the test data
all_result <- data.frame(class_pre=factor(c(0)),
                         class_true=factor(c(0)))
# the record of the estimation risk (firs-order-term only)
result_FirN <- c()

##################### Loop for the cross validation ############################
for (i in 1:10) {
# Randomly divide the data into the train data and the test data
dt = sort(sample(nrow(winedata), nrow(winedata)*.9))
train <- winedata[dt,]
test <- winedata[-dt,]

# Scaling the continuous variables x1 so that
# they are distributed within [0,1]
max_x1 <- apply(train[,-12],2,max)
winedata_con <- as.matrix(train[,-12])%*%diag(1/max_x1)

# Change "quality" so that the range 1:Rq
winedata_dis <- train[12]-2

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

#Find the 10 highest correlated xi's
Corr_s <- cor(xi_emp)
Corr_s[lower.tri(Corr_s,diag=TRUE)] <- 0
# the index of the xi's that will be removed
#mar_corr <- sort(abs(Corr_s),decreasing = TRUE)[20]
#zsi <- unique(which(abs(Corr_s) >= mar_corr, arr.ind=T)[,"col"])
zsi <- unique(which(abs(Corr_s) >= 0.95, arr.ind=T)[,"col"])
n_zsi <- length(zsi)
# the number of xi's in the submodel
n_base_s <- n_base - n_zsi


# Sample mean of xi's in the submodel
eta0 <- apply(xi_emp,2,mean)[-zsi]
# Sample covariance of xi's
Sigma <- cov(xi_emp)[-zsi,-zsi]


## The solution of theta_* using "nleqslv"
# Procedure 
# Step 1: Under the given theta(t), calculate
# eta(t), the mean of xi's, and DdPsi(t), the covariance matrix of xi's
# from the samples generated by MCMC.
# "eta_DdPsi" function with "expo_stan_Db=12_beta.stan" is used.
# Step 2: By Newton-Raphson method (together with another searching method), 
# search the value of theta_* iteratively.
# theta(t+1) = theta(t)-(eta(t)-eta)%*%DdPsi(t)^(-1)
seedchange <- 1234  
theta_z1 <- rep(0,n_base_s) # initial value of theta
result <- nleqslv(theta_z1,fn = function(x) {
  eta_DdPsi(theta_f(x,zsi),seedchange)$eta[-zsi] - eta0},
  jac = function(x) {eta_DdPsi(theta_f(x,zsi),seedchange)$Ddpsi[-zsi,-zsi]} ,
  method = "Newton",
  global = "cline",
  xscalm = "auto",jacobian=TRUE,control=list(trace=1))

theta_f(rep(1,n_base_s),zsi)

theta_star <- result$x # the solution of theta_*
saveRDS(result,"result.rds")

# N^{-1} order term coeff. of the estimation risk
Psin <- ginv(result$jac)
term0 <- sum(diag(Psin%*%Sigma))

## First-order term of ED ##
FirN <- term0/(2*N)
result_FirN <- c(result_FirN,FirN)

#first-order p-n criteria 
alpha <- 0.05
if (FirN < 8*alpha^2) {
  print("the model passed p-n criteria")
}else{
  print("sample size is too small")
}

# Classify each individual in the train data
class_pre_train <- as.vector(
  apply(winedata_con,1,
        function(x1) {
          bayes_dis(x1, theta=theta_f(theta_star,zsi), Rqn=Rqn) 
        }
  )
)

result_t_train <- table(class_pre_train, class_true_train = train[,(p1+1)]) 


# Classify each individual in the test data
max_x1_t <- apply(test[,-(p1+1)],2,max)
test_x1 <- as.matrix(test[,-(p1+1)])%*%diag(1/max_x1_t)

class_pre <- as.vector(
  apply(test_x1,1,
      function(x1) {
        bayes_dis(x1, theta=theta_f(theta_star,zsi), Rqn=Rqn) 
      }
  )
)

result_t <- table(class_pre,class_true = test[,p1+1])

all_result_train <- dplyr::full_join(all_result_train,as.data.frame(result_t_train),
                                     by=c("class_pre_train","class_true_train"))
all_result <- dplyr::full_join(all_result,as.data.frame(result_t),
                               by=c("class_pre","class_true"))
}


all_result <- all_result[-1,]
write.csv(all_result,"ver4_sub2_all_result_corr_95")
count_result <- apply(all_result[c(-1,-2)],1,function(x) sum(x,na.rm=TRUE))
all_result <- dplyr::mutate(all_result, count = count_result)
all_result_count <- all_result[c(1,2,ncol(all_result))]
cross_table <- tidyr::spread(all_result_count, key="class_true",value="count")
write.csv(cross_table,"cross_table.csv")
(cross_table[1,4]+cross_table[2,5]+cross_table[3,6])/
  sum(cross_table[,2:7],na.rm=TRUE)

all_result_train <- all_result_train[-1,]
write.csv(all_result_train,"ver4_sub2_all_result_train_corr_95")
count_result_train <- apply(all_result_train[c(-1,-2)],1,function(x) sum(x,na.rm=TRUE))
all_result_train <- dplyr::mutate(all_result_train, count = count_result_train)
all_result_count_train <- all_result_train[c(1,2,ncol(all_result_train))]
cross_table_train <- tidyr::spread(all_result_count_train, key="class_true_train",value="count")
write.csv(cross_table_train,"cross_table_train.csv")
(cross_table_train[1,4]+cross_table_train[2,5]+cross_table_train[3,6])/
  sum(cross_table[,2:7],na.rm=TRUE)
                            