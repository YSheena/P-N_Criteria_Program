functions{
 // "xi": the funciotn xi's in the exponetial family model
 // Input
 // x1: the vector of the coninuous variabel of "p1"-dimensions
 // x2: the discrete variable of one-dimension, integers from 1:"Rq""
 // Output 
 // all the possible cross between the elements of x=(x1,x2)
 // Note that you can define any xi's as you like 
 vector xi (vector x1, int x2){
  int count=0;
  int n1 = num_elements(x1);
  int d = choose(n1,2)+n1; 
  vector[d] x;
  for (i in 1:(n1-1)){
    for (j in (i+1):n1){
      count += 1;
      x[count]= x1[i]*x1[j];
    }
  }
  for (i in 1:n1){
    count += 1;
    x[count]= x2*x1[i];
  }
  return x; 
 }
 // "base_pdf": The log density of the base measure mu.
 // Input: x=(x1,x2), beta_1, beta_2, Rqn
 // Output: the value of the log d(mu)/dx;
 // Independent Beta distributions for x1, whose parameters are beta_1, beta_2.
 // Multinomial distribution for x2 defined by Rqn(log probabiltiy math func.), 
 // which is independet of Beta ditributions.
 // Note that you can define any xi's as you like 
 real base_pdf(vector x1, int x2, vector beta_1, vector beta_2, vector Rqn){
  return (beta_1-1)'*log(x1)+(beta_2-1)'*log(1-x1)+Rqn[x2];
 }
 // "log_exp_pdf": the log p.d.f of x
 // Input: x=(x1,x2), beta_1, beta_2, theta, Rqn
 // Output: log(g(x;theta))+log(d(mu)/dx)
 real log_exp_pdf(vector x1, int x2, vector theta,
 vector beta_1, vector beta_2, vector Rqn){
  return base_pdf(x1, x2, beta_1,beta_2, Rqn) + (theta)'*xi(x1,x2);
 }
 // "log_exp_con_lpdf": the log p.d.f. of x1
 // Input: x=(x1,x2), beta_1, beta_2, theta, Rqn
 // Output: log(p.d.f. of x1)
 real log_exp_con_lpdf(vector x1, vector theta, 
 vector beta_1, vector beta_2, vector Rqn){
  int n = num_elements(Rqn);
  vector[n] x;
  for (i in 1:n){
   x[i] = log_exp_pdf(x1,i,theta, beta_1, beta_2, Rqn) ;
  }
  return log_sum_exp(x);
 }
 // "prpb_x2_given_x1": probability mass function of x2 given x1
 // Input: x=(x1,x2), beta_1, beta_2, theta, Rqn
 // Output: p.m.f of x2 given x1
 vector prpb_x2_given_x1(vector x1, vector theta,
 vector beta_1, vector beta_2, vector Rqn){
  int n = num_elements(Rqn);
  vector[n] y;
  for (i in 1:n){
   y[i] = log_exp_pdf(x1,i,theta, beta_1, beta_2, Rqn);
  }
  return softmax(y);
 }
}
data{
  int p1;
  int n_base;
  int Rq;
  vector[n_base] theta ;
  vector[p1] beta_1;
  vector[p1] beta_2;
  vector[Rq] Rqn;
}
parameters{
  vector<lower=0,upper=1>[p1] x1; //Step.1-1
}
model{
  target += log_exp_con_lpdf(x1 | theta, beta_1, beta_2, Rqn);
}
// Generation of xi's, each of which is the "theoritically" avaraged over
// the coditional distribution of x2 given x1
generated quantities{
  vector[n_base] xs = rep_vector(0,n_base);
  matrix[n_base,n_base] xsq = rep_matrix(xs,n_base);
  for (s in 1:Rq) {
    vector[n_base] x = xi(x1,s); 
    xs = xs + prpb_x2_given_x1(x1,theta,beta_1,beta_2,Rqn)[s]*x ;// 1-2,1-3
    xsq = xsq + prpb_x2_given_x1(x1,theta,beta_1,beta_2,Rqn)[s]*x*x';
  }
}
