##log-normal study for misspecification paper

install.packages("LaplacesDemon")
install.packages("compositions")
install.packages("FNN")
install.packages("e1071")
install.packages("nortest")

library(LaplacesDemon) ##to use is.positive.definite() function
library(compositions) ##to generate log normal data
library(FNN) ##mutinfo function
library(e1071) ##for skewness
library(nortest) ##for ad test


start_time<- Sys.time()
print(start_time)


##generate random mu
mu_1<-rnorm(750, mean = 1.5, sd = 0.5)
mu_2<-rnorm(750, mean = 4, sd = 2)
mu_3<-rnorm(750, mean = 7, sd = 0.75)
mu_4<-rnorm(750, mean = 10, sd = 0.25)


##generate random sigmas
n <- 100  # Number of matrices
p <- 4  # Dimension
df <- 10  # Degrees of freedom
Sigma <- toeplitz((p:1)/p)  # the matrix parameter of the distribution
# Draw n Wishart distributed matrices
Sigmas <- drop(rWishart(1000, df, Sigma))

for(i in 1:750) {
  mu<- c(mu_1[i], mu_2[i], mu_3[i], mu_4[i])
  mu<- abs(mu)
  Sigma<- Sigmas[, ,i:i]
  Sigma<- abs(Sigma)
  is.positive.definite(Sigma)
  
  #TO CONVERT PARAMETER
  convert_mu_1<- log( (mu)^2/ sqrt( (mu)^2+ diag(Sigma) ))
  
  ##convert parameter
  var11<- log(1+ ( Sigma[1,1] / mu[1]^2) )
  var22<- log(1+ ( Sigma[2,2] / mu[2]^2) )
  var33<- log(1+ ( Sigma[3,3] / mu[3]^2) )
  var44<- log(1+ ( Sigma[4,4] / mu[4]^2) )
  var12<- log(1+ ( Sigma[1,2]/  (mu[1]*mu[2])) )
  var13<- log(1+ ( Sigma[1,3] / (mu[1]*mu[3])) )
  var14<- log(1+ ( Sigma[1,4] / (mu[1]*mu[4])) )
  var23<- log(1+ ( Sigma[2,3] / (mu[2]*mu[3])) )
  var24<- log(1+ ( Sigma[2,4] / (mu[2]*mu[4])) )
  var34<- log(1+ ( Sigma[3,4]/  (mu[3]*mu[4])) )
  
  cov_convert_sigma_1<- matrix( c( var11, var12, var13, var14,
                                   var12, var22, var23, var24,
                                   var13, var23, var33, var34,
                                   var14, var24, var34, var44), nrow = 4, ncol = 4)
  
  
  is.positive.definite(cov_convert_sigma_1)
  
  
  if (is.positive.definite(cov_convert_sigma_1) == FALSE)
  {
    print(i)
  }
  
  else {
    
    a<- i
    #ICA_lognormal<- list()
    #ICA_normal<- list()
    
    
    data_log<-  rlnorm.rplus(2000,convert_mu_1,cov_convert_sigma_1)
    
    T0<-data_log[,1]
    T1<-data_log[,2]
    S0<-data_log[,3]
    S1<-data_log[,4]
    
    delta_T<- T1-T0
    delta_S<- S1-S0
    
    rho<- cor(delta_T, delta_S)
    ICA_normal_2<- rho^2
    
    #mutual_info<-   mutinfo(delta_T, delta_T, k=10, direct=TRUE) it is slow
    mutual_info_2<- mutinfo(delta_T, delta_S, k=10, direct=FALSE)
    mutual_info_2<- mean(mutual_info_2)
    
    #ICA_lognormal_1<- 1-exp(-2*mutual_info)
    ICA_lognormal_2<- 1-exp(-2*mutual_info_2)
    
    diff<- ICA_lognormal_2- ICA_normal_2
    
    
    data_log_2<- cbind(a,T0,T1,S0,S1)
    outputs<- cbind(a, rho, ICA_normal_2, ICA_lognormal_2, diff)
    
    write.table(outputs, file="log_normal_output.txt",  sep = "\t", row.names = FALSE, col.names = FALSE, append=TRUE)
    write.table(data_log_2, file="data_lognormal.txt", sep = "\t", row.names = FALSE, col.names = FALSE, append=TRUE)
    
    
    
  } }

end_time<- Sys.time()
print(end_time)

