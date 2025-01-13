##log-normal study for misspecification paper

# install.packages("LaplacesDemon")
# install.packages("compositions")
# install.packages("FNN")
# install.packages("nortest")

library(LaplacesDemon) ##to use is.positive.definite() function
library(compositions) ##to generate log normal data
library(FNN) ##mutinfo function
library(cubature) ##for integral (adaptIntegrate)
library(MASS) # For rWishart function


##library(utils) ##write csv and write table
#library(pracma) ##for integral2
#library(emdbook) ##for 4 dim normal


start_time<- Sys.time()
print(start_time)

generate_parameters <- function(p = 4) {
  # Generate random mu
  mu <- rbind(
    rnorm(1, mean = runif(1, 1, 2), sd = runif(1, 6, 7)),  # mu_1
    rnorm(1, mean = runif(1, 2.5, 15), sd = runif(1, 1.5, 25)),  # mu_2
    rnorm(1, mean = runif(1, 6.5, 10.5), sd = runif(1, 0.5, 10)),  # mu_3
    rnorm(1, mean = runif(1, 3.5, 15.5), sd = runif(1, 0.2, 3))  # mu_4
  )
  
  # Generate random positive definite Sigma with high variability
  repeat {
    df <- sample((p + 1):(p + 5), 1)  # Small range close to p for high variability
    Sigma_base <- toeplitz(runif(p, 0.5, 1))  # Random base matrix
    if (is.positive.definite(Sigma_base)) {
      Sigmas <- drop(rWishart(1, df, Sigma_base))
      if (is.positive.definite(Sigmas)) {
        break
      }
    }
  }
  
  list(mu = mu, Sigma = Sigmas)
}


# Generate 100 sets of random mu and Sigma
set.seed(123)  # For reproducibility
num_simulations <- 100
simulation_results <- lapply(1:num_simulations, function(i) {
  generate_parameters()
})


df1<- data.frame()
df2<- data.frame()
df3<- data.frame()
df4<- data.frame()

for(i in 1:num_simulations) {
  mu<- c( abs(simulation_results[[i]]$mu)[,1]  )
  Sigma<-simulation_results[[i]]$Sigma 
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
    
    print("ok")
    a<- i
    
    data_log<-  rlnorm.rplus(1e5,convert_mu_1,cov_convert_sigma_1)
    
    T0<-data_log[,1]
    T1<-data_log[,2]
    S0<-data_log[,3]
    S1<-data_log[,4]
    
    delta_T<- T1-T0
    delta_S<- S1-S0
    
    rho<- cor(delta_T, delta_S)
    ICA_normal_2<- rho^2
    
    mutual_info_2<- mutinfo(delta_T, delta_S, k=10, direct=FALSE)
    mutual_info_2<- mean(mutual_info_2)
    ICA_lognormal_2<- 1-exp(-2*mutual_info_2)
    
    diff<- ICA_lognormal_2- ICA_normal_2
    rel_diff<- (ICA_lognormal_2- ICA_normal_2)/ICA_lognormal_2
    
    f_multi_lognormal<- function(X) {
      x <- c(X[1], X[2], X[3],X[4])
      dlnorm.rplus(x,convert_mu_1,cov_convert_sigma_1) }
    
    f_multi_normal<- function(X) {
      x <- c(X[1], X[2], X[3],X[4])
      dmvnorm(x, mu, Sigma, log = FALSE, tol = 1e-06) }
    
    
    f_average<- function(x) (f_multi_lognormal(x)+f_multi_normal(x))/2
    
    lower_limit<- c(1,1,1,1)
    upper_limit = c((mu[1]+3*Sigma[1,1]),(mu[2]+3*Sigma[2,2]),(mu[3]+3*Sigma[3,3]),(mu[4]+3*Sigma[4,4]))
    volume<- (upper_limit[1]-lower_limit[1])*(upper_limit[2]-lower_limit[2])*(upper_limit[3]-lower_limit[3])*(upper_limit[4]-lower_limit[4])
    
    
    ##Kullback Leibler distance 1
    KL_fun1<- function(X) {
      f_multi_lognormal(X)* log(f_multi_lognormal(X)/f_multi_normal(X)) }
    
    
    
    ##Kullback Leibler distance 2
    KL_fun2<- function(X) {
      f_multi_normal(X)* log(f_multi_normal(X)/f_multi_lognormal(X)) }
    
    
    
    ##jensen functions
    js_function1<- function(X) { (f_multi_lognormal(X)* log(f_multi_lognormal(X)/f_average(X))) }
    js_function2<- function(X) { (f_multi_normal(X)* log(f_multi_normal(X)/f_average(X))) }
    
    
    kl1_values_2 <- data.frame()
    kl2_values_2 <- data.frame()
    js1_values_2 <- data.frame()
    js2_values_2 <- data.frame()
    
    for (j in 1:100) {
      
      ##monte carlo integrations
      N <- 10000
      x1_values<- runif(N, min=lower_limit[1], max= upper_limit[1])
      x2_values<- runif(N, min=lower_limit[2], max= upper_limit[2])
      x3_values<- runif(N, min=lower_limit[3], max= upper_limit[3])
      x4_values<- runif(N, min=lower_limit[4], max= upper_limit[4])
      
      kl1_values <- numeric(N)
      kl2_values <- numeric(N)
      js1_values <- numeric(N)
      js2_values <- numeric(N)  # Initialize a vector to store function values
      
      for (j in 1:N) {
        
        kl1_values[j] <- KL_fun1( c(x1_values[j], x2_values[j], x3_values[j], x4_values[j]))
        kl2_values[j] <- KL_fun2( c(x1_values[j], x2_values[j], x3_values[j], x4_values[j]))
        js1_values[j] <- js_function1( c(x1_values[j], x2_values[j], x3_values[j], x4_values[j]))
        js2_values[j] <- js_function2( c(x1_values[j], x2_values[j], x3_values[j], x4_values[j]))
        
        
        kl1_values_3 <- mean(na.omit(kl1_values[is.finite(kl1_values)]))*volume
        kl1_values_3 <- kl1_values_3[kl1_values_3 >= 0]
        
        kl2_values_3 <- mean(na.omit(kl2_values[is.finite(kl2_values)]))*volume
        kl2_values_3 <-  kl2_values_3[kl2_values_3 >= 0]
        
        js1_values_3 <- mean(na.omit(js1_values[is.finite(js1_values)]))*volume
        js1_values_3 <- js1_values_3[js1_values_3 >= 0]
        
        js2_values_3 <- mean(na.omit(js2_values[is.finite(js2_values)] ))*volume
        js2_values_3 <-  js2_values_3[js2_values_3 >= 0] }
      
      
      kl1_values_2 <- rbind(kl1_values_3,kl1_values_2 )
      kl2_values_2 <- rbind(kl2_values_3,kl2_values_2 )
      js1_values_2 <- rbind(js1_values_3,js1_values_2 )
      js2_values_2 <- rbind(js2_values_3,js2_values_2 )
      
    }
    
    
    js2_values_22 <- js2_values_2[js2_values_2 >= 0]
    js1_values_22 <- js1_values_2[js1_values_2 >= 0]
    
    
    KL1<- mean(as.matrix(kl1_values_2))
    KL2<- mean(as.matrix(kl2_values_2))
    js1<- mean(as.matrix(js1_values_22))
    js2<- mean(as.matrix(js2_values_22))
    jensen<- ( js1+ js2 )/2
    
    
    
    output1<- cbind(a, rho, ICA_normal_2, ICA_lognormal_2, diff, rel_diff,KL1, KL2, jensen)
    output2<- cbind(mu[1], mu[2], mu[3], mu[4],convert_mu_1[1],convert_mu_1[2],convert_mu_1[3],convert_mu_1[4])
    output3<- Sigma
    output4<- cov_convert_sigma_1
    
    df1<- rbind(df1,output1)
    df2<- rbind(df2, output2) 
    df3<- rbind(df3, output3) 
    df4<- rbind(df4, output4) 
    
  } }

end_time<- Sys.time()
print(end_time)


##save the data as csv
write.csv(df1,"~/Downloads/all_outputs_lognormal_4dec.csv", row.names = FALSE)
write.csv(df2,"~/Downloads/mu_all_4dec.csv", row.names = FALSE)
write.csv(df3,"~/Downloads/sigmas_lognormal_4dec.csv", row.names = FALSE)
write.csv(df4,"~/Downloads/converted_sigmas_lognormal_4dec.csv", row.names = FALSE)


