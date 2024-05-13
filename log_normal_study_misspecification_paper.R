##log-normal study for misspecification paper

##note:I am going to change/edit this codes

library(LaplacesDemon) ##to use is.positive.definite() function


start_time<- Sys.time()
print(start_time)


##generate random mu 
mu_1<-rnorm(100, mean = 1.5, sd = 0.5)
mu_2<-rnorm(100, mean = 4, sd = 2)
mu_3<-rnorm(100, mean = 7, sd = 0.75)
mu_4<-rnorm(100, mean = 10, sd = 0.25)


##generate random sigmas
n <- 100  # Number of matrices
p <- 4  # Dimension
df <- 10  # Degrees of freedom
Sigma <- toeplitz((p:1)/p)  # the matrix parameter of the distribution
# Draw n Wishart distributed matrices
Sigmas <- drop(rWishart(200, df, Sigma))



for(i in 1:500) {
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
    
    ICA_lognormal<- list()
    ICA_normal<- list()
    
    for( j in 1:3000) {
      data_log<-  rlnorm.rplus(1e5,convert_mu_1,cov_convert_sigma_1)
      
      T0<-data_log[,1]
      T1<-data_log[,2]
      S0<-data_log[,3]
      S1<-data_log[,4]
      
      delta_T<- T1-T0
      delta_S<- S1-S0
      
      rho<- cor(delta_T, delta_S)
      ICA_normal_2<- rho^2
      ICA_normal<- append(ICA_normal_2, ICA_normal)
      
      #mutual_info<-   mutinfo(delta_T, delta_T, k=10, direct=TRUE)
      mutual_info_2<- mutinfo(delta_T, delta_S, k=10, direct=FALSE)
      mutual_info_2<- mean(mutual_info_2)
      
      #ICA_lognormal_1<- 1-exp(-2*mutual_info)
      ICA_lognormal_2<- 1-exp(-2*mutual_info_2)
      ICA_lognormal<- append( ICA_lognormal_2, ICA_lognormal)
      write(mean(as.numeric(rho)), file="rho_all_09november.txt", append=TRUE)
      write(mean(as.numeric(ICA_normal_2)), file="ica_normal_2_all.txt", append=TRUE)
      write(mean(as.numeric(ICA_lognormal_2)), file="ica_lognormal_2_all.txt", append=TRUE)
      
      
    } }
  
  diff<- ICA_lognormal_2- ICA_normal_2
  
  
  write(mean(as.numeric(ICA_normal)), file="log_normal_ica.txt", append=TRUE)
  write(mean(as.numeric(ICA_lognormal)), file="normal_ica.txt", append=TRUE)
  
}


end_time<- Sys.time()
print(end_time)


ICAs_lognormals<-read.delim(file= '/Users/gokcedeliorman/log_normal_ica.txt', header = FALSE, sep = "\t", dec = ".")
ICAs_normals<-read.delim(file= '/Users/gokcedeliorman/normal_ica.txt', header = FALSE, sep = "\t", dec = ".")

