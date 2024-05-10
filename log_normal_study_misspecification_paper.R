##log-normal study for misspecification paper

##note:I am going to change/edit this codes

mu_1<- c(1.2, 4, 1.3, 1.5)
Sigma_1<- rbind(c(1.0780992, -0.3676667,  0.2402575, -0.4867786),
                c(-0.3676667,  1.0623604, -0.4322828,  0.3046945),
                c(0.2402575 ,-0.4322828 , 1.0594279, -0.2495095),
                c(-0.4867786 , 0.3046945, -0.2495095,  1.0704080))



# ##TO CONVERT PARAMETER
# convert_mu_1<- log( (mu_1)^2/ sqrt( (mu_1)^2+ diag(Sigma_1) ))
# 
# ##convert parameter
# var11<- log(1+ ( Sigma_1[1,1] / mu_1[1]^2) )
# var22<- log(1+ ( Sigma_1[2,2] / mu_1[2]^2) )
# var33<- log(1+ ( Sigma_1[3,3] / mu_1[3]^2) )
# var44<- log(1+ ( Sigma_1[4,4] / mu_1[4]^2) )
# var12<- log(1+ ( Sigma_1[1,2]/  (mu_1[1]*mu_1[2])) )
# var13<- log(1+ ( Sigma_1[1,3] / (mu_1[1]*mu_1[3])) )
# var14<- log(1+ ( Sigma_1[1,4] / (mu_1[1]*mu_1[4])) )
# var23<- log(1+ ( Sigma_1[2,3] / (mu_1[2]*mu_1[3])) )
# var24<- log(1+ ( Sigma_1[2,4] / (mu_1[2]*mu_1[4])) )
# var34<- log(1+ ( Sigma_1[3,4]/  (mu_1[3]*mu_1[4])) )
# 
# cov_convert_sigma_1<- matrix( c( var11, var12, var13, var14,
#                                  var12, var22, var23, var24,
#                                  var13, var23, var33, var34,
#                                  var14, var24, var34, var44), nrow = 4, ncol = 4)
# 

A = matrix(c(-1, 1, 0, 0,
             0, 0, -1, 1), ncol = 4, byrow = TRUE)

Sigma_Delta = A %*% Sigma_1 %*% t(A)
rho_Delta = Sigma_Delta[1, 2] / sqrt(Sigma_Delta[1, 1] * Sigma_Delta[2, 2])

# beta<-convert_mu_1[2]-convert_mu_1[1]
# alpha<-convert_mu_1[4]-convert_mu_1[3]

beta<-mu_1[2]-mu_1[1]
alpha<-mu_1[4]-mu_1[3]

# ##variances
# vT0T0<- cov_convert_sigma_1[1,1]
# vT1T1<- cov_convert_sigma_1[2,2]
# vS1S1<- cov_convert_sigma_1[4,4]
# vS0S0<- cov_convert_sigma_1[3,3]
# vT1S1<- cov_convert_sigma_1[2,4]
# vT0S0<- cov_convert_sigma_1[1,3]

##variances
vT0T0<- Sigma_1[1,1]
vT1T1<- Sigma_1[2,2]
vS1S1<- Sigma_1[4,4]
vS0S0<- Sigma_1[3,3]

#vT1S1<- Sigma_1[2,4]
#vT0S0<- Sigma_1[1,3]

##correlations 
# cT0T1<- cov_convert_sigma_1[1,2]/ (sqrt(vT0T0*vT1T1 ))
# cT0S1<- cov_convert_sigma_1[1,4]/ (sqrt(vT0T0*vS1S1 ))
# cT1S0<- cov_convert_sigma_1[2,3]/ (sqrt(vT1T1*vS0S0 ))
# cS0S1<- cov_convert_sigma_1[3,4]/ (sqrt(vS0S0*vS1S1 ))


##ESTİMATED correlations 
cT0T1<- Sigma_1[1,2]/ (sqrt(Sigma_1[1,1]*Sigma_1[2,2] ))
cT0S1<- Sigma_1[1,4]/ (sqrt(Sigma_1[1,1]*Sigma_1[4,4] ))
cT1S0<- Sigma_1[2,3]/ (sqrt(Sigma_1[2,2]*Sigma_1[3,3] ))
cS0S1<- Sigma_1[3,4]/ (sqrt(Sigma_1[3,3]*Sigma_1[4,4] ))


T0S0<-seq(-1, 1, by=.1)
T1S1<-seq(-1, 1, by=.1)
combins <- expand.grid(T0S0, T1S1)
i<-412
##to calculate ICA_normal using formula rho_delta using different metrics
for (i in 1: nrow(combins))  
{
  
  # pay<-(sqrt(vT0T0*vS0S0)*combins[i,1])+(sqrt(vT1T1*vS1S1)*combins[i,2])-(sqrt(vT1T1*vS0S0)*cT1S0)-(sqrt(vT0T0*vS1S1)*cT0S1)
  # payda<-sqrt((vT0T0+vT1T1-2*sqrt(vT0T0*vT1T1)*cT0T1)*(vS0S0+vS1S1-2*sqrt(vS0S0*vS1S1)*cS0S1))
  # ICAformul<-pay/payda
  
  
  pay<-(sqrt(Sigma_1[1,1]*Sigma_1[3,3])*combins[i,1])+(sqrt(Sigma_1[2,2]*Sigma_1[4,4])*combins[i,2])-(sqrt(Sigma_1[2,2]*Sigma_1[3,3])*cT1S0)-(sqrt(Sigma_1[1,1]*Sigma_1[4,4])*cT0S1)
  payda<-sqrt((Sigma_1[1,1]+Sigma_1[2,2]-2*sqrt(Sigma_1[1,1]*Sigma_1[2,2])*cT0T1)*(Sigma_1[3,3]+Sigma_1[4,4]-2*sqrt(Sigma_1[3,3]*Sigma_1[4,4])*cS0S1))
  ICAformul<-pay/payda
  
  
  sigmadelta_s<-Sigma_1[3,3]+Sigma_1[4,4]-(2*cS0S1*sqrt(Sigma_1[3,3]*Sigma_1[4,4]))
  sigmadelta_t<-Sigma_1[1,1]+Sigma_1[2,2]-(2*cT0T1*sqrt(Sigma_1[1,1]*Sigma_1[2,2]))
  
  # sigmadelta_s<-vS0S0+vS1S1- (2*cS0S1*sqrt(vS0S0*vS1S1))
  # sigmadelta_t<-vT0T0+vT1T1-2*cT0T1*sqrt(vT0T0*vT1T1)
  
  cormatrix<- matrix(c(1, cT0T1, combins[i, 1],cT0S1,
                       cT0T1,1,cT1S0,combins[i, 2],
                       combins[i, 1],cT1S0,1,cS0S1,
                       cT0S1,combins[i, 2],cS0S1, 1 ), nrow = 4, ncol=4)
  
  covmatrix<- matrix(c(vT0T0, cT0T1*sqrt(vT0T0*vT1T1), sigmadelta_t, cT0S1*sqrt(vT0T0*vS1S1),
                       cT0T1*sqrt(vT0T0*vT1T1),vT1T1,cT1S0*sqrt(vT1T1*vS0S0),sigmadelta_s,
                       sigmadelta_t,cT1S0*sqrt(vT1T1*vS0S0),vS0S0,cS0S1*sqrt(vS0S0*vS1S1),
                       cT0S1*sqrt(vT0T0*vS1S1), sigmadelta_s, cS0S1*sqrt(vS0S0*vS1S1), vS1S1 ), nrow = 4, ncol=4)
  
  covmatrix<- matrix(c(vT0T0, Sigma_1[1,2], sigmadelta_t, cT0S1*sqrt(vT0T0*vS1S1),
                       Sigma_1[1,2],vT1T1,cT1S0*sqrt(vT1T1*vS0S0),sigmadelta_s,
                       sigmadelta_t,cT1S0*sqrt(vT1T1*vS0S0),vS0S0,cS0S1*sqrt(vS0S0*vS1S1),
                       Sigma_1[1,4], sigmadelta_s, cS0S1*sqrt(vS0S0*vS1S1), vS1S1 ), nrow = 4, ncol=4)
  
  
  result<- is.positive.definite(cormatrix)
  result<- is.positive.definite(covmatrix)
  
  print(result) 
  print(i) 
  
  A = matrix(c(-1, 1, 0, 0,
               0, 0, -1, 1), ncol = 4, byrow = TRUE)
  
  cov_sigma_t_s<- ICAformul/( sqrt(sigmadelta_s*sigmadelta_t))
  Sigma_Delta_2<- matrix(c(sigmadelta_t, cov_sigma_t_s, cov_sigma_t_s, sigmadelta_s), nrow = 2, ncol = 2)
  
  Sigma_Delta_2 = A %*% covmatrix %*% t(A)
  Sigma_Delta_2[2,1]<- Sigma_Delta_2[1,2]
  is.positive.definite(Sigma_Delta_2)
  rho_Delta_2 = Sigma_Delta_2[1, 2] / sqrt(Sigma_Delta_2[1, 1] * Sigma_Delta_2[2, 2])
  
  if (result==TRUE){  write(ICAformul, file="log_normal_ICA_normal.txt", append=TRUE) }
  
}


alphas_betas<- c(beta, alpha)
alphas_betas<- c(mu_1[2]-mu_1[1], mu_1[4]-mu_1[3])

convert_alpha_beta<- log( (alphas_betas)^2/ sqrt( (alphas_betas)^2+ diag(Sigma_Delta_2) ))

##convert parameter
var11<- log(1+ (  Sigma_Delta_2[1,1] / alphas_betas[1]^2) )
var22<- log(1+ (  Sigma_Delta_2[2,2] / alphas_betas[2]^2) )
var12<- log(1+ (  Sigma_Delta_2[1,2] / (alphas_betas[1]*alphas_betas[2])) )

var11<- log(1+ (  Sigma_Delta_2[1,1] / beta^2) )
var22<- log(1+ (  Sigma_Delta_2[2,2] / alpha^2) )
var12<- log(1+ (  Sigma_Delta_2[1,2] / (alpha*beta)) )

convert_sigma_delta<- matrix(c(var11, var12, var12, var22), nrow = 2, ncol=2)
is.positive.definite(convert_sigma_delta)

data_1<- rlnorm.rplus(1e5,convert_alpha_beta,convert_sigma_delta)

data_1<- rlnorm.rplus(1e5,log(alphas_betas),log(Sigma_Delta_2))
delta_T<- data_1[,1]
delta_S<- data_1[,2]


# cT0S0<- vT0S0/(sqrt(vT0T0)*sqrt(vS0S0))
# cT1S1<-vT1S1/(sqrt(vT1T1)*sqrt(vS1S1))
# cTS<- rho_Delta

data_1<- rlnorm.rplus(3000,convert_mu_1,cov_convert_sigma_1)
T0<-data_1[,1]
T1<-data_1[,2]
S0<-data_1[,3]
S1<-data_1[,4]
treatment<-rbinom(3000,1,0.5)
data<-data.frame(T0,T1,S0,S1, treatment)


##select 2 coloumn for observed data
newdata<-matrix(NA,nrow= 3000, ncol=3)
for (i in 1:3000) {
  if(data$treatment[i]==0)
  {
    newdata[i,]<-c(data$T0[i], data$S0[i], data$treatment[i])
  }
  else
    newdata[i,]<-cbind(data$T1[i], data$S1[i], data$treatment[i])
}

newdata<-as.data.frame(newdata)
placebo<-newdata[newdata$V3 == "0",]
experimental<-newdata[newdata$V3  == "1",]

##basic measures
Tr<-newdata$V1
S<-newdata$V2
T0<- placebo$V1
T1<- experimental$V1
S0<- placebo$V2
S1<- experimental$V2
var(Tr,Tr)
##means
mT0<-mean(T0)
mT1<-mean(T1)
mS0<-mean(S0)
mS1<-mean(S1)

beta<-mT1-mT0
alpha<-mS1-mS0

##variances
vT0T0<-var(T0,T0)
vT1T1<-var(T1,T1)
vS1S1<-var(S1,S1)
vS0S0<-var(S0,S0)
vT1S1<-var(T1,S1)
vT0S0<-var(T0,S0)

##correlations
cT0S0<-cor(T0,S0)
cT1S1<-cor(T1,S1)
cTS<-cor(Tr,S)

##ICA normal using formula rho_delta 
payhigh<-(sqrt(vT0T0*vS0S0)*cT0S0)+(sqrt(vT1T1*vS1S1)*cT1S1)-(sqrt(vT1T1*vS0S0)*metric_1[,4])-(sqrt(vT0T0*vS1S1)*metric_1[,3])
paydahigh<-sqrt((vT0T0+vT1T1-2*sqrt(vT0T0*vT1T1)*metric_1[,1])*(vS0S0+vS1S1-2*sqrt(vS0S0*vS1S1)*metric_1[,6]))
ICAformulhigh<-payhigh/paydahigh

##variances
vT0T0<- Sigma_1[1,1]
vT1T1<- Sigma_1[2,2]
vS1S1<- Sigma_1[4,4]
vS0S0<- Sigma_1[3,3]
vT1S1<- Sigma_1[2,4]
vT0S0<- Sigma_1[1,3]

##correlations
cT0S0<- vT0S0/(sqrt(vT0T0)*sqrt(vS0S0))
cT1S1<-vT1S1/(sqrt(vT1T1)*sqrt(vS1S1))
cTS<- rho_Delta


sigmadelta_s<-var(S0,S0)+var(S1,S1)-2*metric_1[,6]*sqrt(var(S0,S0)*var(S1,S1))
sigmadelta_t<-var(T0,T0)+var(T1,T1)-2*metric_1[,1]*sqrt(var(T0,T0)*var(T1,T1))
cormatrixhigh<- matrix(c(1, metric_1[,1], cT0S0,metric_1[,3],
                         metric_1[,1],1,metric_1[,4],cT1S1,
                         cT0S0,metric_1[,4],1,metric_1[,6],
                         metric_1[,3],cT1S1,metric_1[,6], 1 ), nrow = 4)


var11<-log(1+ ( sigmadelta_t / (beta)^2) )
var22<-log(1+ ( sigmadelta_s / (alpha)^2) )
var12<- log(1+ ( (cormatrixhigh[2,3]*(sqrt(var11*var22))) /  (beta+alpha)^2) )

sigma_2_covert<- matrix( c(var11, var12, var12, var22), nrow=2, ncol=2)

mean_T<- log( (beta)^2/ sqrt( (beta)^2+ sigmadelta_t ))
mean_S<- log( (alpha)^2/ sqrt( (alpha)^2+ sigmadelta_s ))
mu_2_covert<- cbind(mean_T, mean_S)

data_delta<- rlnorm.rplus(1e5, mu_2_covert ,sigma_2_covert )
delta_T<- data_delta[,1]
delta_S<- data_delta[,2]


delta_T<- rlnorm.rplus(1e5, mean_T , var11)
delta_S<- rlnorm.rplus(1e5, mean_S , var22)

mean(delta_T)
var(delta_T)
mean(delta_S)
var(delta_S)

plot(density(delta_T))
plot(density(delta_S))


mutual_info<- mutinfo(delta_T, delta_S, k=10, direct = FALSE)
mutual_info_2<- mutinfo(delta_T, delta_S, k=10, direct = TRUE)
mutual_info<- mean(mutual_info)

1-exp(-2*mutual_info)
1-exp(-2*mutual_info_2)

