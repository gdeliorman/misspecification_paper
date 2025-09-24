library(Surrogate) ##to get data and to calculate ICA normal
library(LaplacesDemon) ##to check positive definite matrix
library(dplyr) ##for merge
library(utils) #for write csv
library(fitHeavyTail) ##to estimate degree of freedom for PANSS and BPRS

 
data("Schizo")
##NA rows
Schizo<- Schizo[-c(405,  705, 1358, 1719, 2111),]
placebo<-Schizo[Schizo$Treat == "-1",]
experimental<-Schizo[Schizo$Treat == "1",]

True_endpoint<-Schizo$PANSS+103
Surrogate_endp<-Schizo$BPRS+56

T0<- placebo$PANSS
T1<- experimental$PANSS
S0<- placebo$BPRS
S1<- experimental$BPRS


##mean
mT0<-mean(T0)
mT1<-mean(T1)
mS0<-mean(S0)
mS1<-mean(S1)

beta<-mT1-mT0
alpha<- mS1-mS0
muu<- c(mT0, mT1, mS0, mS1)
mudelta<-c(beta,alpha)

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
cTS<-cor(True_endpoint,Surrogate_endp)

##check normality:

##histograms
hist(True_endpoint, # histogram
     col="peachpuff", # column color
     border="black",
     prob = TRUE, # show densities instead of frequencies
     xlab = "True endpoint", xlim = c(-100, 100),
     main = "")
lines(density(True_endpoint), # density plot
      lwd = 2, # thickness of line
      col = "chocolate3")

hist(True_endpoint, # histogram
     col="peachpuff", # column color
     border="black",
     prob = TRUE, # show densities instead of frequencies
     xlab = "PANSS", xlim = c(-100, 100),
     main = "")
lines(density(True_endpoint), # density plot
      lwd = 2, # thickness of line
      col = "chocolate3")

hist(Surrogate_endp, # histogram
     col="peachpuff", # column color
     border="black",
     prob = TRUE, # show densities instead of frequencies
     xlab = "Surrogate endpoint", xlim = c(-70, 70),
     main = "")
lines(density(Surrogate_endp), # density plot
      lwd = 2, # thickness of line
      col = "chocolate3")

hist(Surrogate_endp, # histogram
     col="peachpuff", # column color
     border="black",
     prob = TRUE, # show densities instead of frequencies
     xlab = "BPRS", xlim = c(-70, 70),
     main = "")
lines(density(Surrogate_endp), # density plot
      lwd = 2, # thickness of line
      col = "chocolate3")


hist(T1, # histogram
     col="peachpuff", # column color
     border="black",
     prob = TRUE, # show densities instead of frequencies
     xlab = "T1", xlim = c(-100, 100),
     main = "")
lines(density(T1), # density plot
      lwd = 2, # thickness of line
      col = "chocolate3")


hist(T0, # histogram
     col="peachpuff", # column color
     border="black",
     prob = TRUE, # show densities instead of frequencies
     xlab = "T0", xlim = c(-100, 100),
     main = "")
lines(density(T0), # density plot
      lwd = 2, # thickness of line
      col = "chocolate3")


box()
grid()

hist(S0, # histogram
     col="peachpuff", # column color
     border="black",
     prob = TRUE, # show densities instead of frequencies
     xlab = "S0", xlim = c(-70, 70),
     main = "")
lines(density(S0), # density plot
      lwd = 2, # thickness of line
      col = "chocolate3")


box()
grid()

hist(S1, # histogram
     col="peachpuff", # column color
     border="black",
     prob = TRUE, # show densities instead of frequencies
     xlab = "S1", xlim = c(-70, 70),
     main = "")
lines(density(S1), # density plot
      lwd = 2, # thickness of line
      col = "chocolate3")


box()
grid()

##Approximately normal


##qq plots
#par(mfrow = c(1, 2))
qqnorm(True_endpoint, main='True endpoint')
qqnorm(True_endpoint, main='PANSS')
qqline(True_endpoint)

qqnorm(Surrogate_endp, main='Surrogate endpoint')
qqnorm(Surrogate_endp, main='BPRS')
qqline(Surrogate_endp)

qqnorm(T0, main='T0')
qqline(T0)

qqnorm(T1, main='T1')
qqline(T1)

qqnorm(S0, main='S0')
qqline(S0)

qqnorm(S1, main='S1')
qqline(S1)


##Shapiro-Wilk normality test
shapiro.test(True_endpoint)
shapiro.test(Surrogate_endp)

shapiro.test(T0)
shapiro.test(T1)
shapiro.test(S0)
shapiro.test(S1)

new_Schizo<- cbind(Schizo$PANSS, Schizo$BPRS)
ks.test(new_Schizo, "pnorm")
ks.test(new_Schizo, "plnorm")


ks.test(True_endpoint, "pnorm")
ks.test(Surrogate_endp, "pnorm")
ks.test(True_endpoint,Surrogate_endp)
ks.test(True_endpoint,Surrogate_endp, "plnorm")



ks.test(T0, "pnorm")
ks.test(T1, "pnorm")
ks.test(S0, "pnorm")
ks.test(S1, "pnorm")

ks.test(True_endpoint, "lnorm")
ks.test(Surrogate_endp, "lnorm")

library(EnvStats)
gofTest(T0, distribution = "lnorm")


library(MASS)
library(logspline)
library(fitdistrplus)

fitdistr(True_endpoint, "log-normal")
fitdistr(True_endpoint, "lognormal")
fitdistr(True_endpoint, "normal")

descdist(True_endpoint, discrete = FALSE)
descdist(Surrogate_endp, discrete = FALSE)

descdist(T0, discrete = FALSE)
descdist(T1, discrete = FALSE)
descdist(S0, discrete = FALSE)
descdist(S1, discrete = FALSE)


fitdistr(True_endpoint, "t")
fitdistr(Surrogate_endp, "t")

nn<- rnorm(4000, mean=4, sd=2)
descdist(nn, discrete = FALSE)
ks.test(nn, "pnorm")
shapiro.test(nn)

ll<-rlnorm(1e4, meanlog = 0, sdlog = 1)
descdist(ll, discrete = FALSE)

ks.test(True_endpoint, Surrogate_endp)
ks.test(new_Schizo, "pnorm")


##ICA under normality

##ICA grid={-1, -0.95,...} 41^4=2825761 6143 pdf matrices
#ICA<-ICA.ContCont(cT0S0, cT1S1, vT0T0, vT1T1, vS0S0, vS1S1, T0T1=seq(-1, 1, by=.05), 
      #           T0S1=seq(-1, 1, by=.05), T1S0=seq(-1, 1, by=.05), S0S1=seq(-1, 1, by=.05))

#ICA grid2={-1, -0.90,...} 21^4=194481 343 pdf matrices 
ICA_2<-ICA.ContCont(cT0S0, cT1S1, vT0T0, vT1T1, vS0S0, vS1S1, T0T1=seq(-1, 1, by=.1), 
                 T0S1=seq(-1, 1, by=.1), T1S0=seq(-1, 1, by=.1), S0S1=seq(-1, 1, by=.1))

ICA<- ICA_2

rho_delta<- ICA$ICA
Ica_normal<- rho_delta^2

mean(Ica_normal)
max(Ica_normal)
min(Ica_normal)
s1<-summary(Ica_normal)
s1
quantile(Ica_normal)


hist(Ica_normal, xlim = c(0,1), ylim=c(0,200), main=bquote("Histogram of ICA"[N]), xlab=bquote("ICA"[N]))
grid()
box()

##ICA under t-distributed 
new_Schizo<- cbind(Schizo$PANSS, Schizo$BPRS)
sigma_t<- cov(new_Schizo)

##estimate degree of freedom \nu
fit_log<-fit_mvt(new_Schizo)
#fit_mvt(new_Schizo, nu = "iterative")
nu<- fit_log$nu

library(fitHeavyTail)
nu_2<- nu_POP_estimator(Xc = new_Schizo, nu = 8.12, Sigma = sigma_t)
#library(MixMatrix)

Ica_t<- ICA_t(
  #df=nrow(Schizo), if it is too big then ICA_t=ICA_N
  df=nu_2,
  T0S0=cT0S0,
  T1S1=cT1S1,
  T0T0 = vT0T0,
  T1T1 = vT1T1,
  S0S0 = vS0S0,
  S1S1 = vS1S1,
  T0T1 = seq(-1, 1, by = 0.1),
  T0S1 = seq(-1, 1, by = 0.1),
  T1S0 = seq(-1, 1, by = 0.1),
  S0S1 = seq(-1, 1, by = 0.1)
)


Ica_t_dist<- Ica_t$ICA_t

mean(Ica_t_dist)
max(Ica_t_dist)
min(Ica_t_dist)
s2<-summary(Ica_t_dist)
s2
quantile(Ica_t_dist)
sd(Ica_t_dist)

hist(Ica_t_dist, xlim = c(0,1), ylim=c(0,200), main=bquote("Histogram of ICA"[t]), xlab=bquote("ICA"[t]))
grid()
box()

##Relative effect from the data:
RE<- beta/alpha
RE

##Adjusted Assosiaction (AA) 
model_s<- lm(BPRS  ~ Treat, data=Schizo)
model_t<- lm(PANSS ~ Treat, data=Schizo)
beta_model<- model_t$coefficients[2]
alpha_model<- model_s$coefficients[2]
RE_model<-  as.numeric(beta_model/alpha_model)
RE_model

res_t<- model_t$residuals
res_s<- model_s$residuals
AA_model<- var(res_t, res_s)/sqrt(var(res_t)*var(res_s))

##AA= gamma= sigma_st/sqrt(sigma_ss*sigma_tt)
Sur <- Single.Trial.RE.AA(Dataset=Schizo, Surr=BPRS, True=PANSS, Treat=Treat, Pat.ID=Id)
summary(Sur)
plot(Sur)

AA<- Sur$AA
mean(Sur$AA.Boot.Samples)

##ICA log normal: 
mu<- cbind(mT0, mT1, mS0, mS1)
mu_2<- abs(mu)
Sigma<- matrix(c(vT0T0, NA, vT0S0, NA,
                 NA, vT1T1, NA, vT1S1,
                 vT0S0,NA, vS0S0, NA,
                 NA, vT1S1, NA, vS1S1), nrow = 4, ncol=4)

convert_mu_1<- log( (mu)^2/ sqrt( (mu)^2+ diag(Sigma) ))

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

T0T1<- seq(-1, 1, by=.05)
T0S1<- seq(-1, 1, by=.05)
T1S0<- seq(-1, 1, by=.05)
S0S1<- seq(-1, 1, by=.05)
T0T0<- cov_convert_sigma_1[1,1]
T1T1<- cov_convert_sigma_1[2,2]
S0S0<- cov_convert_sigma_1[3,3]
S1S1<- cov_convert_sigma_1[4,4]

combins <- expand.grid(T0T1, T0S1, T1S0, S0S1)
library(compositions) 
library(FNN)
for (i in 1: nrow(combins)) {

cov_convert_sigma_1[1,2]<-cov_convert_sigma_1[2,1]<- combins[i,1] * (sqrt(T0T0)*sqrt(T1T1))
cov_convert_sigma_1[4,1]<-cov_convert_sigma_1[1,4]<- combins[i,2] * (sqrt(T0T0)*sqrt(S1S1))
cov_convert_sigma_1[3,2]<-cov_convert_sigma_1[2,3]<- combins[i,3] * (sqrt(S0S0)*sqrt(T1T1))
cov_convert_sigma_1[3,4]<-cov_convert_sigma_1[4,3]<- combins[i,4] * (sqrt(S0S0)*sqrt(S1S1)) 

Sigma[1,2]<-Sigma[2,1]<- combins[i,1] * (sqrt(Sigma[1,1])*sqrt(Sigma[2,2]))
Sigma[4,1]<-Sigma[1,4]<- combins[i,2] * (sqrt(Sigma[1,1])*sqrt(Sigma[4,4]))
Sigma[3,2]<-Sigma[2,3]<- combins[i,3] * (sqrt(Sigma[3,3])*sqrt(Sigma[2,2]))
Sigma[3,4]<-Sigma[4,3]<- combins[i,4] * (sqrt(Sigma[3,3])*sqrt(Sigma[4,4])) 

if(is.positive.definite(cov_convert_sigma_1) == TRUE){
print(i)
data_log<-  rlnorm.rplus(1e5,convert_mu_1,cov_convert_sigma_1)


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
a<-i

diff<- ICA_lognormal_2- ICA_normal_2
outputs<- cbind(a, rho, ICA_normal_2, ICA_lognormal_2, diff)
write.table(outputs, file="log_normal_output_4aug.txt",  sep = "\t", row.names = FALSE, col.names = FALSE, append=TRUE)


}

}

data_lognormal<- read.table(file="log_normal_output_4aug.txt")
data_lognormal<- read.table(file="/Users/gokcedeliorman/log_normal_output_4aug.txt")
ica_log_normal<- as.numeric(unlist(data_lognormal[4]))
hist(ica_log_normal, xlim = c(0,1), ylim=c(0,1200), main=bquote("Histogram of ICA"[L]), xlab=bquote("ICA"[L]))
grid()
box()
summary(ica_log_normal)
mean(ica_log_normal)
min(ica_log_normal)
max(ica_log_normal)

quantile(ica_log_normal)
sd(ica_log_normal)
s3<-summary(ica_log_normal)
s3

##compare ICA_N and ICA_t
boxplot(Ica_normal ,Ica_t_dist, ica_log_normal, AA_model,
        names = c(expression(ICA[N]), expression(ICA[t]), expression(ICA[L]), expression(AA (gamma)) ),  ylim=c(0,1)
        #, border = "brown", col="orange")
)

PTE_roboust<- 0.864
PTE_model<- 0.849
PTE_freedman<- 0.847
PTE<- 0.849

boxplot(Ica_normal ,Ica_t_dist, ica_log_normal, AA_model, PTE,
        names = c(expression(ICA[N]), expression(ICA[t]), expression(ICA[L]), expression(AA (gamma)), "PTE" ),  ylim=c(0,1)
        #, border = "brown", col="orange")
)

summary_table <- bind_rows(s1, s2, s3)
print(summary_table)

###all results histograms and density plots
hist(Ica_t_dist, xlim = c(0,1), ylim=c(0,200), main=bquote("Histogram of ICA"[t]), xlab=bquote("ICA"[t]))
grid()
box()

hist(Ica_normal, xlim = c(0,1), ylim=c(0,200), main=bquote("Histogram of ICA"[N]), xlab=bquote("ICA"[N]))
grid()
box()

hist(ica_log_normal, xlim = c(0,1), ylim=c(0,1200), main=bquote("Histogram of ICA"[L]), xlab=bquote("ICA"[L]))
grid()
box()


hist(Ica_normal, breaks = 15, freq = FALSE, 
     xlim = c(0,1), ylim = c(0, 35),     # y ekseni daha uzun
     col = rgb(0,0.5,1,0.4),   
     main = "Histogram of ICA", 
     xlab = "ICA", ylab = "")

hist(Ica_t_dist, breaks = 15, freq = FALSE, 
     col = rgb(1,0.5,0,0.4), add = TRUE)

hist(ica_log_normal, breaks = 15, freq = FALSE, 
     col = rgb(0,0.6,0,0.3), add = TRUE)

legend("topleft", 
       legend = c(expression(ICA[N]), expression(ICA[t]), expression(ICA[L])), 
       fill = c(rgb(0,0.5,1,0.4), rgb(1,0.5,0,0.4), rgb(0,0.6,0,0.3)), 
       bty = "n")

grid()
box()


par(mfrow = c(1,1))
df <- data.frame(
  value = c(Ica_normal, Ica_t_dist, ica_log_normal),
  group = factor(rep(c("ICA_N", "ICA_t", "ICA_L"), 
                     times = c(length(Ica_normal), length(Ica_t_dist), length(ica_log_normal))))
)

# Facet histogram
ggplot(df, aes(x = value, fill = group)) +
  geom_histogram(bins = 15, alpha = 0.7, color = "black") +
  facet_wrap(~ group, nrow = 1, scales = "fixed") +
  xlim(0,1) +
  labs(title = "Histograms of ICA", x = "ICA", y = "") +
  theme_minimal()


plot(density(Ica_normal),
     xlim = c(0,1),     # density için de sınırlama yapabilirsin
     ylim = c(0, 50), 
     main = "",
     xlab = "ICA", 
     col = "blue", lwd = 2, lty = 1)

lines(density(Ica_t_dist), col = "red", lwd = 2, lty = 2)   # farklı çizgi tipi
lines(density(ica_log_normal), col = "darkgreen", lwd = 2, lty = 1)

legend("topleft", 
       legend = c(expression(ICA[N]), expression(ICA[t]), expression(ICA[L])), 
       col = c("blue", "red", "darkgreen"), 
       lty = c(1,2,1), lwd = 2, bty = "n")
grid()
box()

###all results saving
save(Ica_normal, file = "Ica_normal_schizo.rda")
Ica_normal_schizo<-load("Ica_normal_schizo.rda")
Ica_normal_schizo<-get(Ica_normal_schizo)

save(Ica_t_dist, file = "Ica_t_schizo.rda")
Ica_t_schizo<-load("Ica_t_schizo.rda")
Ica_t_schizo<-get(Ica_t_schizo)

save(ica_log_normal, file = "Ica_lognormal_schizo.rda")
Ica_lognormal_schizo<-load("Ica_lognormal_schizo.rda")
Ica_lognormal_schizo<-get(Ica_lognormal_schizo)




