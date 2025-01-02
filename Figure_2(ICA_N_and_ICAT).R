##paper figure 2

install.packages("pracma")
library(pracma) ##for psi function

zeta_fun <- function(nu) 2*log( (gamma(nu/2)*gamma(1/2))/(sqrt(pi)*gamma( (1+nu)/2 ))*sqrt(nu/2))- (nu/2)*psi(nu/2)+(nu+1)*psi((1+nu)/2)- ((nu+2)/2)*psi((nu+2)/2)

##x=y line
x<- seq(0,1, by=0.01)
y<- seq(0,1, by=0.01)
reg<- lm(x ~y)


##\nu=3
nu<-3
zeta_fun(nu)

ICA_N3<- seq(0,1, by=0.2)
mut_n3<- -0.5*log(1-ICA_N3)
mut_t3<- mut_n3+ zeta_fun(nu)
ICA_t3<- 1-exp(-2*mut_t3)
##test
ICA_t3_2<- 1-exp(-2*mut_n3-2*zeta_fun(nu))

##\nu=4

nu<-4
ICA_N4<- seq(0,1, by=0.2)
mut_n4<- -0.5*log(1-ICA_N4)
mut_t4<- mut_n4+ zeta_fun(nu)
ICA_t4<- 1-exp(-2*mut_t4)
##test
ICA_t4_2<- 1-exp(-2*mut_n4-2*zeta_fun(nu))


##\nu=5

nu<-5
ICA_N5<- seq(0,1, by=0.2)
mut_n5<- -0.5*log(1-ICA_N5)
mut_t5<- mut_n5+ zeta_fun(nu)
ICA_t5<- 1-exp(-2*mut_t5)
##test
ICA_t5_2<- 1-exp(-2*mut_n5-2*zeta_fun(nu))

##\nu=7

nu<-7 
ICA_N7<- seq(0,1, by=0.2)
mut_n7<- -0.5*log(1-ICA_N7)
mut_t7<- mut_n7+ zeta_fun(nu)
ICA_t7<- 1-exp(-2*mut_t7)
##test
ICA_t7_2<- 1-exp(-2*mut_n7-2*zeta_fun(nu))


##plot
par(mfrow=c(2,2))

plot(ICA_N3, ICA_t3, type = "l", col="red", main=(expression(nu~"=3")), xlab = bquote(ICA[t]), ylab =( (bquote(ICA[N])))  )
abline(reg, col="black")

plot(ICA_N4, ICA_t4, type = "l", col="red", main=(expression(nu~"=4")), xlab = bquote(ICA[t]), ylab =( (bquote(ICA[N])))  )
abline(reg, col="black")

plot(ICA_N5, ICA_t5, type = "l", col="red", main=(expression(nu~"=5")), xlab = bquote(ICA[t]), ylab =( (bquote(ICA[N])))  )
abline(reg, col="black")

plot(ICA_N7, ICA_t7, type = "l", col="red", main=(expression(nu~"=7")), xlab = bquote(ICA[t]), ylab =( (bquote(ICA[N])))  )
abline(reg, col="black")


