##Figure 3


##reading data
data_lognormal<- read.table(file="data_lognormal.txt")
colnames(full_data_lognormal) <- c("a", "T0", "T1", "S0", "S1")


##reading output
lognormal_output<- read.table(file="log_normal_output.txt")
colnames(log_normal_output) <- c("a","rho", "ICA_N", "ICA_L", "diff")
log_normal_output<- log_normal_output[1:200,]


##merge data and output
full<- merge(data_lognormal, lognormal_output, by="a")

###density plots example:
diff_min<- which(full$a==79) 
diff_min_data<- full[diff_min,]

##density plot of delta T and delta S
par(mfrow=c(2,2))
plot(density(diff_min_data$T0), xlab = "", main=(expression(T[0])), ylab = "density")
plot(density(diff_min_data$T1), xlab = "",  main=(expression(T[1])), ylab = "density")
plot(density(diff_min_data$S0), xlab = "",  main=(expression(S[0])), ylab = "density")
plot(density(diff_min_data$S1),  xlab = "",  main=(expression(S[1])), ylab = "density")
mtext( expression(bold("d= -0.00079")),side = 3, line = - 1.5, outer = TRUE)



##for differences graph (Figure 3 d=ICA_L-ICA_N)
##add numbers
numbers<- seq(1, 200, by=1)
data_full_unique<- cbind(numbers, log_normal_output)
data_full_unique<- as.data.frame(data_full_unique)
x<-  as.numeric(data_full_unique2$numbers)

par(mfrow=c(1,1))
plot(x, log_normal_output$diff, pch=20, main= bquote("d=" ~"ICA"[L] ~ "-" ~ "ICA"[N]), col="blue", ylab=bquote("ICA"[L] ~ "-" ~ "ICA"[N]), xlab=bquote(mu ~ "and"~  Sigma ~ "parameter pairs"))
abline(h=0, lty=2)

max(log_normal_output$diff,)
min(log_normal_output$diff,)
summary(log_normal_output$diff,)




