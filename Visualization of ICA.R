##Visualization of ICA values in Schizophrenia data

##read ICA values
Ica_normal_schizo<-load("Ica_normal_schizo.rda")
Ica_normal_schizo<-get(Ica_normal_schizo)

Ica_t_schizo<-load("Ica_t_schizo.rda")
Ica_t_schizo<-get(Ica_t_schizo)

Ica_lognormal_schizo<-load("Ica_lognormal_schizo.rda")
Ica_lognormal_schizo<-get(Ica_lognormal_schizo)

##summarize 
s1<-summary(Ica_normal_schizo)
s2<-summary(Ica_t_schizo)
s3<-summary(Ica_lognormal_schizo)
summary_table <- bind_rows(s1, s2, s3)
print(summary_table)

##histograms 
hist(Ica_normal_schizo, xlim = c(0,1), ylim=c(0,200), main=bquote("Histogram of ICA"[N]), xlab=bquote("ICA"[N]))
grid()
box()

hist(Ica_t_schizo, xlim = c(0,1), ylim=c(0,200), main=bquote("Histogram of ICA"[t]), xlab=bquote("ICA"[t]))
grid()
box()



hist(Ica_lognormal_schizo, xlim = c(0,1), ylim=c(0,1200), main=bquote("Histogram of ICA"[L]), xlab=bquote("ICA"[L]))
grid()
box()

##densities
plot(density(Ica_normal_schizo),
     xlim = c(0,1),     # density için de sınırlama yapabilirsin
     ylim = c(0, 50), 
     main = "",
     xlab = "ICA", 
     col = "blue", lwd = 2, lty = 1)

lines(density(Ica_t_schizo), col = "red", lwd = 2, lty = 2)   # farklı çizgi tipi
lines(density(Ica_lognormal_schizo), col = "darkgreen", lwd = 2, lty = 1)

legend("topleft", 
       legend = c(expression(ICA[N]), expression(ICA[t]), expression(ICA[L])), 
       col = c("blue", "red", "darkgreen"), 
       lty = c(1,2,1), lwd = 2, bty = "n")
grid()

##histograms of ICA values in a single plot
hist(Ica_normal_schizo, breaks = 15, freq = FALSE, 
     xlim = c(0,1), ylim = c(0, 35),     # y ekseni daha uzun
     col = rgb(0,0.5,1,0.4),   
     main = "Histogram of ICA", 
     xlab = "ICA", ylab = "")

hist(Ica_t_schizo, breaks = 15, freq = FALSE, 
     col = rgb(1,0.5,0,0.4), add = TRUE)

hist(Ica_lognormal_schizo, breaks = 15, freq = FALSE, 
     col = rgb(0,0.6,0,0.3), add = TRUE)

legend("topleft", 
       legend = c(expression(ICA[N]), expression(ICA[t]), expression(ICA[L])), 
       fill = c(rgb(0,0.5,1,0.4), rgb(1,0.5,0,0.4), rgb(0,0.6,0,0.3)), 
       bty = "n")

grid()
box()

df <- data.frame(
  value = c(Ica_normal_schizo, Ica_t_schizo, Ica_lognormal_schizo),
  group = factor(rep(c("ICA_N", "ICA_t", "ICA_L"), 
                     times = c(length(Ica_normal), length(Ica_t_dist), length(ica_log_normal))))
)


ggplot(df, aes(x = value, fill = group)) +
  geom_histogram(bins = 15, alpha = 0.7, color = "black") +
  facet_wrap(~ group, nrow = 1, scales = "fixed") +
  xlim(0,1) +
  labs(title = "Histograms of ICA", x = "ICA", y = "") +
  theme_minimal()
# + theme(panel.border = element_rect(color = "black", fill = NA, size = 1)  # Adds a box )


