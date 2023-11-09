# this script generates Fig. 3.2


default <- c(0.9510, 0.9530, 0.5568, 0.5656)
i. <- c(0.9507, 0.9497, 0.5683, 0.5987)
iii. <- c(0.9520, 0.9560, 0.5585, 0.5674)
iv. <- c(0.9393, 0.9433, 0.5612, 0.5703)
v. <- c(0.9497, 0.9527, 0.5583, 0.5672)

rirs <- rbind(default, i., iii., iv., v.)

size <- 1.5

pdf("model_vio_rirs.pdf", width=7, height=7)

plot(x = rirs[,3], y = rirs[,1], xlim = c(0.5565, 0.599), ylim = c(0.939, 0.957), 
     col = c(1,2,4,5,6), pch = 1, cex = size, xlab = "Confidence interval width", 
     ylab = "Empirical coverage frequency", 
     main = "Model violation simulation results: \nRandom Coefficient Model vs. Aggregated Regressions")

points(x = rirs[,4], y = rirs[,2], col = c(1,2,4,5,6), pch = 19, cex = size)

thick <- 2

lines(x = c(rirs[1,3], rirs[1,4]), y = c(rirs[1,1], rirs[1,2]), col = 1, lwd = thick)
lines(x = c(rirs[2,3], rirs[2,4]), y = c(rirs[2,1], rirs[2,2]), col = 2, lwd = thick)
lines(x = c(rirs[3,3], rirs[3,4]), y = c(rirs[3,1], rirs[3,2]), col = 4, lwd = thick)
lines(x = c(rirs[4,3], rirs[4,4]), y = c(rirs[4,1], rirs[4,2]), col = 5, lwd = thick)
lines(x = c(rirs[5,3], rirs[5,4]), y = c(rirs[5,1], rirs[5,2]), col = 6, lwd = thick)

legend("bottomright", legend = c("RIRS", "Agg. Regr.", NA, "Default", "n_meas small", 
                                  "varying n_meas", "Exp(1) errors", 
                                 "AR(1) errors"), ncol = 1, pch = c(1, 19, NA, NA, NA, NA, NA, NA),
       pt.cex = size, text.col = c("azure4","azure4",NA,c(1,2,4,5,6)), col = c("azure4","azure4",NA,c(1,2,4,5,6)), 
       lty = c(NA, NA,NA, rep(1,5)), lwd = c(NA, NA,NA, rep(thick, 5)))

dev.off()


















