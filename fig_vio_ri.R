# this script generates Fig. 3.1

default <- c(0.9460, 0.9493, 0.0321, 0.0348)
i. <- c(0.9530, 0.9573, 0.0681, 0.1643)
ii. <- c(0.9497, 0.9507, 0.0721, 0.0849)
iii. <- c(0.9533, 0.9580, 0.0121, 0.0339)
iv. <- c(0.9503, 0.9457, 0.0640, 0.0694)
v. <- c(0.9500, 0.9527, 0.0507, 0.0552)

ri <- rbind(default, i., ii., iii., iv., v.)

size <- 1.5

pdf("model_vio_ri.pdf", width=7, height=7)

plot(x = ri[,3], y = ri[,1], xlim = c(0.01, 0.17), ylim = c(0.945, 0.958), 
     col = 1:6, pch = 1, cex = size, xlab = "Confidence interval width", 
     ylab = "Empirical coverage frequency", 
     main = "Model violation simulation results: \nRandom Intercept Model vs. Aggregated Regressions")

points(x = ri[,4], y = ri[,2], col = 1:6, pch = 19, cex = size)

thick <- 2

lines(x = c(ri[1,3], ri[1,4]), y = c(ri[1,1], ri[1,2]), col = 1, lwd = thick)
lines(x = c(ri[2,3], ri[2,4]), y = c(ri[2,1], ri[2,2]), col = 2, lwd = thick)
lines(x = c(ri[3,3], ri[3,4]), y = c(ri[3,1], ri[3,2]), col = 3, lwd = thick)
lines(x = c(ri[4,3], ri[4,4]), y = c(ri[4,1], ri[4,2]), col = 4, lwd = thick)
lines(x = c(ri[5,3], ri[5,4]), y = c(ri[5,1], ri[5,2]), col = 5, lwd = thick)
lines(x = c(ri[6,3], ri[6,4]), y = c(ri[6,1], ri[6,2]), col = 6, lwd = thick)

legend("bottomright", legend = c("RI", "Agg. Regr.", NA, "Default", "n_meas small", 
                                 "n_ind small", "varying n_meas", "Exp(1) errors", 
                                 "AR(1) errors"), ncol = 1, pch = c(1, 19, NA,NA, NA, NA, NA, NA, NA),
       pt.cex = size, text.col = c("azure4","azure4",NA,1:6), col = c("azure4","azure4",NA,1:6), 
       lty = c(NA, NA,NA, rep(1,6)), lwd = c(NA, NA,NA, rep(thick, 6)))

dev.off()


















