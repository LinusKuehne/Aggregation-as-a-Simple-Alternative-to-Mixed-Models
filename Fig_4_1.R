# here, we fit a RI model on data generated from a RIRS model
# this is used to generate Fig. 4.1

library(lme4)
library(dplyr)


# structure of the data
n.ind <- 50
n.meas <- 10

# population intercept
a <- 1

# population slope
b <- 2

# error distribution
sd.error <- 0.5

# random intercept distribution
sd.ri <- 2

# random slope distribution
sd.rs.vec <- seq(from= 0, to = 1, by= 0.02)

# number of simulations
B <- 500


set.seed(1)


# function to simulate data from a RIRS model
RIRS.simulate <- function(n.ind, n.meas, a, b, sd.error, sd.ri, sd.rs){
  n.total <- n.ind * n.meas
  
  #structure of the data set
  individual <- rep(1:n.ind, each = n.meas)
  
  # days of sleep deprivation 
  x <- runif(n=n.total, min = 0, max = 10)
  
  # errors
  error <- rnorm(n.total, mean = 0, sd = sd.error)
  
  # random intercept
  gamma.0 <- rnorm(n.ind, mean = 0, sd = sd.ri)
  
  # random slope
  gamma.1 <- rnorm(n.ind, mean = 0, sd = sd.rs)
  
  # linear RI model
  y <- a + b*x + gamma.1[individual]*x + gamma.0[individual] + error
  d <- data.frame(x,y,individual)
  
  return(d)
}



# single simulation run 
# generates data set, fits both models on same data set, and computes relevant estimates
combined.sim.run <- function(n.ind, n.meas, a, b, sd.error, sd.ri, sd.rs){
  
  data <- RIRS.simulate(n.ind, n.meas, a, b, sd.error, sd.ri, sd.rs)
  
  
  # RI model ##########################################################
  RI.model <- lmer(y ~ x + (1 | individual), data = data)
  
  conf.int.ri <- suppressMessages(confint(RI.model, level = 0.95))
  left.ri <- conf.int.ri[4,1]
  right.ri <- conf.int.ri[4,2]
  
  mean.ri <- summary(RI.model)$coeff[2,1]
  #####################################################################
  
  c.ri <- 0
  if((left.ri<=b) && (right.ri>=b)){c.ri <- 1}
  
  
  
  
  # aggregated OLS ######################################################
  slopes <- numeric(n.ind)
  
  # linear model for each individual individually
  for(i in 1:n.ind){
    data.i <- data[data$individual == i,]
    model.i <- lm(y ~ x, data = data.i)
    slope.i <- model.i$coefficients[2]
    slopes[i] <- slope.i
  }
  
  test <- t.test(slopes, conf.level = 0.95)
  
  mean.agg <- test$estimate[1]
  
  left.agg <- test$conf.int[1]
  right.agg <- test$conf.int[2]
  #######################################################################
  
  c.agg <- 0
  if((left.agg<=b) && (right.agg>=b)){c.agg <- 1}
  
  
  
  output <- data.frame(c(c.ri, c.agg), c(left.ri, left.agg), c(right.ri, right.agg), 
                       c(mean.ri, mean.agg), row.names = c("ri", "agg"))
  return(output)
}





#store the results
results.run <- data.frame(width.ri = sd.rs.vec, width.agg = sd.rs.vec, 
                          coverage.ri = sd.rs.vec, coverage.agg = sd.rs.vec)





for(j in 1:length(sd.rs.vec)){
  sd.rs <- sd.rs.vec[j]
  
  
  # counts how many times the conf.int covered the mean
  count.ri <- 0
  count.agg <- 0
  
  # store the widths
  width.ri <- numeric(B)
  width.agg <- numeric(B)
  
  
  
  for(i in 1:B){
    res.df <- combined.sim.run(n.ind, n.meas, a, b, sd.error, sd.ri, sd.rs)
    
    cover.ri <- res.df[1,1]
    count.ri <- count.ri + cover.ri
    
    cover.agg <- res.df[2,1]
    count.agg <- count.agg + cover.agg
    
    width.ri[i] <- res.df[1,3] - res.df[1,2]
    
    width.agg[i] <- res.df[2,3] - res.df[2,2]
    
    
    # progress update
    if(i%%100 == 0){print(i+j)}
  }
  
  
  results.run$width.ri[j] <- mean(width.ri)
  results.run$width.agg[j] <- mean(width.agg)
  results.run$coverage.ri[j] <- count.ri/B
  results.run$coverage.agg[j] <- count.agg/B
  
}







save(results.run, file = "data_ri_for_rirs.rda")

pdf("studyplot_ri_for_rirs.pdf", width=7, height=7)



par(mar = c(5, 4, 4, 4) + 0.3)
plot(sd.rs.vec, results.run$width.ri, type = "l", col = "black", xlab = "Standard deviation random slope", 
     ylab = "Average width of the slope confidence interval",
     main = "RI model fitted to RIRS data", lwd = 2,
     ylim = c(0,0.6))
points(sd.rs.vec, results.run$width.agg, type = "l", col = "red", lwd = 2)

par(new = TRUE)

plot(sd.rs.vec, results.run$coverage.ri, type = "l", lty = 2, lwd = 2, axes = FALSE,
     bty = "n", xlab = "", ylab = "", col = "black", ylim = c(0.4, 1.02))
points(sd.rs.vec, results.run$coverage.agg, type = "l", lty = 2, lwd = 2, col = "red")

axis(side=4, at = pretty(range(results.run$coverage.ri)))
mtext("Confidence interval coverage frequency", side=4, line=3)

legend("topleft", legend = c("CI width RI model", "CI width Agg. Regr."), 
       col = c("black", "red"), lwd = 2, lty = c(1,1), cex=0.8)

legend("topright", legend = c("Coverage RI model", "Coverage Agg. Regr."), 
       col = c("black", "red"), lwd = 2, lty = c(2,2), cex=0.8)


dev.off()






