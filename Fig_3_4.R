# this script generates Fig. 3.4 (both plots)


library(lme4)
library(dplyr)


# structure of the data
n.ind.min <- 10
n.ind.max <- 70
n.ind.vec <- seq(from=n.ind.min, to=n.ind.max, by=5)

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
sd.rs <- 1

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
  
  
  
  
  # RIRS model ##########################################################
  RIRS.model <- lmer(y ~ x + (1 + x | individual), data = data)
  
  conf.int.rirs <- suppressMessages(confint(RIRS.model, level = 0.95))
  left.rirs <- conf.int.rirs[6,1]
  right.rirs <- conf.int.rirs[6,2]
  
  mean.rirs <- summary(RIRS.model)$coeff[2,1]
  #######################################################################
  
  c.rirs <- 0
  if((left.rirs<=b) && (right.rirs>=b)){c.rirs <- 1}
  
  
  
  
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
  
  
  
  output <- data.frame(c(c.rirs, c.agg), c(left.rirs, left.agg), c(right.rirs, right.agg), 
                       c(mean.rirs, mean.agg), row.names = c("rirs", "agg"))
  return(output)
}






#store the results
results.run <- data.frame(width.rirs = n.ind.vec, width.agg = n.ind.vec, 
                          coverage.rirs = n.ind.vec, coverage.agg = n.ind.vec)





for(j in 1:length(n.ind.vec)){
  n.ind <- n.ind.vec[j]
  
  
  # counts how many times the conf.int covered the mean
  count.rirs <- 0
  count.agg <- 0
  
  # store the widths
  width.rirs <- numeric(B)
  width.agg <- numeric(B)
  
  
  
  for(i in 1:B){
    res.df <- combined.sim.run(n.ind, n.meas, a, b, sd.error, sd.ri, sd.rs)
    
    cover.rirs <- res.df[1,1]
    count.rirs <- count.rirs + cover.rirs
    
    cover.agg <- res.df[2,1]
    count.agg <- count.agg + cover.agg
    
    width.rirs[i] <- res.df[1,3] - res.df[1,2]
    
    width.agg[i] <- res.df[2,3] - res.df[2,2]
    
    
    
    # progress update
    if(i%%100 == 0){print(i+n.ind)}
  }
  
  
  results.run$width.rirs[j] <- mean(width.rirs)
  results.run$width.agg[j] <- mean(width.agg)
  results.run$coverage.rirs[j] <- count.rirs/B
  results.run$coverage.agg[j] <- count.agg/B
  
}



save(results.run, file = "data_rirs_vio2_study.rda")

pdf("studyplot_vio2_rirs.pdf", width=7, height=7)


plot(n.ind.vec, results.run$width.rirs, type = "l", col = "black", xlab = "n_ind", 
     ylab = "Average width of the slope confidence interval",
     main = "Comparison of confidence interval widths for varying n_ind", lwd = 2)
points(n.ind.vec, results.run$width.agg, type = "l", col = "red", lwd = 2)
legend("topright", legend = c("Random Coefficient Model", "Aggregated Regressions"), 
       col = c("black", "red"), lwd = 2, cex=0.8)


dev.off()







#plot with log-log axes
pdf("studyplot_vio2_rirs_loglog.pdf", width=7, height=7)

plot(n.ind.vec, results.run$width.rirs, type = "l", col = "black", xlab = "n_ind", 
     ylab = "Average width of the slope confidence interval",
     main = "Comparison of confidence interval widths for varying n_ind", lwd = 2, log = "xy")
points(n.ind.vec, results.run$width.agg, type = "l", col = "red", lwd = 2)
legend("topright", legend = c("Random Coefficient Model", "Aggregated Regressions"), 
       col = c("black", "red"), lwd = 2, cex=0.8)

dev.off()


