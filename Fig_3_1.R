# this script generates Fig. 3.1

library(lme4)
library(dplyr)


# structure of the data
n.ind <- 50

n.meas.min <- 3
n.meas.max <- 30
n.meas.vec <- n.meas.min:n.meas.max

# population intercept
a <- 1

# population slope
b <- 2

# error distribution
sd.error <- 0.5

# random intercept distribution
sd.ri <- 2

# number of simulations
B <- 300




set.seed(1)

# function to simulate data from an RI model
RI.simulate <- function(n.ind, n.meas, a, b, sd.error, sd.ri){
  n.total <- n.ind * n.meas
  
  #structure of the data set
  individual <- rep(1:n.ind, each = n.meas)
  
  # days of sleep deprivation 
  x <- runif(n=n.total, min = 0, max = 10)
  
  # errors
  error <- rnorm(n.total, mean = 0, sd = sd.error)
  
  # random intercept
  gamma <- rnorm(n.ind, mean = 0, sd = sd.ri)
  
  # linear RI model
  y <- a + b*x + gamma[individual] + error
  d <- data.frame(x,y,individual)
  
  return(d)
}





# single simulation run 
# generates data set, fits both models on same data set, and computes relevant estimates
combined.sim.run <- function(n.ind, n.meas, a, b, sd.error, sd.ri){
  
  data <- RI.simulate(n.ind, n.meas, a, b, sd.error, sd.ri)
  
  
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
results.run <- data.frame(width.ri = n.meas.vec, width.agg = n.meas.vec, 
                     coverage.ri = n.meas.vec, coverage.agg = n.meas.vec)







for(j in 1:length(n.meas.vec)){
  n.meas <- n.meas.vec[j]
  
  
  # counts how many times the conf.int covered the mean
  count.ri <- 0
  count.agg <- 0
  
  # store the widths
  width.ri <- numeric(B)
  width.agg <- numeric(B)
  

  
  for(i in 1:B){
    res.df <- combined.sim.run(n.ind, n.meas, a, b, sd.error, sd.ri)
    
    cover.ri <- res.df[1,1]
    count.ri <- count.ri + cover.ri
    
    cover.agg <- res.df[2,1]
    count.agg <- count.agg + cover.agg
    
    width.ri[i] <- res.df[1,3] - res.df[1,2]
    
    width.agg[i] <- res.df[2,3] - res.df[2,2]
    
    
    # progress update
    if(i%%100 == 0){print(i+n.meas)}
  }
  
  
  results.run$width.ri[j] <- mean(width.ri)
  results.run$width.agg[j] <- mean(width.agg)
  results.run$coverage.ri[j] <- count.ri/B
  results.run$coverage.agg[j] <- count.agg/B
  
}


save(results.run, file = "data_ri_vio1_study.rda")

pdf("studyplot_vio1.pdf", width=7, height=7)
plot(n.meas.vec, results.run$width.ri, type = "l", col = "black", xlab = "n_meas", 
     ylab = "Average width of the slope confidence interval",
     main = "Comparison of confidence interval widths for varying n_meas", lwd = 2)
points(n.meas.vec, results.run$width.agg, type = "l", col = "red", lwd = 2)
legend("topright", legend = c("Random Intercept Model", "Aggregated Regressions"), 
       col = c("black", "red"), lwd = 2, cex=0.8)

dev.off()