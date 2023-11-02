# model violation 3: n.meas is different for each individual
# this script does the simulation and generates the corresponding plot in Fig. A.3


library(lme4)
library(dplyr)


# structure of the data
n.ind <- 50

set.seed(1)
n.meas <- sample(x = c(5,10,25,50,100,200), size = n.ind, replace = TRUE)

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
B <- 3000


# function to simulate data from a RIRS model
RIRS.simulate <- function(n.ind, n.meas, a, b, sd.error, sd.ri, sd.rs){
  n.total <- sum(n.meas)
  
  #structure of the data set
  individual <- rep(1:n.ind, times = n.meas)
  
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






# counts how many times the conf.int covered the mean
count.rirs <- 0
count.agg <- 0

# store the ci endpoints
ci.left.rirs <- numeric(B)
ci.left.agg <- numeric(B)
ci.right.rirs <- numeric(B)
ci.right.agg <- numeric(B)

# store the means
est.means.rirs <- numeric(B)
est.means.agg <- numeric(B)

# store the widths
width.rirs <- numeric(B)
width.agg <- numeric(B)

set.seed(3)

for(i in 1:B){
  res.df <- combined.sim.run(n.ind, n.meas, a, b, sd.error, sd.ri, sd.rs)
  
  est.means.rirs[i] <- res.df[1,4]
  est.means.agg[i] <- res.df[2,4]
  
  cover.rirs <- res.df[1,1]
  count.rirs <- count.rirs + cover.rirs
  
  cover.agg <- res.df[2,1]
  count.agg <- count.agg + cover.agg
  
  width.rirs[i] <- res.df[1,3] - res.df[1,2]
  
  width.agg[i] <- res.df[2,3] - res.df[2,2]
  
  ci.left.rirs[i] <- res.df[1,2]
  ci.left.agg[i] <- res.df[2,2]
  ci.right.rirs[i] <- res.df[1,3]
  ci.right.agg[i] <- res.df[2,3]
  
  # progress update
  if(i%%10 == 0){print(i)}
}

cat("The successful coverage frequency of rirs is ", count.rirs/B)
cat("The successful coverage frequency of agg is ", count.agg/B)

cat("The average width of a CI using rirs is ", mean(width.rirs))
cat("The average width of a CI using agg is ", mean(width.agg))


# determine number of necessary simulation runs (B)
# goal: 95% CI with max +/- 1% deviation from 0.95
binom.test(x=count.rirs, n=B, p=0.95)$conf.int


alpha.bonferroni <- 0.05/12
t.test(x = width.rirs, y = width.agg, paired = TRUE, conf.level = 1-alpha.bonferroni)
# difference is width.rirs - width.agg




results.rirs.vio3 <- data.frame(left.rirs.ci = ci.left.rirs, right.rirs.ci = ci.right.rirs, left.agg.ci = ci.left.agg,
                                right.agg.ci = ci.right.agg, mean.rirs = est.means.rirs, mean.agg = est.means.agg)

save(results.rirs.vio3, file = "data_rirs_vio3.rda")








n.plot <- min(B, 50)

pdf("ci_rirs_plot_vio3.pdf", width=7, height=7)

plot(est.means.rirs[1:n.plot], (1:n.plot+0.15), pch = 16, col = "chartreuse3",
     xlim = c(1.4,2.6),
     ylim = c(0.5,n.plot+0.5),
     xlab = "Slope estimate and confidence interval width",
     ylab = "Index of the confidence interval",
     main = paste(n.plot, "confidence intervals (violation 3)")
)

points(est.means.agg[1:n.plot], (1:n.plot)-0.15, pch = 16, col = "blue1")

# plot rirs ci
for (i in 1:n.plot){
  if(between(b, ci.left.rirs[i], ci.right.rirs[i])){
    segments(ci.left.rirs[i], i+0.15, ci.right.rirs[i], i+0.15, lwd = 2, col = "chartreuse3") #plot CI's that contain the mean in black
  } else {
    segments(ci.left.rirs[i], i+0.15, ci.right.rirs[i], i+0.15, col = "red", lwd = 2) #plot CI's that don't contain the mean in red
  }
}


for (i in 1:n.plot){
  if(between(b, ci.left.agg[i], ci.right.agg[i])){
    segments(ci.left.agg[i], i-0.15, ci.right.agg[i], i-0.15, lwd = 2, col = "blue1") #plot CI's that contain the mean in black
  } else {
    segments(ci.left.agg[i], i-0.15, ci.right.agg[i], i-0.15, col = "red", lwd = 2) #plot CI's that don't contain the mean in red
  }
}

abline(v=b, col = "darkmagenta", lwd = 2) #plot a vertical line at the population mean


legend("topright", legend = c("CI RIRS", "CI Agg. Regr.", "True slope", "No coverage"), 
       col = c("chartreuse3", "blue1", "darkmagenta", "red"), lwd = 2, cex=0.8)

dev.off()
