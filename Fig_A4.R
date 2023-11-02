# This code generates Fig. A.4

library(dplyr)

set.seed(1)

# structure of the data
n.ind <- 4
n.meas <- 40

# population intercept
a <- 1

# population slope
b <- 2

# error distribution
sd.error <- 0.5

# random intercept distribution
sd.ri <- 2

# random slope distribution
sd.rs <- 0.2





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

# generate synthetic data set
dat <- RIRS.simulate(n.ind, n.meas, a, b, sd.error, sd.ri, sd.rs)

pdf("almost_RI.pdf", width=7, height=7)
plot(dat$x, dat$y, col = dat$individual, pch = 16, xlab = "x", ylab = "y")
dev.off()



