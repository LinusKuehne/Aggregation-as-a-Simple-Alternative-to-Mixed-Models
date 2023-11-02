# this script generates Fig. 1.1 and Fig. 1.2


library(lme4)

# function to simulate data from a RI model
RI.simulate <- function(n.subject, n.meas.per.subject, a, b, sd.error, sd.ri){
  n.total <- n.subject * n.meas.per.subject
  
  #structure of the data set
  subject <- rep(1:n.subject, each = n.meas.per.subject)
  
  # days of sleep deprivation 
  x = rnorm(n.total)
  
  # errors
  error <- rnorm(n.total, mean = 0, sd = sd.error)
  
  # random intercept
  gamma <- rnorm(n.subject, sd = sd.ri)
  
  # linear RI model
  y <- a + b*x + gamma[subject] + error
  d <-  data.frame(x,y,subject)
  
  return(d)
}

# structure of the data
n.subject <- 3
n.meas.per.subject <- 40


# population intercept
a <- 1

# population slope
b <- 2

# error distribution
sd.error <- 1.5

# random intercept distribution
sd.ri <- 2

set.seed(6)
data <- RI.simulate(n.subject, n.meas.per.subject, a, b, sd.error, sd.ri)


par(mfrow = c(1,2))



# generate Fig. 1.1 ######################

plot(data$x, data$y, pch = 16, xaxt='n', yaxt='n', col = 1, main = "Synthetic data", xlab = "", ylab = "")
plot(data$x, data$y, pch = 16, xaxt='n', yaxt='n', col = 1, 
     main = "Synthetic data with regression line", xlab = "", ylab = "")

linmod <- lm(y ~ x, data = data)
abline(linmod, , lwd = 2, col = "black")
legend(x=-2, y=8, legend=c("Linear model"), col="black", lty=1:2, cex=0.8, lwd=2,y.intersp=2)
