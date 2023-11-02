# this script generates Fig. 1.2


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




# generate Fig. 1.2 ######################

plot(data$x, data$y, col=(data$subject+1), pch = 16, xaxt='n', yaxt='n', main = "Synthetic data with indicated individuals", xlab = "", ylab = "")
legend(x=-2, y=8, legend=c("Individual 1", "Individual 2", "Individual 3"), pch = 16, col=c(2,3,4),cex=0.8,y.intersp=2)



# plot data set with different regression lines
plot(data$x, data$y, col=(data$subject+1), pch = 16, xaxt='n', yaxt='n', main = "Random intercept model fits", xlab = "", ylab = "")
abline(linmod, , lwd = 2, col = "black")
legend(x=-2, y=8, legend=c("Linear model", "Individual 1 (RI)", "Individual 2 (RI)", "Individual 3 (RI)"), col=c("black",2,3,4), lty=1, cex=0.8, lwd=2,y.intersp=2)


# fit random intercept model
RI.model <- lmer(y ~ x + (1 | subject), data = data)
summary(RI.model)


intercept.RI <- summary(RI.model)$coefficients[1,1]
slope.RI <- summary(RI.model)$coefficients[2,1]


for(i in 1:n.subject){
  randef <- ranef(RI.model)$subject[i,1]
  abline(a = (intercept.RI + randef), b = slope.RI, col = i+1, , lwd = 2)
}




