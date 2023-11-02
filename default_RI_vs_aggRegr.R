# this script is used for the "default" simulation in Section 2.1
# we compare the random intercept model with the aggregated regression approach


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

# number of simulations
B <- 3000



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




# counts how many times the conf.int covered the mean
count.ri <- 0
count.agg <- 0

# store the ci endpoints
ci.left.ri <- numeric(B)
ci.left.agg <- numeric(B)
ci.right.ri <- numeric(B)
ci.right.agg <- numeric(B)

# store the means
est.means.ri <- numeric(B)
est.means.agg <- numeric(B)

# store the widths
width.ri <- numeric(B)
width.agg <- numeric(B)

set.seed(1)

for(i in 1:B){
  res.df <- combined.sim.run(n.ind, n.meas, a, b, sd.error, sd.ri)
  
  est.means.ri[i] <- res.df[1,4]
  est.means.agg[i] <- res.df[2,4]
  
  cover.ri <- res.df[1,1]
  count.ri <- count.ri + cover.ri
  
  cover.agg <- res.df[2,1]
  count.agg <- count.agg + cover.agg
  
  width.ri[i] <- res.df[1,3] - res.df[1,2]
  
  width.agg[i] <- res.df[2,3] - res.df[2,2]
  
  ci.left.ri[i] <- res.df[1,2]
  ci.left.agg[i] <- res.df[2,2]
  ci.right.ri[i] <- res.df[1,3]
  ci.right.agg[i] <- res.df[2,3]
  
  # progress update
  if(i%%100 == 0){print(i)}
}

cat("The successful coverage frequency of ri is ", count.ri/B)
cat("The successful coverage frequency of agg is ", count.agg/B)

cat("The average width of a CI using ri is ", mean(width.ri))
cat("The average width of a CI using agg is ", mean(width.agg))


# determine number of necessary simulation runs (B)
# goal: 95% CI with max +/- 1% deviation from 0.95
binom.test(x=count.ri, n=B, p=0.95)$conf.int


t.test(x = width.ri, y = width.agg, paired = TRUE)
# difference is width.ri - width.agg


results.ri <- data.frame(left.ri.ci = ci.left.ri, right.ri.ci = ci.right.ri, left.agg.ci = ci.left.agg,
                         right.agg.ci = ci.right.agg, mean.ri = est.means.ri, mean.agg = est.means.agg)

save(results.ri, file = "data_ri.rda")






# here we generate a plot of the first 50 confidence intervals (Fig. 2.1)


n.plot <- min(B, 50)

pdf("ci_ri_plot.pdf", width=7, height=7)

plot(est.means.ri[1:n.plot], (1:n.plot+0.15), pch = 16, col = "chartreuse3",
     xlim = c(1.95,2.05),
     ylim = c(0.5,n.plot+0.5),
     xlab = "Slope estimate and confidence interval width",
     ylab = "Index of the confidence interval",
     main = paste(n.plot, "Confidence intervals")
)

points(est.means.agg[1:n.plot], (1:n.plot)-0.15, pch = 16, col = "blue1")

# plot ri ci
for (i in 1:n.plot){
  if(between(b, ci.left.ri[i], ci.right.ri[i])){
    segments(ci.left.ri[i], i+0.15, ci.right.ri[i], i+0.15, lwd = 2, col = "chartreuse3") #plot CI's that contain the mean in black
  } else {
    segments(ci.left.ri[i], i+0.15, ci.right.ri[i], i+0.15, col = "red", lwd = 2) #plot CI's that don't contain the mean in red
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



legend("topright", legend = c("CI RI", "CI Agg. OLS", "True slope", "No coverage"), 
       col = c("chartreuse3", "blue1", "darkmagenta", "red"), lwd = 2, cex=0.8)

dev.off()