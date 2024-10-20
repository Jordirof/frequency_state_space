# Author: Jordi Rof
# Purpose: 
#   This script shows how to obtain system-wide frequency readings from local frequency readings.
#   The modelling technique of is state-space, where system-level frequency is considered unobservable.
#   The data used is synthetic, constructed from real Frequency data downloaded from Elexon.
#   See: https://bmrs.elexon.co.uk/rolling-system-frequency
# Dependencies: MARSS, MASS, ggplot2, rstudioapi

# Set wd
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Import data (downloaded from Elexon, see above)
freq <- read.csv("data/RollingSystemFrequency-2024-10-12T13_00_00.000Z-2024-10-19T13_00_00.000Z.csv")
freq <- freq[1:1000,] # 1000 data points are selected

# Generate data - disturbance same var, zero covar, different means
## create the variance covariance matrix
sigma <- diag(5)/10
## create the mean vector
mu<-c(0, -10, 5, -5, 10)/100
## Generate error (set seed to ensure it can be replicated)
## Add the frequency to it
set.seed(123); data <- MASS::mvrnorm(n=nrow(freq), mu=mu, Sigma=sigma) + as.vector(freq[,2])
## Plot and check properties
cor(data) # Observations are not strongly correlated
var(data) # Var-covar very close to the disturbance generated
plot(data.frame(data)) # Very low correlation confirmed visually

# Model 1 
## Multivariate model estimating only one state
mod.list1 <- list(
  # U - Drift is set to zero, frequency does not have a trend
  U = matrix(0),
  # x0 - Initial state set to frequency target
  x0 = matrix(50),
  # B - We have no ex ante information of the AR coefficient value, we estimate it
  B = matrix("b"), 
  # Q - This is the weakest assumption, but we need to fix either this of 'B' to have a system with solution
  Q = matrix(1/20),
  # Z - This is ensuring that we account for the relationship between state and space in each series
  Z = matrix(c(1,1,1,1,1)),
  # A - Drift in Yt. Estimated. It could also be set to zero
  A = "unconstrained",
  # R - Variance-covariance matrix of residuals. If there were 2 y series, matrix(list("r1",0,0,"r2"),2,2).
  R = "diagonal and unequal",
  # The tinitx element tells MARSS whether the initial state for x is at t=1 (tinitx=1) or t=0 (tinitx=0)
  tinitx = 0
)

# Model 2
## Multivariate model estimating only one state
mod.list2 <- list(
  # U - Drift is set to zero, frequency does not have a trend
  U = matrix(0),
  # x0 - Initial state set to frequency target
  x0 = matrix(50),
  # B - We have no ex ante information of the AR coefficient value, we estimate it
  B = matrix("b"), 
  # Q - This is the weakest assumption, but we need to fix either this of 'B' to have a system with solution
  Q = matrix(1/20),
  # Z - This is ensuring that we account for the relationship between state and space in each series
  Z = matrix(c(1,1,1,1,1)),
  # A - Drift in Yt. Estimated. It could also be set to zero
  A = "unconstrained",
  # R - Variance-covariance matrix of residuals. If there were 2 y series, matrix(list("r",0,0,"r"),2,2).
  R = "diagonal and equal",
  # The tinitx element tells MARSS whether the initial state for x is at t=1 (tinitx=1) or t=0 (tinitx=0)
  tinitx = 0
)

# Model 3
## Multivariate model estimating only one state
mod.list3 <- list(
  # U - Drift is set to zero, frequency does not have a trend
  U = matrix(0),
  # x0 - Initial state set to frequency target
  x0 = matrix(50),
  # B - We have no ex ante information of the AR coefficient value, we estimate it
  B = matrix(1), 
  # Q - This is the weakest assumption, but we need to fix either this of 'B' to have a system with solution
  Q = matrix("q"),
  # Z - This is ensuring that we account for the relationship between state and space in each series
  Z = matrix(c(1,1,1,1,1)),
  # A - Drift in Yt. Estimated. It could also be set to zero
  A = "unconstrained",
  # R - Variance-covariance matrix of residuals. If there were 2 y series, matrix(list("r",0,0,"r"),2,2).
  R = "diagonal and equal",
  # The tinitx element tells MARSS whether the initial state for x is at t=1 (tinitx=1) or t=0 (tinitx=0)
  tinitx = 0
)

# Fit models
## Fit the model 1
fit1 <- MARSS::MARSS(t(data),model = mod.list1)
#ggplot2::autoplot(fit1) # Plots isf and diagnostics
#ggplot2::autoplot(MARSS::forecast(fit1,h=100)) # plots fcast

## Fit the model 2
fit2 <- MARSS::MARSS(t(data),model = mod.list2)
#ggplot2::autoplot(fit2) # Plots isf and diagnostics
#ggplot2::autoplot(MARSS::forecast(fit2,h=100)) # plots fcast

## Fit the model 3
fit3 <- MARSS::MARSS(t(data),model = mod.list3)
#ggplot2::autoplot(fit2) # Plots isf and diagnostics
#ggplot2::autoplot(MARSS::forecast(fit2,h=100)) # plots fcast

# compare predictions
freq$state1 <- as.vector(fit1$states)
freq$state2 <- as.vector(fit2$states)
freq$state3 <- as.vector(fit3$states)
rownames(freq) <- freq$MeasurementTime
freq$MeasurementTime <- NULL
## Correlation
cor(freq)
## RMSE
colSums(((freq[,2:4] - freq[,1])^(2)/nrow(freq))^(1/2))

# Compare selection metrics
aic <- c(fit1$AIC, fit2$AIC, fit3$AIC)
loglik <- c(fit1$logLik, fit2$logLik, fit3$logLik)

# Compare some of the parameters to their true values
## Compare drifts
### Collect data
parameter <- rep(paste0("A",1:5),4)
model <- c(rep("Actual",length(mu)), rep("M1",length(mu)), rep("M2",length(mu)),rep("M3",length(mu))) 
value <- c(mu, as.vector(fit1$par$A), as.vector(fit2$par$A), as.vector(fit3$par$A))
data_c1 <- data.frame(value,parameter,model)
### Create chart
ggplot2::ggplot(data_c1, ggplot2::aes(fill=model, y=value, x=parameter)) + 
  ggplot2::geom_bar(position="dodge", stat="identity")
## Compare drifts
### Collect data
parameter <- rep(paste0("R",1:5),4)
value <- c(diag(sigma), as.vector(fit1$par$R), rep(as.vector(fit2$par$R),5), rep(as.vector(fit3$par$R),5))
data_c2 <- data.frame(value,parameter,model)
### Create chart
ggplot2::ggplot(data_c2, ggplot2::aes(fill=model, y=value, x=parameter)) + 
  ggplot2::geom_bar(position="dodge", stat="identity")

# Plot the different models
## Correlation plot
plot(freq[,1:4])
## Time series comparison
freq$Period <- c(1:nrow(freq))
freqrs <- reshape2::melt(freq[c(1,4,ncol(freq))], id="Period")
ggplot2::ggplot(freqrs, ggplot2::aes(x = Period, y = value, group=variable)) +
  ggplot2::geom_line() +
  ggplot2::theme(legend.position = "bottom")

