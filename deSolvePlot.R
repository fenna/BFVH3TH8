##############################################
# 
# Script from deSolve tutorial section 10.1/10.3
# Author tutorial: 
#     Karline Soeteart, Thomas Petzoldt, 
#     R. Woodrwo Setzer
# published 2010
# 
# This script demonstrates several plotting 
# options in combination with deSolve modelling
# 
##############################################

library(deSolve)

##############################################
# 10.1 plotting multiple scenario's
##############################################

# model
combustion <- function(t,y, parms){
  list(y^2 * (1-y))
}

# parameters
yini <- 0.01
times <- 0:200

# solve with different initial conditions for y
out  <- ode(times = times, y = yini,   parms = 0, func = combustion)
out2 <- ode(times = times, y = yini*2, parms = 0, func = combustion)
out3 <- ode(times = times, y = yini*3, parms = 0, func = combustion)
out4 <- ode(times = times, y = yini*4, parms = 0, func = combustion)

#plot the result all in one
plot(out, out2, out3, out4, main="combustion")
legend("bottomright", lty = 1:4, col = 1:4, legend = 1:4, title = "yini*i")



##############################################
# 10.2 plotting output with observations
##############################################

# we make use of the ccl4data set
head(ccl4data)
# and the ccl4model
help("ccl4model")

# we select subset of data animal "A" from ccl4data 
obs <- subset (ccl4data, animal == "A", c(time, ChamberConc))
names(obs) <- c("time", "CP")
head(obs)

# assign parameter values
parms <- c(0.182, 4.0, 4.0, 0.08, 0.04, 0.74, 0.05, 0.15, 0.32, 16.17, 
           281.48, 13.3, 16.17, 5.487, 153.8, 0.04321671,
           0.40272550, 951.46, 0.02, 1.0, 3.80000000)
#assign y initial
yini <- c(AI = 21, AAM = 0, AT = 0, AF = 0, AL = 0, CLT = 0, AM = 0)
times <- seq(0, 6, by = 0.05)

#call model
out <- ccl4model(times = times, y = yini, parms = parms)

#tweak parameters a bit (scenario 2 and 3)
par2 <- parms
par2[1] <- 0.1
par3 <- parms
par3[1] <- 0.05

#and run model again
out2 <- ccl4model(times = times, y = yini, parms = par2)
out3 <- ccl4model(times = times, y = yini, parms = par3)

# plot AI, MASS and CP for all scenario's
plot(out, out2, out3, which = c("AI", "MASS", "CP"),
     col = c("black", "red", "green"), lwd = 2,
     obs = obs, obspar = list(pch = 18, col = "blue", cex = 1.2))
legend("topright", lty = c(1,2,3,NA), pch= c(NA, NA, NA, 10), cex = 0.4,
       col = c("black", "red", "green", "blue"), lwd = 2,
       legend = c("par1", "par2", "par3", "obs"))

# get data for a specific measurement
obs2 <- data.frame(time = 6, MASS = 12)

# plot variable in common with observations
plot(out, out2, out3, lwd = 2,
     obs = list(obs, obs2),
     obspar = list(pch = c(16, 18), col = c("blue", "black"),
                   cex = c(1.2, 2))
    )


##############################################
# 10.3 plotting summary histograms
##############################################

# plotting histograms of all output variables
hist(out, col = grey(seq(0, 1, by = 0.1)), mfrow = c(3, 4))





