library(deSolve)

#--------------------------------------------------
# Example from deSolve Tutorial
# 1. A simple ODE: chaos in the atmosphere
#--------------------------------------------------

parameters <- c(a = -8/3, b=-10, c=28)
state <- c(X=1, Y=1, Z=1)

Lorenz <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    #rate of change
    dX <-a*X + Y*Z
    dY <-b * (Y-Z)
    dZ <- -X*Y + c*Y - Z
    #return tge rate of change
    list(c(dX,dY,dZ))
  }) #end with as.list
}

# run model for 100 days and give output at 0.01 daily intervals
# use seq() to create the time sequence
times <- seq(0,100,  by = 0.01)

#use ode() as default integration routine
#ode() input state variabel y, times, model function and parameters. 
#returns matrix with stat variables at requested output times

# ode = ordenairy differential equations
out <- ode(y = state, times = times, func = Lorenz, parms = parameters)
head(out)

#make nice plot
par(mar=c(5,5,2,5))
par(oma = c(0,0,3,0))
plot(out, xlab = "time", ylab = "-")
plot(out[, "X"], out[, "Z"], pch = ".")
mtext(outer = TRUE, side = 3, "Lorenz model", cex = 1.5)

#-------------------------------------------------------------------------
# 2. Solvers for initial value problems of ordinary differential equations
#-------------------------------------------------------------------------

# you can change the default integration method by using method

outc <- ode(state, times, Lorenz, parameters, method = "radau", atol = 1e-4, rtol=1e-4)

# using different integration routines and printing the time it takes
print(system.time(out1 <- rk4 (state, times, Lorenz, parameters)))
print(system.time(out2 <- lsode (state, times, Lorenz, parameters)))
print(system.time(out3 <- lsoda (state, times, Lorenz, parameters)))
print(system.time(out4 <- lsodes (state, times, Lorenz, parameters)))
print(system.time(out5 <- daspk (state, times, Lorenz, parameters)))
print(system.time(out6 <- vode (state, times, Lorenz, parameters)))

# to get an idea of all Runge-Kutta methods
help(rkMethod)

#In numerical analysis, the Runge–Kutta methods are a 
#family of implicit and explicit iterative methods, which includes the well-known 
#routine called the Euler Method, used in temporal discretization for the approximate 
#solutions of ordinary differential equations

#you also can define your own method

func <- function(t, x, parms){
  with(as.list(c(parms, x)),{
    dP <- a * P     - b * C * P
    dC <- b * P * C - c * C
    res <- c(dP, dC)
    list(res)
  })
}
# new method
rKnew <- rkMethod(ID = "midpoint", 
                  varstep = FALSE,
                  A       = c(0, 1/2),
                  b1      = c(0,1),
                  c       = c(0, 1/2),
                  stage   = 2,
                  Qerr    = 1)
out7 <-ode(y = c(P=2, C=1), 
           times = 0:100, 
           func, 
           parms = c(a=0.1, b=0.1, c=0.1), 
           method = rKnew)

head(out7)

#0.01 seconds step
times <- seq(0,40,  by = 0.01)
out <- ode(y = state, times = times, func = Lorenz, parms = parameters, method = "euler")
head(out)
# 1 second step
out <- ode(y = state, times = 0:40, func = Lorenz, parms = parameters, method = "euler")
head(out)
# prefered 1 second but smaller steps
out <- ode(y = state, times = 0:40, func = Lorenz, parms = parameters, method = "euler", hini=0.01)
head(out)


# running diagnostics
diagnostics(out1)
summary(out1)





