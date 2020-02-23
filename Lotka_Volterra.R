###################################################
# 
# according deSolve paper
# predator-prey model Lotka-Volterra
#
###################################################

###################################################
### preliminaries
###################################################
par(mfrow=c(2,2))
library(deSolve)


###################################################
### 0-D Prey-Consumer model
###################################################

# variabeles of interest:
# C = consumer concentration
# P = Prey concentration

# parameter values
pars <- c(rI = 0.2,  # ingestion rate of consumer
          rG = 1.0,  # growth rate
          rM = 0.2,  # consumers mortality rate
          AE = 0.5,  # assimilation efficiency
          K = 10)    # carrying capacity

# model function
LVmod <- function(time, state, pars){
  with(as.list(c(state, pars)), {
    IngestC <- rI * P * C
    GrowthP <- rG * P * (1 - P/K)
    MortC   <- rM * C
    
    dP <- GrowthP - IngestC
    dC <- IngestC * AE - MortC
    
    return(list(c(dP, dC)))
  })
}

# stop criteria: stop when change is less then 1e-4 (10 ^ -4)
rootfun <- function(time, state, pars){
  dstate <- unlist(LVmod(time, state, pars))
  sum(abs(dstate)) - 1e-4
}

# initial state and timeframe
yini <- c(P = 1, C = 2)
times <- seq(0, 200, by = 1)

###################################################
# simulate without stopcriteria
###################################################
print(system.time(
  out <- ode(times = times, y = yini, parms = pars, func = LVmod)
))

# results
head(out, n = 3)
tail(out, n = 3)

# plot
matplot(out[,"time"], out[,2:3], type = "l", xlab = "time", ylab = "Conc",
        main = "Lotka-Volterra", lwd = 2)
legend("topright", c("prey", "predator"), col = 1:2, lty = 1:2)

###################################################
# simulate with stopcriteria
###################################################
print(system.time(
  out <- ode(times = times, y = yini, parms = pars, func = LVmod, rootfun = rootfun)
))

# results
head(out, n = 3)
tail(out, n = 3)


# plot
matplot(out[,"time"], out[,2:3], type = "l", xlab = "time", ylab = "Conc",
        main = "Lotka-Volterra with root", lwd = 2)

###################################################
# Consumer and prey dispersing on a 1-D grid
# dispersion in x direction
###################################################

# extra parameters
R <- 20     # total length of surface, m
N <- 1000   # devided by 1000 boxes
dx <- R/N   # box size in x direction
Da <- 0.05  # m2/d, dispersion coefficient

# model with flux correction
# by imposing P[1] and P[N] at the upper and lower boundaries, 
# we effectively impose a zero-gradient (or a zero-flux) boundary condition

LVmodDis <- function(time, state, pars, N, Da, dx){
  with(as.list(c(state, pars)), {
    
    P <- state[1:N]
    C <- state[-(1:N)]
    
    ## Dispersive fluxes; zero-gradient boundaries
    FluxP <- -Da * diff(c(P[1], P, P[N]))/dx
    FluxC <- -Da * diff(c(C[1], C, C[N]))/dx
    
    ## Biology: Lotka-Volterra dynamics
    IngestC <- rI * P * C
    GrowthP <- rG * P * (1 - P/K)
    MortC   <- rM * C
    
    ## Rate of change = -Flux gradient + Biology
    dP <- - diff(FluxP)/dx + GrowthP - IngestC
    dC <- - diff(FluxC)/dx + IngestC * AE - MortC
    
    return(list(c(dP, dC)))
  })
}

# initial value and time sequence
yini <- rep(0, 2*N)
# define start value in each box
yini[500:501] <- yini[1500:1501] <- 10
times <-seq(0, 200, by = 1)

print(system.time(out1D <- ode.1D(y = yini, times = times, func = LVmodDis,
                                  parms = pars, nspec = 2, N = N, dx = dx, Da = Da)))

# The matrix out has in its first column the time sequence, 
# followed by 1000 columns with prey concentrations, one for each box, 
# followed by 1000 columns with consumer concentrations. 
head(out)

Z   <- out1D[,2:(N + 1)]      
filled.contour(x = times, z = Z, y = seq(0, R, length=N),
               color = gray.colors, xlab = "Time, days", ylab= "Distance, m",
               main = "Prey density")



###################################################
# Consumer and prey dispersing on a 2-D grid
# Dispersion in x and y direction
###################################################

# extra parameters
R <- 20     # total lenght of surface, m
N <- 50     # devided by 1000 boxes
dx <- R/N   # box size x direction
dy <- R/N   # box size y direction
Da <- 0.05  # m2/d, dispersion coefficient
NN <- N * N 

# model
LVmod2D <- function (time, state, parms, N, Da, dx, dy) {
  P<- matrix(nr = N, nc = N, state[1:NN])
  C<- matrix(nr = N, nc = N, state[-(1:NN)])
  with (as.list(parms), {
    dP    <- rG * P *(1 - P/K) - rI * P *C
    dC    <- rI * P * C * AE - rM * C
    
    zero  <- numeric(N)
    ## Fluxes in x-direction; zero fluxes near boundaries
    FluxP <- rbind(zero, -Da * (P[-1,] - P[-N,])/dx, zero)
    FluxC <- rbind(zero, -Da * (C[-1,] - C[-N,])/dx, zero)
    
    dP    <- dP - (FluxP[-1,] - FluxP[-(N+1),])/dx
    dC    <- dC - (FluxC[-1,] - FluxC[-(N+1),])/dx
    
    ## Fluxes in y-direction
    FluxP <- cbind(zero, -Da * (P[,-1] - P[,-N])/dy, zero)
    FluxC <- cbind(zero, -Da * (C[,-1] - C[,-N])/dy, zero)
    
    dP    <- dP - (FluxP[,-1] - FluxP[,-(N+1)])/dy
    dC    <- dC - (FluxC[,-1] - FluxC[,-(N+1)])/dy
    
    return(list(c(as.vector(dP), as.vector(dC))))
  })
}

# yini now matrix
yini <- rep(0, 2 * N * N)
cc <- c((NN/2):(NN/2 + 1) + N/2, (NN/2):(NN/2 + 1) - N/2)
yini[cc] <- yini[NN + cc] <- 10
times  <- seq(0, 200, by = 1)


# run simulation
print(system.time(
  out <- ode.2D(y = yini, times = times, func = LVmod2D,
                parms = pars, dimens = c(N, N), N = N,
                dx = dx, dy = dy, Da = Da, ynames = FALSE, lrw = 440000)))

# plot 2-D image
P   <- out[,2:(N + 1)]
par(mfrow=c(2,2))
par(oma=c(0,0,2,0))
xx <- seq(0, R, dx)
yy <- seq(0, R, dy)
Col <- gray.colors
image(x=xx, y=yy, z=matrix(nr=N,nc=N,out[1,-1]),  zlim=c(0,10), col=(Col(100)),
      main="initial", xlab="x", ylab="y")
image(x=xx, y=yy, z=matrix(nr=N,nc=N,out[21,-1]), zlim=c(0,10), col=(Col(100)),
      main="20 days", xlab="x", ylab="y")
image(x=xx, y=yy, z=matrix(nr=N,nc=N,out[31,-1]), zlim=c(0,10), col=(Col(100)),
      main="30 days", xlab="x", ylab="y")
image(x=xx, y=yy, z=matrix(nr=N,nc=N,out[41,-1]), zlim=c(0,10), col=(Col(100)),
      main="40 days", xlab="x", ylab="y")
mtext(side=3, outer=TRUE, cex=1.25, 
      "Lotka-Volterra Prey concentration on 2-D grid")
