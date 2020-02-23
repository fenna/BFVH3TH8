#------------------------------------------------------------------
# Example from deSolve tutorial
# 3. Partial differential equations
#------------------------------------------------------------------


# modelling of Aphids (pest insect) slowly diffuse and grow on a row of plants
# model equation
Aphid <- function(t, APHIDS, parameters) {
  deltax     <- c(0.5, rep(1, numboxes -1), 0.5)
  Flux       <- -D * diff(c(0, APHIDS, 0)) / deltax  # Flux = -D*(delta_N/delta_x)
  dAPHIDS    <- - diff(Flux) / delx + APHIDS * r
  
  # return value 
  list(dAPHIDS)
}

# model parameters
D         <- 0.3  # m2/day   diffusion rate
r         <- 0.01 # /day    net growth rate
delx      <- 1    # m       thickness of a box
numboxes  <- 60
Distance  <- seq(from = 0.5, by = delx, length.out = numboxes)

# initiation conditions
APHIDS        <- rep(0, times = numboxes)
APHIDS[30:31] <- 1
state         <- c(APHIDS = APHIDS)

# The model runs for 200 days 
times = seq(0, 200, by = 1)

# Run ode.D1 to integrate 1-dimensional problem compromizing one or many spieces
print(system.time(out <- ode.1D(state, times, Aphid, parms = 0, nspec = 1, names = "Aphid")))

# print time followed by densities of Aphids (1:5)
head(out[,1:5])

# plot the solution for the 1-dimensional aphid model 
image( out, method = "filled.contour", grid = Distance, 
       xlab = "time, days",
       ylab = "Distance on plant, m",
       main = "Aphid density on a row of plants")

# plot dsolution of the Aphid model - plotted with matplot.1D

data <- cbind(dist = c(0, 10, 20, 30, 40, 50, 60),
              Aphid = c(0, 0.1, 0.25, 0.5, 0.25, 0.1, 0))

par(mfrow = c(1, 2))
matplot.1D(out, grid = Distance, type = "l", mfrow = NULL,
           subset = time %in% seq(0, 200, by = 10),
           obs = data, obspar = list(pch = 18, cex = 2, col = "red"))

par(mfrow = c(1, 2))
plot.1D(out, grid = Distance, type = "l", mfrow = NULL,
        subset = time == 100,
        obs = data, obspar = list(pch = 18, cex = 2, col = "red"))
