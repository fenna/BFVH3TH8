####################################
#  add 10 - 10% of volume example
####################################

library(deSolve)

# define volume to add and percentage to decrease
parameters <- c(addVolume = 10, pV = 0.1) 

# define function 
volume <- function(t,y,parms){
  with(as.list(c(parms)),{
         dY <- addVolume - pV * (y+addVolume)
         return(list(c(dY)))
       }
       )
}

#initial state
state <- c(Volume = 0)

#define time sequence you want to run the model
times <- seq(0, 100,  by = 1)

# run simulation using continuous approach
out  <- ode(times = times, y = state,   parms = parameters, func = volume, method = "euler")
head(out)
# plot results
plot(out)
