##################################################
## Project: Week 4: Detailed analysis of the model
## Script purpose: Model when drug treatment is stopped
## Date: 2020-05-25
## Author: "James Gray, Naomi Hindriks"
##################################################

new.state <- state
new.state["MPL"] <- parameters["MPL"]

new.parms <- parameters
new.parms <- new.parms[-which(names(new.parms) == "MPL")]

new.dynamic.model.equations <- function(t,state,parms){
  with(as.list(c(state, parms)),{
    
    # Calculate change in mRNA
    d.mRNA <- synthesized.mRNA * (1 - DRN/(ic.50 + DRN)) - decayed.mRNA.fraction * (mRNA)
    
    # Calculate change in free receptors in cytoplasm
    d.R <- (synthesized.receptor.frac * mRNA 
            + recycle.frac * decay.DRN.frac * DRN 
            - complex.formation.factor * MPL * R 
            - receptor.decay.frac * R)
    
    # Calculate change in receptor coplex in cytoplasm
    d.DR <- complex.formation.factor * MPL * R - DR.transport.frac * DR
    
    # Calculate change in receptor complex in nucleus
    d.DRN <- DR.transport.frac * DR - decay.DRN.frac * DRN
    
    d.MPL <- 0
    
    # Return calculated values
    return(list(c(d.mRNA, d.R, d.DR, d.DRN, d.MPL)))
  }
  )
}

rootfun <- function(time, state, pars){
  dstate <- unlist(new.dynamic.model.equations(time, state, pars))
  rep(sum(abs(dstate)) - 0.1, times = 2)
}

eventfun <- function(time, state, pars) {
  state["MPL"] <- 0
  state
}

out  <- ode(
  # Run the simulation for this amount of times
  times = seq(0, 200,  by = 0.1), 
  
  # Set initial values to run simutation with
  y = new.state,
  
  # Set parameters to run simulation with
  parms = new.parms,
  
  # Use the dynamic.model.equations function  to calculate the change in the 
  # initial values (mRNA, R, DR, DRN) for every time step
  func = new.dynamic.model.equations, 
  
  rootfun = rootfun,
  
  events = list(func = eventfun, root = TRUE, maxroot = 2)
)

max.time <- ceiling(attributes(out)$troot[2])

par(
  mfrow = c(3, 2),
  mar = c(5, 6, 4, 2) + 0.1,
  oma = c(0, 0, 3, 0)
)

plot(
  out,
  xlim = c(0, max.time),
  lwd = 2,
  mfrow = NULL,
  main = titles,
  xlab = "tijd (uur)",
  ylab = ylabs,
  col = color.palet[8]
)

mtext(
  outer = TRUE, 
  side = 3, 
  "Concentraties nadat medicatie stopt", 
  cex = 1.35
)
