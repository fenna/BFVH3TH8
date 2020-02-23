###########################################################
# Differential equations for fift-generation model of corticosteroid receptor dynamics
# ---------------------------------------------------------
# Course: Thema 8
# -----------------------

# Remove all objects
rm(list=ls(all=TRUE))

# Load deSolve library
library(deSolve)

###################
# Specify functions
# -----------------

# User specified function: 

# four-equations dynamica system
# --------------------------------------------------------------------
fifthGenerationModel <- function(t, x, parms) 
{
  Rm <- x[1]           # mRNA for the receptor
  R <- x[2]           # receptor density
  DR <- x[3]           # drug-receptor comples
  DRN <- x[4]         # nuclear activated drug-receptor complex
  
  with(as.list(parms),
  {
     dRm <- ksyn_rm*(1-DRN/(IC50 + DRN)) - kdgr_rm*Rm
     dR <- ksyn_r*Rm + rf*kre*DRN - kon*D*R - kdgr_r*R
     dDR <- kon*D*R - kt*DR
     dDRN <- kt*DR - kre*DRN
 
     # store results in res
     res<-c(dRm,dR,dDR,dDRN)
     list(res)
  })
}


# User specified function: 
# Initialize the system
# could also be a vector
# ----------
systemInit<- function() 
{
    # initialize 
    xstart <- c(Rm=4.74,R=267,DR=0,DRN=0)
    return(xstart)
}

###################
# Script
# ------

# Observation time series (using function seq)
ObservationTimes<-seq(0,24*3,by=1)

# Load model parameters
MyParms<-c(
  ksyn_rm=2.9,
  IC50=26.2,
  kon=0.00329,
  kt=0.63,
  kre=0.57,
  rf=0.49,
  kdgr_r=0.0572,
  kdgr_rm=0.612,
  ksyn_r=0,
  D=20*1000/374.471    #(20/374.471)*1000
)

# call lsoda and store result in out
out <- ode(
   y=systemInit(), 
   times=ObservationTimes, 								
   func=fifthGenerationModel, 
   parms=MyParms
) 

# Convert output to dataframe for easy post-processing using $ notation
out<-as.data.frame(out)
out$tot <- out$R + out$DR + out$DRN

pdf("week4_noSynthesis.pdf")
# Make plot of concentration of receptor mRNa and free receptor density
par(mfrow=c(2,3))
  plot(out$time,out$Rm,ylim = c(0,5), xlab="Time",ylab="receptor mRNA",type="l",lwd=2)
  plot(out$time,out$R, ylim = c(0,300), xlab="Time",ylab="free receptor density",type="l",lwd=2)
  plot(out$time,out$DR, ylim = c(0,50), xlab="Time",ylab="crug-receptor complex",type="l",lwd=2)
  plot(out$time,out$DRN, ylim = c(0,50), xlab="Time",ylab="activated receptor complex",type="l",lwd=2)
  plot(out$time,out$tot, ylim = c(0,300), xlab="Time",ylab="total receptor density",type="l",lwd=2)

dev.off()

head(out)
tail(out)