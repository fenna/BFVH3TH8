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
ObservationTimes<-seq(0,24*2,by=1)

factors <- c(1/5,1/2,1,2,5)
ksyn_rm=2.9*factors
kdgr_rm=ksyn_rm/4.74

# Load model parameters
MyParms_ksyn_rm<-cbind(
  ksyn_rm=ksyn_rm,
  IC50=rep(26.2,times=length(factors)),
  kon=rep(0.00329,times=length(factors)),
  kt=rep(0.63,times=length(factors)),
  kre=rep(0.57,times=length(factors)),
  rf=rep(0.49,times=length(factors)),
  kdgr_r=rep(0.0572,times=length(factors)),
  kdgr_rm=kdgr_rm,
  ksyn_r=rep(3.22,times=length(factors)),
  D=rep(20*1000/374.471,times=length(factors))   #(20/374.471)*1000
)


# simulations for different ksyn_rm

for (i in 1:length(factors)){
# call lsoda and store result in out
  out <- ode(
     y=systemInit(), 
     times=ObservationTimes,   							
    func=fifthGenerationModel, 
    parms=MyParms_ksyn_rm[i,]
  )  
  # Convert output to dataframe for easy post-processing using $ notation
  out<-as.data.frame(out)
  out$ksyn_rm <- MyParms_ksyn_rm[i,c("ksyn_rm")]
  
  if (i>1)
   {  
       simset_ksyn_rm<-rbind(simset_ksyn_rm,out)
   }
   else
   {
       simset_ksyn_rm<-out
   }
}

simset_ksyn_rm$tot <- simset_ksyn_rm$R + simset_ksyn_rm$DR + simset_ksyn_rm$DRN


pdf("week4_ksyn_rm_factors.pdf")

# Make plots of concentration of receptor mRNa and free receptor density
par(mfrow=c(2,3))
  plot(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[1,c("ksyn_rm")]],simset_ksyn_rm$Rm[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[1,c("ksyn_rm")]],ylim = c(0,5), xlab="Time",ylab="receptor mRNA",type="l",lwd=2, col="orange")
  lines(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[2,c("ksyn_rm")]],simset_ksyn_rm$Rm[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[2,c("ksyn_rm")]],ylim = c(0,5), type="l",lwd=2, col="red")
  lines(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[3,c("ksyn_rm")]],simset_ksyn_rm$Rm[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[3,c("ksyn_rm")]],ylim = c(0,5), type="l",lwd=2, col="black")
  lines(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[4,c("ksyn_rm")]],simset_ksyn_rm$Rm[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[4,c("ksyn_rm")]],ylim = c(0,5), type="l",lwd=2, col="darkblue")
  lines(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[5,c("ksyn_rm")]],simset_ksyn_rm$Rm[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[5,c("ksyn_rm")]],ylim = c(0,5), type="l",lwd=2, col="lightblue")

  plot(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[1,c("ksyn_rm")]],simset_ksyn_rm$R[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[1,c("ksyn_rm")]], ylim = c(0,300), xlab="Time",ylab="free receptor density",type="l",lwd=2, col="orange")
  lines(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[2,c("ksyn_rm")]],simset_ksyn_rm$R[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[2,c("ksyn_rm")]], ylim = c(0,300), type="l",lwd=2, col="red")
  lines(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[3,c("ksyn_rm")]],simset_ksyn_rm$R[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[3,c("ksyn_rm")]], ylim = c(0,300), type="l",lwd=2, col="black")
  lines(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[4,c("ksyn_rm")]],simset_ksyn_rm$R[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[4,c("ksyn_rm")]], ylim = c(0,300), type="l",lwd=2, col="darkblue")
  lines(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[5,c("ksyn_rm")]],simset_ksyn_rm$R[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[5,c("ksyn_rm")]], ylim = c(0,300), type="l",lwd=2, col="lightblue")

  plot(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[1,c("ksyn_rm")]],simset_ksyn_rm$DR[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[1,c("ksyn_rm")]], ylim = c(0,60), xlab="Time",ylab="drug-receptor complex",type="l",lwd=2, col="orange")
  lines(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[2,c("ksyn_rm")]],simset_ksyn_rm$DR[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[2,c("ksyn_rm")]], ylim = c(0,60), type="l",lwd=2, col="red")
  lines(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[3,c("ksyn_rm")]],simset_ksyn_rm$DR[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[3,c("ksyn_rm")]], ylim = c(0,60), type="l",lwd=2, col="black")
  lines(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[4,c("ksyn_rm")]],simset_ksyn_rm$DR[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[4,c("ksyn_rm")]], ylim = c(0,60), type="l",lwd=2, col="darkblue")
  lines(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[5,c("ksyn_rm")]],simset_ksyn_rm$DR[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[5,c("ksyn_rm")]], ylim = c(0,60), type="l",lwd=2, col="lightblue")

  plot(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[1,c("ksyn_rm")]],simset_ksyn_rm$DRN[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[1,c("ksyn_rm")]], ylim = c(0,60), xlab="Time",ylab="activated drug-receptor complex",type="l",lwd=2, col="orange")
  lines(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[2,c("ksyn_rm")]],simset_ksyn_rm$DRN[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[2,c("ksyn_rm")]], ylim = c(0,60), type="l",lwd=2, col="red")
  lines(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[3,c("ksyn_rm")]],simset_ksyn_rm$DRN[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[3,c("ksyn_rm")]], ylim = c(0,60), type="l",lwd=2, col="black")
  lines(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[4,c("ksyn_rm")]],simset_ksyn_rm$DRN[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[4,c("ksyn_rm")]], ylim = c(0,60), type="l",lwd=2, col="darkblue")
  lines(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[5,c("ksyn_rm")]],simset_ksyn_rm$DRN[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[5,c("ksyn_rm")]], ylim = c(0,60), type="l",lwd=2, col="lightblue")

  plot(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[1,c("ksyn_rm")]],simset_ksyn_rm$tot[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[1,c("ksyn_rm")]], ylim = c(0,300), xlab="Time",ylab="total receptor density",type="l",lwd=2, col="orange")
  lines(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[2,c("ksyn_rm")]],simset_ksyn_rm$tot[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[2,c("ksyn_rm")]], ylim = c(0,300), type="l",lwd=2, col="red")
  lines(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[3,c("ksyn_rm")]],simset_ksyn_rm$tot[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[3,c("ksyn_rm")]], ylim = c(0,300), type="l",lwd=2, col="black")
  lines(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[4,c("ksyn_rm")]],simset_ksyn_rm$tot[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[4,c("ksyn_rm")]], ylim = c(0,300), type="l",lwd=2, col="darkblue")
  lines(simset_ksyn_rm$time[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[5,c("ksyn_rm")]],simset_ksyn_rm$tot[simset_ksyn_rm$ksyn_rm==MyParms_ksyn_rm[5,c("ksyn_rm")]], ylim = c(0,300), type="l",lwd=2, col="lightblue")

dev.off()



head(simset_ksyn_rm)
tail(simset_ksyn_rm)