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

# Load model parameters
MyParms_kon<-cbind(
  ksyn_rm=rep(2.9,times=length(factors)),
  IC50=rep(26.2,times=length(factors)),
  kon=0.00329*factors,
  kt=rep(0.63,times=length(factors)),
  kre=rep(0.57,times=length(factors)),
  rf=rep(0.49,times=length(factors)),
  kdgr_r=rep(0.0572,times=length(factors)),
  kdgr_rm=rep(0.612,times=length(factors)),
  ksyn_r=rep(3.22,times=length(factors)),
  D=rep(20*1000/374.471,times=length(factors))   #(20/374.471)*1000
)

MyParms_kre<-cbind(
  ksyn_rm=rep(2.9,times=length(factors)),
  IC50=rep(26.2,times=length(factors)),
  kon=rep(0.00329,times=length(factors)),
  kt=rep(0.63,times=length(factors)),
  kre=0.57*factors,
  rf=rep(0.49,times=length(factors)),
  kdgr_r=rep(0.0572,times=length(factors)),
  kdgr_rm=rep(0.612,times=length(factors)),
  ksyn_r=rep(3.22,times=length(factors)),
  D=rep(20*1000/374.471,times=length(factors))   #(20/374.471)*1000
)

# simulations for different kon

for (i in 1:length(factors)){
# call lsoda and store result in out
  out <- ode(
     y=systemInit(), 
     times=ObservationTimes,   							
    func=fifthGenerationModel, 
    parms=MyParms_kon[i,]
  )  
  # Convert output to dataframe for easy post-processing using $ notation
  out<-as.data.frame(out)
  out$kon <- MyParms_kon[i,c("kon")]
  
  if (i>1)
   {  
       simset_kon<-rbind(simset_kon,out)
   }
   else
   {
       simset_kon<-out
   }
}

simset_kon$tot <- simset_kon$R + simset_kon$DR + simset_kon$DRN

# simulations for different kre

for (i in 1:length(factors)){
# call lsoda and store result in out
  out <- ode(
     y=systemInit(), 
     times=ObservationTimes,     						
    func=fifthGenerationModel, 
    parms=MyParms_kre[i,]
  )  
  # Convert output to dataframe for easy post-processing using $ notation
  out<-as.data.frame(out)
  out$kre <- MyParms_kre[i,c("kre")]
  
  if (i>1)
   {  
       simset_kre<-rbind(simset_kre,out)
   }
   else
   {
       simset_kre<-out
   }
}

simset_kre$tot <- simset_kre$R + simset_kre$DR + simset_kre$DRN


pdf("week4_kon_factors.pdf")

# Make plots of concentration of receptor mRNa and free receptor density
par(mfrow=c(2,3))
  plot(simset_kon$time[simset_kon$kon==MyParms_kon[1,c("kon")]],simset_kon$Rm[simset_kon$kon==MyParms_kon[1,c("kon")]],ylim = c(0,5), xlab="Time",ylab="receptor mRNA",type="l",lwd=2, col="orange")
  lines(simset_kon$time[simset_kon$kon==MyParms_kon[2,c("kon")]],simset_kon$Rm[simset_kon$kon==MyParms_kon[2,c("kon")]],ylim = c(0,5), type="l",lwd=2, col="red")
  lines(simset_kon$time[simset_kon$kon==MyParms_kon[3,c("kon")]],simset_kon$Rm[simset_kon$kon==MyParms_kon[3,c("kon")]],ylim = c(0,5), type="l",lwd=2, col="black")
  lines(simset_kon$time[simset_kon$kon==MyParms_kon[4,c("kon")]],simset_kon$Rm[simset_kon$kon==MyParms_kon[4,c("kon")]],ylim = c(0,5), type="l",lwd=2, col="darkblue")
  lines(simset_kon$time[simset_kon$kon==MyParms_kon[5,c("kon")]],simset_kon$Rm[simset_kon$kon==MyParms_kon[5,c("kon")]],ylim = c(0,5), type="l",lwd=2, col="lightblue")

  plot(simset_kon$time[simset_kon$kon==MyParms_kon[1,c("kon")]],simset_kon$R[simset_kon$kon==MyParms_kon[1,c("kon")]], ylim = c(0,300), xlab="Time",ylab="free receptor density",type="l",lwd=2, col="orange")
  lines(simset_kon$time[simset_kon$kon==MyParms_kon[2,c("kon")]],simset_kon$R[simset_kon$kon==MyParms_kon[2,c("kon")]], ylim = c(0,300), type="l",lwd=2, col="red")
  lines(simset_kon$time[simset_kon$kon==MyParms_kon[3,c("kon")]],simset_kon$R[simset_kon$kon==MyParms_kon[3,c("kon")]], ylim = c(0,300), type="l",lwd=2, col="black")
  lines(simset_kon$time[simset_kon$kon==MyParms_kon[4,c("kon")]],simset_kon$R[simset_kon$kon==MyParms_kon[4,c("kon")]], ylim = c(0,300), type="l",lwd=2, col="darkblue")
  lines(simset_kon$time[simset_kon$kon==MyParms_kon[5,c("kon")]],simset_kon$R[simset_kon$kon==MyParms_kon[5,c("kon")]], ylim = c(0,300), type="l",lwd=2, col="lightblue")

  plot(simset_kon$time[simset_kon$kon==MyParms_kon[1,c("kon")]],simset_kon$DR[simset_kon$kon==MyParms_kon[1,c("kon")]], ylim = c(0,130), xlab="Time",ylab="drug-receptor complex",type="l",lwd=2, col="orange")
  lines(simset_kon$time[simset_kon$kon==MyParms_kon[2,c("kon")]],simset_kon$DR[simset_kon$kon==MyParms_kon[2,c("kon")]], ylim = c(0,130), type="l",lwd=2, col="red")
  lines(simset_kon$time[simset_kon$kon==MyParms_kon[3,c("kon")]],simset_kon$DR[simset_kon$kon==MyParms_kon[3,c("kon")]], ylim = c(0,130), type="l",lwd=2, col="black")
  lines(simset_kon$time[simset_kon$kon==MyParms_kon[4,c("kon")]],simset_kon$DR[simset_kon$kon==MyParms_kon[4,c("kon")]], ylim = c(0,130), type="l",lwd=2, col="darkblue")
  lines(simset_kon$time[simset_kon$kon==MyParms_kon[5,c("kon")]],simset_kon$DR[simset_kon$kon==MyParms_kon[5,c("kon")]], ylim = c(0,130), type="l",lwd=2, col="lightblue")

  plot(simset_kon$time[simset_kon$kon==MyParms_kon[1,c("kon")]],simset_kon$DRN[simset_kon$kon==MyParms_kon[1,c("kon")]], ylim = c(0,100), xlab="Time",ylab="activated drug-receptor complex",type="l",lwd=2, col="orange")
  lines(simset_kon$time[simset_kon$kon==MyParms_kon[2,c("kon")]],simset_kon$DRN[simset_kon$kon==MyParms_kon[2,c("kon")]], ylim = c(0,100), type="l",lwd=2, col="red")
  lines(simset_kon$time[simset_kon$kon==MyParms_kon[3,c("kon")]],simset_kon$DRN[simset_kon$kon==MyParms_kon[3,c("kon")]], ylim = c(0,100), type="l",lwd=2, col="black")
  lines(simset_kon$time[simset_kon$kon==MyParms_kon[4,c("kon")]],simset_kon$DRN[simset_kon$kon==MyParms_kon[4,c("kon")]], ylim = c(0,100), type="l",lwd=2, col="darkblue")
  lines(simset_kon$time[simset_kon$kon==MyParms_kon[5,c("kon")]],simset_kon$DRN[simset_kon$kon==MyParms_kon[5,c("kon")]], ylim = c(0,100), type="l",lwd=2, col="lightblue")

  plot(simset_kon$time[simset_kon$kon==MyParms_kon[1,c("kon")]],simset_kon$tot[simset_kon$kon==MyParms_kon[1,c("kon")]], ylim = c(0,300), xlab="Time",ylab="total receptor density",type="l",lwd=2, col="orange")
  lines(simset_kon$time[simset_kon$kon==MyParms_kon[2,c("kon")]],simset_kon$tot[simset_kon$kon==MyParms_kon[2,c("kon")]], ylim = c(0,300), type="l",lwd=2, col="red")
  lines(simset_kon$time[simset_kon$kon==MyParms_kon[3,c("kon")]],simset_kon$tot[simset_kon$kon==MyParms_kon[3,c("kon")]], ylim = c(0,300), type="l",lwd=2, col="black")
  lines(simset_kon$time[simset_kon$kon==MyParms_kon[4,c("kon")]],simset_kon$tot[simset_kon$kon==MyParms_kon[4,c("kon")]], ylim = c(0,300), type="l",lwd=2, col="darkblue")
  lines(simset_kon$time[simset_kon$kon==MyParms_kon[5,c("kon")]],simset_kon$tot[simset_kon$kon==MyParms_kon[5,c("kon")]], ylim = c(0,300), type="l",lwd=2, col="lightblue")

dev.off()

pdf("week4_kre_factors.pdf")

# Make plots of concentration of receptor mRNa and free receptor density
par(mfrow=c(2,3))
  plot(simset_kre$time[simset_kre$kre==MyParms_kre[1,c("kre")]],simset_kre$Rm[simset_kre$kre==MyParms_kre[1,c("kre")]],ylim = c(0,5), xlab="Time",ylab="receptor mRNA",type="l",lwd=2, col="orange")
  lines(simset_kre$time[simset_kre$kre==MyParms_kre[2,c("kre")]],simset_kre$Rm[simset_kre$kre==MyParms_kre[2,c("kre")]],ylim = c(0,5), type="l",lwd=2, col="red")
  lines(simset_kre$time[simset_kre$kre==MyParms_kre[3,c("kre")]],simset_kre$Rm[simset_kre$kre==MyParms_kre[3,c("kre")]],ylim = c(0,5), type="l",lwd=2, col="black")
  lines(simset_kre$time[simset_kre$kre==MyParms_kre[4,c("kre")]],simset_kre$Rm[simset_kre$kre==MyParms_kre[4,c("kre")]],ylim = c(0,5), type="l",lwd=2, col="darkblue")
  lines(simset_kre$time[simset_kre$kre==MyParms_kre[5,c("kre")]],simset_kre$Rm[simset_kre$kre==MyParms_kre[5,c("kre")]],ylim = c(0,5), type="l",lwd=2, col="lightblue")

  plot(simset_kre$time[simset_kre$kre==MyParms_kre[1,c("kre")]],simset_kre$R[simset_kre$kre==MyParms_kre[1,c("kre")]], ylim = c(0,300), xlab="Time",ylab="free receptor density",type="l",lwd=2, col="orange")
  lines(simset_kre$time[simset_kre$kre==MyParms_kre[2,c("kre")]],simset_kre$R[simset_kre$kre==MyParms_kre[2,c("kre")]], ylim = c(0,300), type="l",lwd=2, col="red")
  lines(simset_kre$time[simset_kre$kre==MyParms_kre[3,c("kre")]],simset_kre$R[simset_kre$kre==MyParms_kre[3,c("kre")]], ylim = c(0,300), type="l",lwd=2, col="black")
  lines(simset_kre$time[simset_kre$kre==MyParms_kre[4,c("kre")]],simset_kre$R[simset_kre$kre==MyParms_kre[4,c("kre")]], ylim = c(0,300), type="l",lwd=2, col="darkblue")
  lines(simset_kre$time[simset_kre$kre==MyParms_kre[5,c("kre")]],simset_kre$R[simset_kre$kre==MyParms_kre[5,c("kre")]], ylim = c(0,300), type="l",lwd=2, col="lightblue")

  plot(simset_kre$time[simset_kre$kre==MyParms_kre[1,c("kre")]],simset_kre$DR[simset_kre$kre==MyParms_kre[1,c("kre")]], ylim = c(0,60), xlab="Time",ylab="drug-receptor complex",type="l",lwd=2, col="orange")
  lines(simset_kre$time[simset_kre$kre==MyParms_kre[2,c("kre")]],simset_kre$DR[simset_kre$kre==MyParms_kre[2,c("kre")]], ylim = c(0,60), type="l",lwd=2, col="red")
  lines(simset_kre$time[simset_kre$kre==MyParms_kre[3,c("kre")]],simset_kre$DR[simset_kre$kre==MyParms_kre[3,c("kre")]], ylim = c(0,60), type="l",lwd=2, col="black")
  lines(simset_kre$time[simset_kre$kre==MyParms_kre[4,c("kre")]],simset_kre$DR[simset_kre$kre==MyParms_kre[4,c("kre")]], ylim = c(0,60), type="l",lwd=2, col="darkblue")
  lines(simset_kre$time[simset_kre$kre==MyParms_kre[5,c("kre")]],simset_kre$DR[simset_kre$kre==MyParms_kre[5,c("kre")]], ylim = c(0,60), type="l",lwd=2, col="lightblue")

  plot(simset_kre$time[simset_kre$kre==MyParms_kre[1,c("kre")]],simset_kre$DRN[simset_kre$kre==MyParms_kre[1,c("kre")]], ylim = c(0,130), xlab="Time",ylab="activated drug-receptor complex",type="l",lwd=2, col="orange")
  lines(simset_kre$time[simset_kre$kre==MyParms_kre[2,c("kre")]],simset_kre$DRN[simset_kre$kre==MyParms_kre[2,c("kre")]], ylim = c(0,130), type="l",lwd=2, col="red")
  lines(simset_kre$time[simset_kre$kre==MyParms_kre[3,c("kre")]],simset_kre$DRN[simset_kre$kre==MyParms_kre[3,c("kre")]], ylim = c(0,130), type="l",lwd=2, col="black")
  lines(simset_kre$time[simset_kre$kre==MyParms_kre[4,c("kre")]],simset_kre$DRN[simset_kre$kre==MyParms_kre[4,c("kre")]], ylim = c(0,130), type="l",lwd=2, col="darkblue")
  lines(simset_kre$time[simset_kre$kre==MyParms_kre[5,c("kre")]],simset_kre$DRN[simset_kre$kre==MyParms_kre[5,c("kre")]], ylim = c(0,130), type="l",lwd=2, col="lightblue")

  plot(simset_kre$time[simset_kre$kre==MyParms_kre[1,c("kre")]],simset_kre$tot[simset_kre$kre==MyParms_kre[1,c("kre")]], ylim = c(0,300), xlab="Time",ylab="total receptor density",type="l",lwd=2, col="orange")
  lines(simset_kre$time[simset_kre$kre==MyParms_kre[2,c("kre")]],simset_kre$tot[simset_kre$kre==MyParms_kre[2,c("kre")]], ylim = c(0,300), type="l",lwd=2, col="red")
  lines(simset_kre$time[simset_kre$kre==MyParms_kre[3,c("kre")]],simset_kre$tot[simset_kre$kre==MyParms_kre[3,c("kre")]], ylim = c(0,300), type="l",lwd=2, col="black")
  lines(simset_kre$time[simset_kre$kre==MyParms_kre[4,c("kre")]],simset_kre$tot[simset_kre$kre==MyParms_kre[4,c("kre")]], ylim = c(0,300), type="l",lwd=2, col="darkblue")
  lines(simset_kre$time[simset_kre$kre==MyParms_kre[5,c("kre")]],simset_kre$tot[simset_kre$kre==MyParms_kre[5,c("kre")]], ylim = c(0,300), type="l",lwd=2, col="lightblue")

dev.off()

head(simset_kon)
tail(simset_kon)

head(simset_kre)
tail(simset_kre)