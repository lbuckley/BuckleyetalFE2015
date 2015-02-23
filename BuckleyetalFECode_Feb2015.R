#CODE FOR
#Buckley, L. B., Ehrenberger, J. C., Angilletta, M. J. (2015), Thermoregulatory behaviour limits local adaptation of thermal niches and confers sensitivity to climate change. Functional Ecology. doi: 10.1111/1365-2435.12406

library(truncnorm) 
#library(TeachingDemos)

#source biophysical model
setwd("C:\\Users\\Buckley\\Google Drive\\BuckleyLab\\LabMeetings\\BuckleyetalFE2015\\")
source("BiophysModel_20Apr.R")
#source zenith angle calculating function
source("ZenithAngleFunction.R")
#source dtr function
source("DTRfunction_Feb2013.R")

#Performance Curve Function from Deutsch et al. 2008
TPC= function(T,Topt=33,CTmin=10.45, CTmax=42.62){
  F=rep(NA, length(T))
  sigma= (Topt-CTmin)/4
  ind=which(T<=Topt)
  F[ind]= exp(-((T[ind]-Topt)/(2*sigma))^2) 
  ind=which(T>Topt)
  F[ind]= 1- ((T[ind]-Topt)/(Topt-CTmax))^2
  return(F)
}

#Performance curve function from Sears and Angilletta
# Create and assign variables and parameters
#temp      #body temperature in degrees Celsius
#shift     #Horizontal shift in thermal performance curve, which determines the thermal optimum and critical thermal limits
#breadth   #A parameter that determines the thermal breadth of performace
#skew      #A parameter that influences the skew of the thermal performance curve
#tolerance   #A parameter that sets the maximal breadth of the thermal performance curve

#FUNCTIONS FOR OPTIMIZING
#if aran=1, include araneous factor
TPC.beta= function(T, shift=-1, breadth=0.1, aran=0, tolerance= 43, skew=0.7){ 
T = T + 273.15 #Convert temperature in degrees Celsius to Kelvin 
shift= shift + 273.15 #Convert temperature in degrees Celsius to Kelvin         
z=rep(0.01, length(T))
sel= which(T-shift>=0 & T-shift<=tolerance)
z[sel]= ((T[sel]-shift)/tolerance)^(skew/breadth-1)*(1-(T[sel]-shift)/tolerance)^((1-skew)/breadth-1) / beta(skew/breadth,(1-skew)/breadth) 
if(aran==1) z[sel]=z[sel]*exp(-0.6/(T[sel]*8.61734*10^(-5)))*10^10 #add scaling factor
return(z)
}

z.beta.sum= function(temp, shift,breadth, aran=0, tolerance= 43, skew=0.7){
temp[1] = temp[1] + 273.15 #Convert temperature in degrees Celsius to Kelvin 
shift= shift + 273.15 #Convert temperature in degrees Celsius to Kelvin
z=0.01         
if(temp[1]-shift>=0 & temp[1]-shift<=tolerance) z= ((temp[1]-shift)/tolerance)^(skew/breadth-1)*(1-(temp[1]-shift)/tolerance)^((1-skew)/breadth-1) / beta(skew/breadth,(1-skew)/breadth) #*exp(-0.6/(temp[1]*8.61734*10^(-5)))
if(aran==1) z=z*exp(-0.6/(temp[1]*8.61734*10^(-5)))*10^10
z=z*temp[2]
return(z)
}

z.beta.prod= function(temp, shift,breadth, aran=0, tolerance= 43, skew=0.7){
#scale instances of temp
temp[2]=temp[2]/10
temp[1] = temp[1] + 273.15 #Convert temperature in degrees Celsius to Kelvin 
shift= shift + 273.15 #Convert temperature in degrees Celsius to Kelvin
z=0.01         
if(temp[1]-shift>=0 & temp[1]-shift<=tolerance) z= ((temp[1]-shift)/tolerance)^(skew/breadth-1)*(1-(temp[1]-shift)/tolerance)^((1-skew)/breadth-1) / beta(skew/breadth,(1-skew)/breadth) #*exp(-0.6/(temp[1]*8.61734*10^(-5)))
if(aran==1) z=z*exp(-0.6/(temp[1]*8.61734*10^(-5)))*10^10
z=z^temp[2] 
if(z==0) z=10^-20
return(z)
}

TPC.beta.sum= function(bc, temps=temps, aran=0){
shift= bc[1]
breadth= bc[2]

sum(apply(temps, MARGIN=1, FUN=z.beta.sum, shift=shift, breadth=breadth, aran=0), na.rm=TRUE)
}

TPC.beta.prod= function(bc, temps=temps, aran=0){
shift= bc[1]
breadth= bc[2]

mean(log(unlist(apply(temps, MARGIN=1, FUN=z.beta.prod, shift=shift, breadth=breadth, aran=0))), na.rm=TRUE)
}


#find CT min and max and Topt
TPC.beta.minmax= function(T, shift, beta, aran=0, tolerance=43, skew=0.7){
CTs= rep(NA,3)
v= TPC.beta(T,shift,beta, aran, tolerance=tolerance, skew=skew)
CTs[2]= T[which.max(v)]
v1= v/ v[which.max(v)]*100
perf= which(v1>1)
CTs[1]= T[perf[1]]
CTs[3]= T[max(perf)]
return(CTs)
}

#THERMOREG FUNCTIONS
TR= function(hr, Topt, TeS.day1, Te.day1){
all.temps= seq( TeS.day1[hr], Te.day1[hr], 0.1) 
T.TR= all.temps[which.min(abs(all.temps-Topt))]
return(T.TR)
}

noTR= function(hr, TeS.day1, Te.day1){
runif(100,TeS.day1[hr], Te.day1[hr] )
}

costTR= function(hr, TPC.cTR1, TPC.TR1, T.cTR1, T.TR1){
T.cTR1[which(TPC.cTR1[,hr]<TPC.TR1[hr]),hr]= T.TR1[hr]
return(T.cTR1[,hr])
}

#### ADD SCENARIOS TO AVOID EXTREMES


noTR.noExt.unif= function(hr, TeS.day1, Te.day1,CTmin=10.45, CTmax=42.62){

lowlim= TeS.day1[hr]
if(Te.day1[hr]>CTmin) lowlim= max(TeS.day1[hr],CTmin)
uplim= Te.day1[hr]
if(TeS.day1[hr]<CTmax) uplim= min(Te.day1[hr], CTmax)

runif(100, lowlim, uplim )
}

#ADD OTHER DISTRIBUTIONS
#MODAL
noTR.noExt.modal= function(hr, TeS.day1, Te.day1,CTmin=10.45, CTmax=42.62){

lowlim= TeS.day1[hr]
if(Te.day1[hr]>CTmin) lowlim= max(TeS.day1[hr],CTmin)
uplim= Te.day1[hr]
if(TeS.day1[hr]<CTmax) uplim= min(Te.day1[hr], CTmax)

rtruncnorm(100, a=lowlim, b=uplim, mean = (lowlim+uplim)/2, sd = (uplim-lowlim)/6)
}

#BIMODAL
noTR.noExt.bimodal= function(hr, TeS.day1, Te.day1,CTmin=10.45, CTmax=42.62){

lowlim= TeS.day1[hr]
if(Te.day1[hr]>CTmin) lowlim= max(TeS.day1[hr],CTmin)
uplim= Te.day1[hr]
if(TeS.day1[hr]<CTmax) uplim= min(Te.day1[hr], CTmax)

mu1 <- (lowlim+uplim)/3   
mu2 <-  2*(lowlim+uplim)/3 
sd1 <- (uplim-lowlim)/10 #12 #change to 10 so some overlap

ns=1000
prop= runif(ns, 0, 1)
n1= length(which(prop<0.5))
x <- c(rtruncnorm(n1, a=lowlim, b=Inf, mu1, sd1), rtruncnorm(ns-n1, a=-Inf, b=uplim, mu2, sd1)) 
x.new <- sample(x, size = ns, replace=TRUE)
}

#----------------

#CURRENTLY USING NUMERIC OPTIMIZATION, COULD ALSO USE OPTIM FUNCTION
#optim(c(20,10), TPC.gaus.optim, T=temps.nTR, control = list(fnscale = -1))

#BIOPHYSICAL MODELING ***************************************
# Calculates operative environmental temperature for Sceloporus undulatus in the US (from Buckley 2008, AmNat)

# read climate data
#10' DATA
#subset to grid cells within range
climdata=read.csv('ClimateData_2160x684_withinScelUndu.csv');

#read surface temperature data
tsurfdata=read.csv("Tsurface_NARRall.csv");
#match to climate data
match1= match(climdata$GridID, tsurfdata$GridID)
tsurfdata= tsurfdata[match1,]

lat= climdata[,3] ; #latitude in degrees
lon= climdata[,2] ; #longitude in degrees

#use surface T data 
Tm= tsurfdata[,2:13]; #mean of an average day of each month
Tr= tsurfdata[,14:25]; #diurnal temperature range of an average day of each month
Wind= climdata[,4]; #mean wind speed m/s
Wind[]=0.1 #assume 0.1 m/s wind
Albedo= climdata[,5:8]; #Albedo percentage for Jan / Apr / Jul / Oct
Elev= climdata[,33]; #mean elevation (m)
TSmax=climdata[,34:45];  #max soil T, mean 14:00hr temperature for five days in the middle of each month (K), from LDAS
TSr=climdata[,46:57];  #range between 02:00 and 14:00hr temperature for five days in the middle of each month (K), from LDAS
Tairm= climdata[,9:20]; #mean of an average day of each month
Tairr= climdata[,21:32]; #diurnal temperature range of an average day of each month

Topt=33 #optimal temperature
#__________________________________________________

#LIZARD DATA
lizdat=read.csv("POPmeans_20Feb2013.csv") 

library(maptools)
#split in quantiles
tb1<-cut(lizdat$meanPBT,breaks=10)
lev<-levels(tb1)
cols<-rainbow(10)

lizdat$TTB= lizdat$meanCTmax - lizdat$meanCTmin
TTB=32.17 #mean TTB= 32.17
PB= 6.4 #performance breadth
#-------------------------------
#abbreviate to AZ and TX
#lizdat=lizdat[c(1,8),]

svl= lizdat[, 8] #mean of max male and female SVL (Ord & Blumstein 2002, Herrel 2002)
svl[]=67 #Use 67 mm, Adult body length of females is 67 mm averaged across all populations
mass= 3.55*10^-5*(svl)^3.00 #Tinkle and Ballinger 1972 (g)

#------------------------------------
# data structures

#storage for Te
Te= array(NA, dim=c( 1,12, 24,nrow(lizdat)))#month, hour, spec
TeS= array(NA, dim=c( 1,12, 24,nrow(lizdat)))#month, hour, spec
Te.day= array(NA, dim=c( 1,12, 24,nrow(lizdat)))#month, hour, spec
TeS.day= array(NA, dim=c( 1,12, 24,nrow(lizdat)))#month, hour, spec

Topts= array(NA, dim=c( nrow(Tm),6,10,2)) #cells, iter, metrics, araneous?
Fmax= array(NA, dim=c( nrow(Tm),6,2)) #cells, iter, araneous?

dimnames(Topts)[[3]]= c("shift.sum","breadth.sum", "CTmin.sum", "Topt.sum", "CTmax.sum", "shift.prod","breadth.prod", "CTmin.prod", "Topt.prod", "CTmax.prod")#------------------------------------
# Calculate Operative Environmental Temperatures 

J.all=c(15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349)

#define combinations
shift=seq(-9,5,1)
breadth=seq(0.07,0.17,0.02)
comb= expand.grid(shift, breadth)
 
#for (cellk in 5000:5100){
for (cellk in 1:nrow(Tm)){ # cell counter
     for (monthk in 5:9){ # month counter	 #JUST MAY TO SEPT

	hourk=1:24
      J= J.all[monthk]  # Julian calendar day   

	#Calculate zenith in degrees
	out=zenith(J, lat[cellk], lon[cellk], hourk)
	psi_d=out[1:24]
	psi_d[psi_d>=89]=89 #set zenith position below horizon to psi=89degrees
      #convert to radians
	psi_rad= psi_d * pi/180            

	#daylength
	sunrise=out[25]
	sunset=out[26]
	dayhrs= hourk[which(hourk<sunset & hourk>sunrise)]
	nighthrs= hourk[-dayhrs]	

	# Calculate the diurnal temperature trend, Parton and Logan 1981
	#surface temp
	Tx= Tm[cellk, monthk] + Tr[cellk, monthk] / 2  # maximum daily temperature
      Tn= Tm[cellk, monthk] - Tr[cellk, monthk] / 2 # minimum daily temperature

	Ta= sapply(hourk, Thour, Tmx=Tx, Tmn=Tn, alpha=alpha, beta=beta, gamma=gamma, tr=sunrise,ts=sunset)

	#air temp
	Tairx= Tairm[cellk, monthk] + Tairr[cellk, monthk] / 2  # maximum daily temperature
      Tairn= Tairm[cellk, monthk] - Tairr[cellk, monthk] / 2 # minimum daily temperature
                
	Taira  = sapply(hourk, Thour, Tmx=Tairx, Tmn=Tairn, alpha=alpha, beta=beta, gamma=gamma, tr=sunrise,ts=sunset)
	
#__________________________________________________
#Calculate Te 
albk= floor (monthk/3.01)+1 # compute season for albedo measurement
rho_S= Albedo[cellk, albk]/100 # (Table 12.2)

# Loop through species
#for(Spec in 1:nrow(lizdat) ){ #species counter
Spec=1    
        
#no CC
Tdat<-biophys(Taira, Ta, WIND=Wind[cellk], SVL=svl[Spec], MASS=mass[Spec], psi=psi_rad, rho_S, elevation=Elev[cellk],J)

#find warmest and coolest temps
Te_a= Tdat[1:24] #Te based on air in sun
TeS_a= Tdat[25:48] #Te based on air in shade
Te_g= Tdat[49:72] #Te based on ground in sun
TeS_g= Tdat[73:96] #Te based on ground in shade
 
hot= which(Te_g>Te_a)
Te_a[hot]=Te_g[hot]
cool= which(TeS_g<TeS_a)
TeS_a[cool]=TeS_g[cool]

#Te[1,monthk,,Spec]  = Te_a
#TeS[1,monthk,,Spec]= TeS_a
Te.day[1,monthk,dayhrs,Spec]= Te_a[dayhrs]
TeS.day[1,monthk,dayhrs,Spec]= TeS_a[dayhrs]
#} #end spec loop

} # end month loop

if( round(cellk/200)==cellk/200) print(cellk)

#Calculate optimum without TR
# 100 interations of choosing uniform value between max and min on a day

temps.TR=NA
temps.nTR=NA
temps.cTR=NA
temps.nTR.noExt.unif=NA
temps.nTR.noExt.modal=NA
temps.nTR.noExt.bimodal=NA

for(mon in 5:9){

hrs= which(!is.na(TeS.day[1,mon,,Spec]))

#WITH TR
#Assume thermoregulate to mean PBT (32) if possible
T.TR=sapply(hrs, FUN=TR, Topt=32, TeS.day1=TeS.day[1,mon,,Spec], Te.day1=Te.day[1,mon,,Spec])   
temps.TR= c(temps.TR, T.TR)

#WITHOUT TR
#noTR(12, TeS.day[1,mon,,Spec], Te.day[1,mon,,Spec]) 
temps=lapply(hrs, FUN=noTR, TeS.day1=TeS.day[1,mon,,Spec], Te.day1=Te.day[1,mon,,Spec])  
temps.mat= sapply(temps,c) #convert to matrix
temps.vec= sapply(temps.mat,c) #convert to vector
temps.nTR=c(temps.nTR, temps.vec ) 

#COST OF THERMOREG, assume 50% performance detriment
TPC.TR= TPC(T.TR)*0.5
T.cTR= temps.mat
TPC.cTR= apply(T.cTR, FUN=TPC, MARGIN=2)

T.cTR.mat=sapply(1:length(hrs), FUN=costTR, TPC.cTR1=TPC.cTR, TPC.TR1=TPC.TR, T.cTR1=T.cTR, T.TR1=T.TR)   
T.cTR.vec= sapply(T.cTR.mat,c) #convert to vector
temps.cTR= c(temps.cTR, T.cTR.vec)

#WITHOUT TR NO EXTREMES
#unif
temps=lapply(hrs, FUN=noTR.noExt.unif, TeS.day1=TeS.day[1,mon,,Spec], Te.day1=Te.day[1,mon,,Spec])  
temps.mat= sapply(temps,c) #convert to matrix
temps.vec= sapply(temps.mat,c) #convert to vector
temps.nTR.noExt.unif=c(temps.nTR.noExt.unif, temps.vec ) 

#modal
temps=lapply(hrs, FUN=noTR.noExt.modal, TeS.day1=TeS.day[1,mon,,Spec], Te.day1=Te.day[1,mon,,Spec])  
temps.mat= sapply(temps,c) #convert to matrix
temps.vec= sapply(temps.mat,c) #convert to vector
temps.nTR.noExt.modal=c(temps.nTR.noExt.modal, temps.vec ) 

#bimodal
temps=lapply(hrs, FUN=noTR.noExt.bimodal, TeS.day1=TeS.day[1,mon,,Spec], Te.day1=Te.day[1,mon,,Spec])  
temps.mat= sapply(temps,c) #convert to matrix
temps.vec= sapply(temps.mat,c) #convert to vector
temps.nTR.noExt.bimodal=c(temps.nTR.noExt.bimodal, temps.vec ) 

}#end mon loop 

#FIT PERFORMANCE CURVE

Topts.shift= matrix(NA, nrow(Tm), 4) 

for(iter in 1:6){

if(iter==1) temp.dist=temps.nTR 
if(iter==2) temp.dist=temps.cTR
if(iter==3) temp.dist=temps.TR
if(iter==4) temp.dist=temps.nTR.noExt.unif
if(iter==5) temp.dist=temps.nTR.noExt.modal
if(iter==6) temp.dist=temps.nTR.noExt.bimodal

#=rnorm(1000, mean = 20, sd = 10) #to test

#bin into temperatures
breaks1= seq(-30.5, 60.5, 1)
breaks2= c(-50, breaks1, 80 ) 
h=hist(temp.dist, breaks= breaks2, plot=FALSE) 
 
temps= cbind(h$mids[2:length(h$mids)], h$counts[2:length(h$mids)])
temps= temps[temps[,2]>0,]

optim= comb[which.max(apply(comb, MARGIN=1, FUN=TPC.beta.sum, temps=temps, aran=0)),]
optim.prod= comb[which.max(apply(comb, MARGIN=1, FUN=TPC.beta.prod, temps=temps, aran=0)),]
Topts[cellk, iter,1:2,1]= as.numeric(optim)
Topts[cellk, iter,3:5,1]= TPC.beta.minmax( seq(-5,50, 0.5), as.numeric(optim[1]), as.numeric(optim[2]) )
Topts[cellk, iter,6:7,1]= as.numeric(optim.prod)
Topts[cellk, iter,8:10,1]=TPC.beta.minmax( seq(-5,50, 0.5), as.numeric(optim.prod[1]), as.numeric(optim.prod[2]))

#optimize with araneous
optim= comb[which.max(apply(comb, MARGIN=1, FUN=TPC.beta.sum, temps=temps, aran=1)),]
optim.prod= comb[which.max(apply(comb, MARGIN=1, FUN=TPC.beta.prod, temps=temps, aran=1)),]
Topts[cellk, iter,1:2,2]= as.numeric(optim)
Topts[cellk, iter,3:5,2]= TPC.beta.minmax(seq(-5,50, 0.5), as.numeric(optim[1]), as.numeric(optim[2]) )
Topts[cellk, iter,6:7,2]= as.numeric(optim.prod)
Topts[cellk, iter,8:10,2]=TPC.beta.minmax( seq(-5,50, 0.5), as.numeric(optim.prod[1]), as.numeric(optim.prod[2]))

#Store max Fitness
Fmax[cellk, iter, 1]= max(apply(comb, MARGIN=1, FUN=TPC.beta.sum, temps=temps, aran=0))
Fmax[cellk, iter, 2]= max(apply(comb, MARGIN=1, FUN=TPC.beta.sum, temps=temps, aran=1))

} #end iter
} # end cell loop


#---------------------------------------------------------
##SIMPLE PLOTS

plot(h, xlim=range(0,60))
TPCs= sapply(temps[,1], FUN=TPC.beta, shift=as.numeric(optim[1]), breadth=as.numeric(optim[2]), skew=0.7, aran=0)
TPCs.prod= sapply(temps[,1], FUN=TPC.beta, shift=as.numeric(optim.prod[1]), breadth=as.numeric(optim.prod[2]), skew=0.7, aran=0)
par(new=TRUE)
plot(temps[,1], TPCs, type="l", xlim=range(0,60))
points(temps[,1], TPCs.prod, type="l", lty="dashed")
 
#Simple map
library(maptools)
cols<-heat.colors(20)
cols<- cols[20:1]

par(mfrow=c(2,3), tcl=-0.5,mar=c(0,0,0,0), cex.lab=1.1, cex.axis=1, lwd=1, las=1, bty="l", pch=15, cex=1);

#plot data
data(wrld_simpl)

breaks.pl= quantile(Topts[, ,4,1], probs = seq(0, 1, length.out=20), na.rm=TRUE)
breaks.pl=unique(breaks.pl)

for(c in 1:6){
plot(wrld_simpl, xlim=c(-125,-65), ylim=c(30,42))
#split in quantiles
tb1<-cut(Topts[, c,4,1],breaks=breaks.pl) #col 4 is sum Topt and col 9 is prod Topt 
lev<-levels(tb1)
points(lon, lat, col=cols[match(tb1,lev)])
}


#---------------
