# Part of the code used in:
# Beckett and Weitz. Code for: A revised dilution methodology and implications for estimates of rates of plankton mortality.
# 
# MIT License

	#load required files
	source("../main/LOADREQUIREMENTS.R") # read in files and packages for dilution experiments.


	#set up basic model parameters

	r1 = 1/24   #plankton growth rate (per hour)
	K1 = 2.2*10^7   #plankton carrying capcity   (per mililitre)
	a1 = 2*10^-6  # attack rate of zooplankton on plankton (ml per zooplankton per hour)
	Dil_levels = seq(0.1,1,0.1) # dilution levels to use -- proportion of whole sea water (WSW in each bottle).
	
	# parameter vector
	ng_params = c(growth= r1, Amplitude = 0 , PERIODICITY = 24, timeshift = 0, carryingCap = K1, GrazingRate = a1)


	#logistic plankton growth - dynamic grazers
	PROJECTINFO = list(filename="Ngbar.RData",modelused="ng")
	

	#CONSTANT ZOOPLANKTON MODEL

	#1)  Low top-down control
	G0=1000 # Top-down control is low
	LowPressure = a1*G0/r1

	#Equilibrium dynamics for low top-down control
	ng_state = list(c( N = K1*(1-(a1*G0)/r1), G = G0 ))[[1]] # steady state equilibrim dynamics
	
	#run dilution experiments
	GG1_b_new = ANALYSIS(NGlogmodel,ng_params,ng_state, host_vector=c(1), dilute_vector=c(2), dilution_levels = Dil_levels, saveon=0, PROJECTINFO) # revised dilution (only zooplankton subjected to dilution -- plankton fit through filter) 
	GG1_b_dil = ANALYSIS(NGlogmodel,ng_params,ng_state, host_vector=c(1), dilute_vector=c(1,2), dilution_levels = Dil_levels, saveon=0, PROJECTINFO) # classic dilution (both plankton and zooplankton subject to dilution)

	truemort1_b   =  -a1*G0  # true mortality rate at steady state dynamics.


	#2)  Intermediate top-down control
	G0=10000
	IntermediatePressure = a1*G0/r1

	#Equilibrium dynamics for intermediate top-down control
	ng_state = list(c( N = K1*(1-(a1*G0)/r1), G = G0 ))[[1]] # steady state equilibrim dynamics

	#run dilution experiments
	GG1_i_new = ANALYSIS(NGlogmodel,ng_params,ng_state, host_vector=c(1), dilute_vector=c(2), dilution_levels = Dil_levels, saveon=0, PROJECTINFO) # revised dilution (only zooplankton subjected to dilution -- plankton fit through filter)
	GG1_i_dil = ANALYSIS(NGlogmodel,ng_params,ng_state, host_vector=c(1), dilute_vector=c(1,2), dilution_levels = Dil_levels, saveon=0, PROJECTINFO) # classic dilution (both plankton and zooplankton subject to dilution)

	truemort1_i   =  -a1*G0 # true mortality rate at steady state dynamics.


	#3)  High top-down control
	G0=20000 #top
	HighPressure = a1*G0/r1

	#Equilibrium dynamics for high top-down control
	ng_state = list(c( N = K1*(1-(a1*G0)/r1), G = G0 ))[[1]]

	#run dilution experiments
	GG1_t_new = ANALYSIS(NGlogmodel,ng_params,ng_state, host_vector=c(1), dilute_vector=c(2), dilution_levels = Dil_levels, saveon=0, PROJECTINFO) # revised dilution (only zooplankton subjected to dilution -- plankton fit through filter)
	GG1_t_dil = ANALYSIS(NGlogmodel,ng_params,ng_state, host_vector=c(1), dilute_vector=c(1,2), dilution_levels = Dil_levels, saveon=0, PROJECTINFO) # classic dilution (both plankton and zooplankton subject to dilution)

	truemort1_t   =  -a1*G0 # true mortality rate at steady state dynamics.


##### PLOTS #######

# Create plots showing the performance of the classical and revised dilution methods to infer zooplankton associated mortality rate for each of the three levels of top-down control:


#Classical dilution method at three example levels of zooplankton pressure


	
	plot(0,0,xaxt="n",yaxt="n",bty="n",xlab="",ylab="",col="white")

	par(mar=c(4,4.8,2,1))
	par(oma=c(0,0,0,0))

	INDEX24H = which(GG1_b_new$Timings==24)

	plot(NA,NA,xlim=c(0.5,3.5),ylim=c(0,24*0.045),xaxt='n',yaxt='n',xaxs='i',yaxs='i',xlab=expression(" Grazing pressure (" *delta [Z] *") "),ylab=expression("Mortality rate (d"^-1*")"),cex.lab=1.8)
	axis(2,cex.axis=1.8) #yaxis
	axis(1,at=1:3,labels=c("Low","Intermediate","High"),cex.axis=1.8)

	width=0.3
	HPRANGE=c(3-width,3+width)
	IPRANGE=c(2-width,2+width)
	LPRANGE=c(1-width,1+width)

	#DRAWEDGES
	BOXCOL = 'grey'
	SSCOL = 'black'
	DILPZCOL = 'brown'
	DILZCOL = 'blue'
	MAXCOL = 'red'

	#Convert from per hour to per day

	TOPd = -GG1_t_dil$All_DY*24
	MIDd = -GG1_i_dil$All_DY*24
	BOTd = -GG1_b_dil$All_DY*24

	baseline_t = -truemort1_t*24
	baseline_i = -truemort1_i*24
	baseline_b = -truemort1_b*24

	RR = r1*24

	
	rect(HPRANGE[1],baseline_t,HPRANGE[2],RR,col="grey48",border='NA')
	
	rect(IPRANGE[1],baseline_i,IPRANGE[2],RR,col=BOXCOL,border='NA')

	rect(LPRANGE[1],baseline_b,LPRANGE[2],RR,col="grey87",border='NA')

	#DRAW STEADY STATES
	CORRECTION = 0.005
	LWDhere=8
	lines(HPRANGE+CORRECTION,rep(baseline_t,2),lwd=LWDhere,lty=3,col=SSCOL)
	lines(IPRANGE+CORRECTION,rep(baseline_i,2),lwd=LWDhere,lty=3,col=SSCOL)
	lines(LPRANGE+CORRECTION,rep(baseline_b,2),lwd=LWDhere,lty=3,col=SSCOL)
	#DRAW DILUTION METHOD @24h
	lines(HPRANGE,rep(TOPd[INDEX24H],2),col=DILPZCOL,lwd=LWDhere)
	lines(IPRANGE,rep(MIDd[INDEX24H],2),col=DILPZCOL,lwd=LWDhere)
	lines(LPRANGE,rep(BOTd[INDEX24H],2),col=DILPZCOL,lwd=LWDhere)
	
	#DRAW WORST CASE (growth rate)
	lines(HPRANGE,rep(RR,2),col=MAXCOL,lwd=LWDhere)
	lines(IPRANGE,rep(RR,2),col=MAXCOL,lwd=LWDhere)
	lines(LPRANGE,rep(RR,2),col=MAXCOL,lwd=LWDhere)

	#DRAWPOINTS
	points(HPRANGE,rep(TOPd[INDEX24H],2),col=DILPZCOL,cex=2,lwd=8)
	points(IPRANGE,rep(MIDd[INDEX24H],2),col=DILPZCOL,cex=2,lwd=8)
	points(LPRANGE,rep(BOTd[INDEX24H],2),col=DILPZCOL,cex=2,lwd=8)


	legend(1.53,24*0.013,c("Baseline rate","Classical estimate",expression(paste("Maximum rate ( ",italic("m=r")," )",sep=""))),col=c(SSCOL,DILPZCOL,MAXCOL),pch=c(NA,1,NA),lty=c(3,1,1),lwd=c(8,8,8),bty='n',pt.cex=2,cex=1.7)

	dev.copy2eps(file="../../figures/classicaldilution_mortalityrates.eps")
	dev.copy2pdf(file="../../figures/classicaldilution_mortalityrates.pdf")
	dev.off()



#Both classical and Z-dilution method at three example levels of zooplankton pressure

	plot(0,0,xaxt="n",yaxt="n",bty="n",xlab="",ylab="",col="white")

	par(mar=c(4,4.8,2,1))
	par(oma=c(0,0,0,0))

	INDEX24H = which(GG1_b_new$Timings==24)

	plot(NA,NA,xlim=c(0.5,3.5),ylim=c(0,24*0.045),xaxt='n',yaxt='n',xaxs='i',yaxs='i',xlab=expression(" Grazing pressure (" *delta [Z] *") "),ylab=expression("Mortality rate (d"^-1*")"),cex.lab=1.8)
	axis(2,cex.axis=1.8) #yaxis
	axis(1,at=1:3,labels=c("Low","Intermediate","High"),cex.axis=1.8)

	width=0.3
	HPRANGE=c(3-width,3+width)
	IPRANGE=c(2-width,2+width)
	LPRANGE=c(1-width,1+width)

	#DRAWEDGES
	BOXCOL = 'grey'
	SSCOL = 'black'
	DILPZCOL = 'brown'
	DILZCOL = 'blue'
	MAXCOL = 'red'

	#Convert from per hour to per day

	TOP = -GG1_t_new$All_DY*24
	MID = -GG1_i_new$All_DY*24
	BOT = -GG1_b_new$All_DY*24

	TOPd = -GG1_t_dil$All_DY*24
	MIDd = -GG1_i_dil$All_DY*24
	BOTd = -GG1_b_dil$All_DY*24

	baseline_t = -truemort1_t*24
	baseline_i = -truemort1_i*24
	baseline_b = -truemort1_b*24

	RR = r1*24


	#lines(rep(HPRANGE[1],2),c(-GG1_t_new$All_DY[INDEX24H],r1))
	#lines(rep(HPRANGE[2],2),c(-GG1_t_new$All_DY[INDEX24H],r1))
	rect(HPRANGE[1],baseline_t,HPRANGE[2],RR,col="grey48",border='NA')
	
	#lines(rep(IPRANGE[1],2),c(-GG1_i_new$All_DY[INDEX24H],r1))
	#lines(rep(IPRANGE[2],2),c(-GG1_i_new$All_DY[INDEX24H],r1))
	rect(IPRANGE[1],baseline_i,IPRANGE[2],RR,col=BOXCOL,border='NA')

	#lines(rep(LPRANGE[1],2),c(-GG1_b_new$All_DY[INDEX24H],r1))
	#lines(rep(LPRANGE[2],2),c(-GG1_b_new$All_DY[INDEX24H],r1))
	rect(LPRANGE[1],baseline_b,LPRANGE[2],RR,col="grey87",border='NA')

	#DRAW STEADY STATES
	CORRECTION = 0.005
	LWDhere=8
	lines(HPRANGE+CORRECTION,rep(baseline_t,2),lwd=LWDhere,lty=3,col=SSCOL)
	lines(IPRANGE+CORRECTION,rep(baseline_i,2),lwd=LWDhere,lty=3,col=SSCOL)
	lines(LPRANGE+CORRECTION,rep(baseline_b,2),lwd=LWDhere,lty=3,col=SSCOL)
	#DRAW DIL @24h
	lines(HPRANGE,rep(TOPd[INDEX24H],2),col=DILPZCOL,lwd=LWDhere)
	lines(IPRANGE,rep(MIDd[INDEX24H],2),col=DILPZCOL,lwd=LWDhere)
	lines(LPRANGE,rep(BOTd[INDEX24H],2),col=DILPZCOL,lwd=LWDhere)
	
	#DRAW PRESSURE DIL @24h
	lines(HPRANGE,rep(TOP[INDEX24H],2),col=DILZCOL,lwd=LWDhere)
	lines(IPRANGE,rep(MID[INDEX24H],2),col=DILZCOL,lwd=LWDhere)
	lines(LPRANGE,rep(BOT[INDEX24H],2),col=DILZCOL,lwd=LWDhere)
	#DRAW WORST CASE (growth rate)
	lines(HPRANGE,rep(RR,2),col=MAXCOL,lwd=LWDhere)
	lines(IPRANGE,rep(RR,2),col=MAXCOL,lwd=LWDhere)
	lines(LPRANGE,rep(RR,2),col=MAXCOL,lwd=LWDhere)

	#DRAWPOINTS
	points(HPRANGE,rep(TOPd[INDEX24H],2),col=DILPZCOL,cex=2,lwd=8)
	points(IPRANGE,rep(MIDd[INDEX24H],2),col=DILPZCOL,cex=2,lwd=8)
	points(LPRANGE,rep(BOTd[INDEX24H],2),col=DILPZCOL,cex=2,lwd=8)

	points(HPRANGE,rep(TOP[INDEX24H],2),col=DILZCOL,pch=2,cex=2,lwd=8)
	points(IPRANGE,rep(MID[INDEX24H],2),col=DILZCOL,pch=2,cex=2,lwd=8)
	points(LPRANGE,rep(BOT[INDEX24H],2),col=DILZCOL,pch=2,cex=2,lwd=8)

	legend(1.53,24*0.013,c("Baseline rate","Classical estimate","Z-dilution estimate",expression(paste("Maximum rate ( ",italic("m=r")," )",sep=""))),col=c(SSCOL,DILPZCOL,DILZCOL,MAXCOL),pch=c(NA,1,2,NA),lty=c(3,1,1,1),lwd=c(8,8,8,8),bty='n',pt.cex=2,cex=1.7)


	dev.copy2eps(file="../../figures/classical_Zdilution_mortalityrates.eps")
	dev.copy2pdf(file="../../figures/classical_Zdilution_mortalityrates.pdf")
	dev.off()


