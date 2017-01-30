# Part of the code used in:
# Beckett and Weitz. Code for: A revised dilution methodology and implications for estimates of rates of plankton mortality.
# 
# MIT License


	#load required files
	source("../main/LOADREQUIREMENTS.R") # read in files and packages for dilution experiments.


	#set up basic model parameters

	r1 = 1/24   #plankton growth rate (per hour)
	K1 = 2.2*10^7   #plankton carrying capcity   (per mililitre)
	a1 = 2*10^-6  # attack rate of zooplankton on plankton (ml per viruses per hour)
	Dil_levels = seq(0.1,1,0.1) # dilution levels to use -- proportion of whole sea water (WSW in each bottle).

	DL = 100 # number of pressure levels
	delta = seq(0,1,length.out=DL)  #pressure levels
	G0 = (r1/a1)*delta   # Zooplankton population size associated with each of the DL pressure levels
	#delta = 0  no zooplankton (low top-down control), delta = 1 -- high zooplankton (high top-down control)
	
	#make storage placeholders for mortality rates through time for each method
	ALL_DY_new = matrix(0,DL,721)  # B&W - revised
	ALL_DY_dil = matrix(0,DL,721)  # L&H - classical

	BIAS_new = matrix(0,DL,721)
	BIAS_dil =  matrix(0,DL,721)
	truemort1_bVIR = c()


	#Run the simulations and collect mortality and bias information.
	for( aa in 1:(DL-1)) { #The last delta is the case when Z = 0 and P=K; and is not calculated as no dynamics occur - causing the linear regression to fail.


		ng_params = c(growth= r1, Amplitude = 0 , PERIODICITY = 24, timeshift = 0, carryingCap = K1, GrazingRate = a1)
		
		PROJECTINFO = list(filename="Ng.RData",modelused="ng")
	

		ng_state = list(c( N = K1*(1-(a1*G0[aa])/r1), G = G0[aa] ))[[1]]
	

		GG1_b_new = ANALYSIS(NGlogmodel,ng_params,ng_state,c(1),c(2),Dil_levels,0,PROJECTINFO) # revised dilution method -- only Z's filtered out
		GG1_b_dil = ANALYSIS(NGlogmodel,ng_params,ng_state,c(1),c(1,2),Dil_levels,0,PROJECTINFO) #classic dilution theory -- 0.45uM filter
		

		ALL_DY_new[aa,]=GG1_b_new$All_DY
		ALL_DY_dil[aa,]=GG1_b_dil$All_DY

		truemort1_b   =  -a1*G0[aa]


		BIAS_new[aa,] = GG1_b_new$All_DY/ (-a1*ng_state[2]) # bias with revised method
		BIAS_dil[aa,] = GG1_b_dil$All_DY/ (-a1*ng_state[2]) # bias with classical method
	
	print(aa)
	}

	#save this data so it doesn't need to be run everytime want to create the plots.
	save(list=ls(all.names=TRUE),file="Zooplankton_ACROSS_CONTROL.RData")


#### PLOTS ######

#classical dilution mortality rate bias

	par(mar=c(4,4.8,2,1))
	par(oma=c(0,0,0,0))

	
	OFFSET = 2*10^-3
	TIMINGS=24
	rrr=1

	YLIM=c(0.1,100)
	TIMEINDEX = which(GG1_b_new$Timings==TIMINGS[rrr])  ##
	plot(NA, NA , type='l',  xlim=c(0,1) ,  xlab=expression(" Grazing pressure (" *delta [Z] *") "), ylab = "Mortality rate bias", log="y", ylim=YLIM, xaxs="i", yaxs="i",cex.lab=1.8,cex.axis=1.8)

	rect(OFFSET,0.75,1-OFFSET,1.25,col='grey50',border='grey50')
	rect(OFFSET,0.9,1-OFFSET,1.1,col='grey',border='grey')

	LWD =5
	lines(delta, BIAS_dil[,TIMEINDEX] ,col='brown' ,lwd=LWD)
	POINTPOS = seq(2,99,8)
	points(delta[POINTPOS],BIAS_dil[POINTPOS,TIMEINDEX] ,col='brown',pch=1,lwd=LWD,cex=1.5)
	

	legend(0.38,50,c("Classical estimate","±10% baseline rate","±25% baseline rate"),pch=c(1,NA,NA),col=c("brown","grey","grey50"),lty=c(1,1,1),bty='n',lwd=c(LWD,8*2,8*2),pt.cex=1.8,cex=1.8)
	
	YPOStxt = YLIM[1]*1.1
	TXTsize = 1.8
	text(0.048,YPOStxt,"Low",pos=3,cex=TXTsize)
	text(0.48,YPOStxt,"Intermediate",pos=3,cex=TXTsize)
	text(0.96-0.02,YPOStxt,"High",pos=3,cex=TXTsize)
	YMAXpos = YLIM[1]*1.1
	lines(rep(0.048,2),c(YLIM[1],YMAXpos))
	lines(rep(0.48,2),c(YLIM[1],YMAXpos))
	lines(rep(0.96,2),c(YLIM[1],YMAXpos))

	dev.copy2eps(file="../../figures/classicaldilution_mort_bias.eps")
	dev.copy2pdf(file="../../figures/classicaldilution_mort_bias.pdf")
	dev.off()


#classical and Z-dilution mortality rate bias

	par(mar=c(4,4.8,2,1))
	par(oma=c(0,0,0,0))

	plot(0,0,xaxt="n",yaxt="n",bty="n",xlab="",ylab="",col="white")
	TIMINGS=24
	rrr=1

	YLIM=c(0.1,100)
	TIMEINDEX = which(GG1_b_new$Timings==TIMINGS[rrr])  ##
	plot(NA, NA , type='l',  xlim=c(0,1) ,  xlab=expression(" Grazing pressure (" *delta [Z] *") "), ylab = "Mortality rate bias", log="y", ylim=YLIM, xaxs="i", yaxs="i",cex.lab=1.8,cex.axis=1.8) 

	rect(OFFSET,0.75,1-OFFSET,1.25,col='grey50',border='grey50')
	rect(OFFSET,0.9,1-OFFSET,1.1,col='grey',border='grey')

	LWD =5
	lines(delta, BIAS_dil[,TIMEINDEX] ,col='brown' ,lwd=LWD)
	POINTPOS = seq(2,99,8)
	points(delta[POINTPOS],BIAS_dil[POINTPOS,TIMEINDEX] ,col='brown',pch=1,lwd=LWD,cex=1.5)
	lines(delta, BIAS_new[,TIMEINDEX],col='blue',lwd=LWD)
	points(delta[POINTPOS],BIAS_new[POINTPOS,TIMEINDEX] ,col='blue',pch=2,lwd=LWD,cex=1.5)
	

	legend(0.38,50,c("Classical estimate","Z-dilution estimate","±10% true rate","±25% true rate"),pch=c(1,2,NA,NA),col=c("brown","blue","grey","grey50"),lty=c(1,1,1,1),bty='n',lwd=c(LWD,LWD,8*2,8*2),pt.cex=1.8,cex=1.8)


	YPOStxt = YLIM[1]*1.1
	TXTsize = 1.8
	text(0.048,YPOStxt,"Low",pos=3,cex=TXTsize)
	text(0.48,YPOStxt,"Intermediate",pos=3,cex=TXTsize)
	text(0.96-0.02,YPOStxt,"High",pos=3,cex=TXTsize)
	YMAXpos = YLIM[1]*1.1
	lines(rep(0.048,2),c(YLIM[1],YMAXpos))
	lines(rep(0.48,2),c(YLIM[1],YMAXpos))
	lines(rep(0.96,2),c(YLIM[1],YMAXpos))

	dev.copy2eps(file="../../figures/classical_Zdilution_mort_bias.eps")
	dev.copy2pdf(file="../../figures/classical_Zdilution_mort_bias.pdf")
	dev.off()



