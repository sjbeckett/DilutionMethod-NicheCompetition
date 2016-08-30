# Part of the code used in:
# Beckett and Weitz. Code for: A revised dilution methodology and implications for estimates of rates of plankton mortality.
# 
# MIT License



TIMEINDEX = which(GG1_b_new$Timings==24) # 24 hours



METH = rbind(GG1_t_dil$Ind_appG[TIMEINDEX,,1],GG1_t_new$Ind_appG[TIMEINDEX,,1],GG1_i_dil$Ind_appG[TIMEINDEX,,1],GG1_i_new$Ind_appG[TIMEINDEX,,1],GG1_b_dil$Ind_appG[TIMEINDEX,,1],GG1_b_new$Ind_appG[TIMEINDEX,,1])*24
MORT = c(-truemort1_t,-truemort1_t,-truemort1_i,-truemort1_i,-truemort1_b,-truemort1_b)*24





GROWTHMORTLINES<-function(r1,MORT,rrr) {
	abline(h= r1*24,lwd=3,col='darkgreen')
	abline(h=MORT[rrr],lwd=3,col='red',lty=1)
}

LinRegPLOT <-function(FRACTIONS,METH,MORT,rrr) {
	QQQ=lm(METH[rrr,]~FRACTIONS)
	SEQ=c(0,1)
	LINMOD = SEQ*QQQ$coefficients[2]+QQQ$coefficients[1]
	lines(SEQ,LINMOD)


	xpos= 0.6	

	text(xpos,24*0.035,bquote(paste('m'['est']* ' = ', .(signif(QQQ$coefficients[2],3)))),pos=4) #estimated mortality rate
	text(xpos,21.5*0.035,bquote(paste('m'['act']* ' = ', .(signif(-MORT[rrr],3)))),pos=4) # actual mortality rate
	Err = 100*abs((-MORT[rrr] - QQQ$coefficients[2])/(-MORT[rrr])) # percentage error in grazing mortality rate estimate
	text(xpos,19*0.035,bquote(paste('error = ', .(signif(Err,3)), '%')),pos=4)

}


FRACTIONS = seq(0.1,1,0.1)
XLIMS = c(-0.01,1.015)
YLIMS =  c(-0.001,r1*1.05*24)

par(mfrow=c(3,2), mar=c(4,5,0.5,0.5), oma=c(0.2,4,4.5,0.2)) # mar: bottom, left , top of all plots, space between plots

rrr=1
plot(FRACTIONS,METH[rrr,],type='p',xlim=XLIMS,ylim=YLIMS,xaxs='i',yaxs='i',xlab="",ylab=expression("App. growth rate (d"^-1 *")"),cex.lab=1.5,pch=20,cex=2)
GROWTHMORTLINES(r1,MORT,rrr)
LinRegPLOT(FRACTIONS,METH,MORT,rrr)

rrr=2
plot(FRACTIONS,METH[rrr,],type='p',xlim=XLIMS,ylim=YLIMS,xaxs='i',yaxs='i',xlab="",ylab="",pch=20,cex=2)
GROWTHMORTLINES(r1,MORT,rrr)
LinRegPLOT(FRACTIONS,METH,MORT,rrr)

rrr=3
plot(FRACTIONS,METH[rrr,],type='p',xlim=XLIMS,ylim=YLIMS,xaxs='i',yaxs='i',xlab="",ylab=expression("App. growth rate (d"^-1 *")"),cex.lab=1.5,pch=20,cex=2)
GROWTHMORTLINES(r1,MORT,rrr)
LinRegPLOT(FRACTIONS,METH,MORT,rrr)

rrr=4
plot(FRACTIONS,METH[rrr,],type='p',xlim=XLIMS,ylim=YLIMS,xaxs='i',yaxs='i',xlab="",ylab="",pch=20,cex=2)
GROWTHMORTLINES(r1,MORT,rrr)
LinRegPLOT(FRACTIONS,METH,MORT,rrr)

rrr=5
plot(FRACTIONS,METH[rrr,],type='p',xlim=XLIMS,ylim=YLIMS,xaxs='i',yaxs='i',xlab="Proportion WSW",ylab=expression("App. growth rate (d"^-1 *")"),cex.lab=1.5,pch=20,cex=2)
GROWTHMORTLINES(r1,MORT,rrr)
LinRegPLOT(FRACTIONS,METH,MORT,rrr)

rrr=6
plot(FRACTIONS,METH[rrr,],type='p',xlim=XLIMS,ylim=YLIMS,xaxs='i',yaxs='i',xlab="Proportion WSW",ylab="",cex.lab=1.5,pch=20,cex=2)
GROWTHMORTLINES(r1,MORT,rrr)
LinRegPLOT(FRACTIONS,METH,MORT,rrr)


par(las=1)
mtext("Classical", NORTH<-3, line=0, adj=0.25, cex=1.2, col="black", outer=TRUE)
mtext("Z-dilution", NORTH<-3, line=0, adj=0.85, cex=1.2, col="black", outer=TRUE)
#par(las=2)#yes but also no.. make it hard!
mtext("HP",line=0,side=3,outer=TRUE,cex=1.2,col='black', adj = -0.07,padj=7)
mtext("IP",line=0,side=3,outer=TRUE,cex=1.2,col='black', adj = -0.07,padj=21.5)
mtext("LP",line=0,side=3,outer=TRUE,cex=1.2,col='black', adj = -0.07,padj=36)






par(fig=c(0.2,1,0,1),oma=c(0,0,0,0),mar=c(0,0,0,0),new=TRUE)
plot(0,0,type="n",bty="n",xaxt="n",yaxt="n")
legend("top",c("Plankton growth rate","Baseline grazing mortality rate","Apparent growth rate estimates","Linear regression fit to data"),lwd=c(2,2,2,1),col=c("darkgreen","red","black","black"),lty=c(1,1,NA,1),pch = c(NA,NA,20,NA),bty='n',ncol=2)



dev.copy2eps(file="../../figures/ZOOAPPGROWTH.eps")
dev.copy2pdf(file="../../figures/ZOOAPPGROWTH.pdf")
dev.off()

