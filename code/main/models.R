# Part of the code used in:
# Beckett and Weitz. Code for: A revised dilution methodology and implications for estimates of rates of plankton mortality.
# 
# MIT License

#Dynamical model for anlaysis

NGlogmodel <- function(t,state,params) { #logistic growth model with grazing 
		with(as.list(c(state, params)), {
		dNdt <-  growth*N*(1-(N/carryingCap)) - GrazingRate*G*N
		dGdt <-  0
		list(c(dNdt,dGdt))
		})
	}
