#MALARIA MODEL R CODE 18 MARCH 2018
version # version.string R version 3.4.2 (2017-09-28) #nickname       Short Summer

# 1. BASIC MODEL

#CLEAR WORKSPACE
rm(list=ls())

#LOAD REQUIRED LIBRARIES
library(deSolve)
library(ggplot2)

#MODEL TIME PERIOD = 10 (in years)
period <- 520 #weeks
times <- seq(0, period, by = (1/7)) # in weeks, by 1 day intervals

#MODEL PARAMETERS
parameters <- c(mui=(1/(54*52)),  # birth rate (1/life expectancy)
                muo=(1/(54*52)),  # death rate (1/life expectancy)
                R0=4.4,              # basic reproduction number
                omega=(1/(5*52)),    # rate of loss of immunity = 1/(average duration of immunity)
                gamma=7/10,          # rate of movement from latent to infectious stage = 1/(average latent period)
                nui=35*7/1000,       # rate of recovery = 1/(average duration of infection)
                report=16/100,       # proportion of all infections that are reported
                amp=0,               # relative amplitude of seasonal forcing
                phi=30,               # week of peak in seasonal forcing
                
                week_net = 1,        #week of commencing intervention (ITNs)
                net_cov = .80,        #ITN coverage
                net_eff = .50        #ITN efficacy in reducing malaria
)

# MODEL INITIAL CONDITIONS
initP<-186000000 # population size
initE<-1 # Exposed
initI<-0 # Infectious
initR<-0 # Recovered

#INITIAL SUSCEPTIBLE COMPARTMENT
initS<-initP-initE-initI-initR # Susceptible (non-immune)

#DEFINE THE MODEL COMPARTMENTS
state <- c(S = initS, E=initE, I = initI,R = initR)

# MODEL FUNCTION, TO SOLVE THE EQUATIONS
malaria<-function(t, state, parameters) 
{
  with(as.list(c(state, parameters)),
       {
         
         # DEFINE VARIABLES
         P <- (S+E+I+R)
         seas<-1+amp*cos(2*pi*(t-phi)/52)
         beta<-R0*(muo+nui)*(gamma+muo)/gamma
         
         net_cov <- (t>=(1 +week_net))*net_cov
         bednet <- (1-net_eff*net_cov)
         
         lam <- beta*seas*I/P
         
         # RATE OF CHANGE
         dS <- mui*P-muo*S-lam*S+omega*R
         dE <- -muo*E+lam*S-gamma*E
         dI <- -muo*I+gamma*E-nui*I
         dR <- -muo*R+nui*I-omega*R
         
         # RETURN THE RATE OF CHANGE
         list(c(dS, dE, dI, dR))
       }
  ) 
  
}

#MODEL OUTPUT
# run the model
out <- ode(y = state, times = times, func = malaria, parms = parameters)
# a simple plot of the model output
plot(out)

# MORE MODEL OUTPUTS
# total population
pop<-out[,"S"]+out[,"E"]+out[,"I"]+out[,"R"]
# weekly incidence
inc <- parameters["report"]*parameters["gamma"]*out[,"E"]


time <-out[,"time"]


# label the axes # make a new panel
par(mfrow=c(1,2))
plot(time,inc,type='l',lwd=3,main = "Malaria Epidemic",xlab = "Time in weeks",ylab="New reported cases per week")

# plot cumulative reports
cumulativereports<- parameters["report"]*out[,"R"]
plot(time,cumulativereports,type='l',lwd=3,main = "Malaria Epidemic",xlab = "Time in weeks",ylab="Cumulative cases") 
