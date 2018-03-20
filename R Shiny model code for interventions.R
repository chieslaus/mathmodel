#MALARIA MODEL R CODE 18 MARCH 2018
version # version.string R version 3.4.2 (2017-09-28) #nickname       Short Summer
#RSHINY INTERVENTION

#CLEAR WORKSPACE
rm(list=ls())

#LOAD REQUIRED LIBRARIES
library(deSolve)
library(ggplot2)
library(shiny)

# Define UI for application that draws a histogram
ui <-  fluidPage(
  
  tabsetPanel(
    id="panels",
    tabPanel(title = strong("Transmission Model"),
             titlePanel("Typology Parameters"),
             
             sidebarPanel(            
               numericInput(inputId="initP", label = "population size", value=186000000, min=1000, max=200000000, step=100),
               sliderInput(inputId="lifeexp", label = "average life expectancy (years)", value = 54, min=40, max=100,step=1),
               sliderInput(inputId="R0", label = "basic reproduction number", value = 4.4, min=0, max=10,step=0.1),
               sliderInput(inputId="dur_imm", label = "duration of immunity (years)", value = 5, min=0, max=20,step=0.1),
               sliderInput(inputId="dur_latency", label = "duration of latency (weeks)", value = 2, min=1, max=5,step=0.1),
               sliderInput(inputId="dur_infect", label = "duration of infectiousness (weeks)", value = 2, min=2, max=8,step=0.1),
               sliderInput(inputId="report", label = "proportion of all infections that are reported", value = 16/100, min=0, max=1,step=0.01),
               sliderInput(inputId="amp", label = "relative amplitude of seasonal forcing", value = 0, min=0, max=1,step=0.01),
               sliderInput(inputId="phi", label = "week of peak in seasonal forcing", value = 30, min=0, max=52,step=1)
               
             ),
             # Show a plot model output
             mainPanel(
               plotOutput(outputId = "tab1") 
             ))
    ,
    tabPanel(title = strong("Intervention 1"),
             titlePanel("Increased bednet coverage"),
             sidebarPanel(            
               sliderInput(inputId="week_net" , label = "weeks after first case when the intervention starts", value=1, min=0, max=72, step=1),
               sliderInput(inputId="net_cov_i", label = "intervention coverage of bednet", value=0.8, min=0, max=1, step=0.01),
               sliderInput(inputId="net_eff", label = "efficacy of bednet in reducing transmission", value=0.5, min=0, max=1, step=0.01)
             ),
             #Show a plot model output
             mainPanel(
               plotOutput(outputId = "tab2")
             )
    ),
    tabPanel(title = strong("Intervention 1"),
             titlePanel("Vaccination"),
             sidebarPanel(            
               sliderInput(inputId="vac_cov_i" , label = "coverage of vaccine campaign", value=0.7, min=0, max=1, step=0.01),
               sliderInput(inputId="campaignweeks", label = "time taken to reach target coverage", value=4, min=1, max=52, step=0.01)
             ),
             #Show a plot model output
             mainPanel(
               plotOutput(outputId = "tab3")
             )
    )
  ))

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
                phi=30               # week of peak in seasonal forcing
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


# Define server logic required to run model and plot outputs
server <- function(input, output) {
  #Create parameters vector - reactive nature
  parametersR <- reactive(c(mui=(1/(input$lifeexp*52)),    # birth
                            muo=(1/(input$lifeexp*52)),    # death
                            R0=input$R0,               # basic reproduction number
                            omega=(1/(input$dur_imm*52)),  # rate of loss of immunity = 1/(average duration of immunity)
                            gamma=1/input$dur_latency,  # rate of movement from latent to infectious stage = 1/(average latent period)
                            nui=1/input$dur_infect,  # rate of recovery = 1/(average duration of infection)
                            report=input$report,         # proportion of all infections that are reported
                            amp=input$amp,              # relative amplitude of seasonal forcing
                            phi=input$phi,               # week of peak in seasonal forcing
                            week_interv = input$week_interv,   # weeks after first case when the intervention starts
                            lat_cov_i = input$lat_cov_i,   # intervention coverage of latrines
                            lat_eff = input$lat_eff,       # efficacy of latrines in reducing transmission                            
                            vac_cov_i = input$vac_cov_i,     # coverage of vaccine campaign
                            campaignweeks = input$campaignweeks   # time taken to reach target coverage))
  ))
  
  #Define intital values  
  initP <- reactive(input$initP)
  initE<-1 # Exposed
  initI<-0 # Infectious
  initR<-0 # Immune
  
  #Create reactive plot function
  output$tab1<-output$tab2<-output$tab3<-output$Plot <- renderPlot({
    
    initS=initP()-initE-initI-initR # Susceptible (non-immune)
    state <- c(S = initS, E=initE, I = initI,R = initR)
    
    #Run model
    outR <- reactive(ode(y = state, times = times, func = malaria, parms = parametersR()))  #reactive object
    out<-outR()  

        


# MORE MODEL OUTPUTS
# total population
pop<-out[,"S"]+out[,"E"]+out[,"I"]+out[,"R"]
# weekly incidence
inc <- parameters["report"]*parameters["gamma"]*out[,"E"]

time <-out[,"time"]



# make a new panel  # label the axes 
par(mfrow=c(1,2))
plot(time,inc,type='l',lwd=3,main = "Malaria Epidemic",xlab = "Time in weeks",ylab="New reported cases per week")

# plot cumulative reports
cumulativereports<- parametersR()["report"]*out[,"R"]
plot(times,cumulativereports,type='l',lwd=3,main = "Malaria Epidemic",xlab = "Time in weeks",ylab="Cumulative cases") 
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

