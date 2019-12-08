#bios 619 final project multiple myeloma smulation
library(dfcrm)
library(clinfun)

set.seed(619)
#phase i (toxicity)
# CRM non-TITE 28 day tox period(?) so it is ok
# need: Probability of selecting each dose as the MTD
#Number/proportion of patients given each dose
#Number of DLTs per dose
#Expected sample size
#Expected study duration

#beta prior i guess for now
#1 param logistic  
prior_tox = c(5,10,15,19)/100
pi_nmax = 35 #max enrollment for phase I (max 4-6 subjects/month)
dlt = 0.2 #target

#Phase I
sim1 = crmsim( PI = c(5,10,15,19)/100, #true probabilities (can change this)
               prior = c(5,10,15,19)/100, #skeleton 
               target = 0.2, #target DLT rate, standard is 1/3
               n = 35, #trial sample size
               x0 = 1, #initial dose
               nsim = 1000,
               mcohort = 5, #number per cohort
               restrict = TRUE, #purely model based is FALSE option, this is no skipping and no esc after toxic
               count = FALSE,
               method = "bayes",
               model = "logistic",
               scale = sqrt(1.34), # default
               seed = 619
)

#results selected dose level 4, 19%

#if time, run competing (aka 3+3 or biased coin)

#goal is to treat 36 people at the MTD
#phase ii simulation: (efficacy)
#Simon two stage, optimuum (for myeloma) 
#25% response no-go, 45% promising
sim2 = ph2simon(pu = 0.25, pa = 0.45, ep1 = 0.05, ep2 = 0.20, nmax = 100)

#

