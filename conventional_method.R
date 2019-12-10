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

#Phase I
sim1 <- crmsim( PI = c(5,10,15,19)/100, #true probabilities (can change this)
               prior = c(5,10,15,19)/100, #skeleton 
               target = 0.2, #target DLT rate, standard is 1/3
               n = 35, #trial sample size
               x0 = 1, #initial dose
               nsim = 1000, #change this to 1000, takes a while to run is why it's set to 10
               mcohort = 5, #number per cohort
               restrict = TRUE, #purely model based is FALSE option, this is no skipping and no esc after toxic
               count = FALSE,
               method = "bayes",
               model = "logistic",
               scale = sqrt(1.34), # default
               seed = 619
)

sim1
#from sim i:
#tox level 4 is selected most often

#goal is to treat 36 people at the MTD
#phase ii simulation: (efficacy)
#Simon two stage, optimuum (for myeloma) 
#25% response no-go, 45% promising
sim2 <- ph2simon(pu = 0.25, pa = 0.45, ep1 = 0.05, ep2 = 0.10, nmax = 100)

#simulating stage 2 outcomes
nsims2 <- 1000
pi_tox <- 0.19
pi_eff_mono <- 0.4 #from STEIN sim of prior efficacy
pi_eff_lvl <- 0.5
pi_eff_quad <- 0.25
n1 <- 22 #stage 1 ss
n <- 57
r1 <- 6
r <- 19

stage2sim <- function(nsims,r1,r,n1,n,prior_tox,prior_eff) {
effs <- rep(0,nsims)
toxs <- rep(0,nsims)
ns <- rep(0,nsims)
result <- rep(NA,nsims2)
for (i in 1:nsims) {
  temp_tox <- 0
  temp_eff <- 0
  temp_n <- 0
  for (j in 1:n1) {
    temp_tox <- sim_dose2(temp_tox,prior_tox,1)
    temp_eff <- sim_dose2(temp_eff,prior_eff,1)
    temp_n <- temp_n + 1
  }
  if(temp_eff < r1) { #stop trial 1: futile
    result[i] <- 1
    toxs[i] <- temp_tox
    effs[i] <- temp_eff
    ns[i] <- temp_n
  } 
  else { #run rest of trial
    for (k in n1:n) {
      temp_tox <- sim_dose2(temp_tox,prior_tox,1)
      temp_eff <- sim_dose2(temp_eff,prior_eff,1)
      temp_n = temp_n + 1
    }
    # 2: full, but not acceptable
    if (temp_eff < r) {
      result[i] <- 2
    } 
    # 3: full, and acceptable
    else {
      result[i] <- 3
    }
    toxs[i] <- temp_tox
    effs[i] <- temp_eff
    ns[i] <- temp_n
  }
}
  temp <- rep(0,3)
  for (i in 1:3) {
    temp[i] <- sum(result == i)/nsims
  }
  info_mat <- matrix(c(temp, mean(toxs), mean(effs), ns[i]),ncol = 6)
  colnames(info_mat) <- c("Early Futility Stop", "Complete, Not Acceptable", 
                          "Complete, Acceptable", "Avg Tox", "Avg Eff", "Avg N")
  rownames(info_mat) <- c("")
  return(info_mat)
}

sim_dose2 <- function(nevent,est_curv,people) {
  cohort <- rbinom(people,1,est_curv)
  nevent <- nevent + sum(cohort)
  return(nevent)
}

s2_mono <- stage2sim(nsims2,r1,r,n1,n,pi_tox,pi_eff_mono)

s2_mono

s2_lvl <- stage2sim(nsims2,r1,r,n1,n,pi_tox,pi_eff_lvl)

s2_lvl

s2_quad <- stage2sim(nsims2,r1,r,n1,n,pi_tox,pi_eff_quad)

s2_quad



