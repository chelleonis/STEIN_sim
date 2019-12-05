#bios 619 final project multiple myeloma smulation
#library(bcrm)
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
doses = c(5,10,15,19)/100
pi_nmax = 35 #max enrollment for phase I (max 4-6 subjects/month)
dlt = 0.2 #target
prior_tox = c(0,0,0,0) 

#Phase I
sim1 = crmsim( PI = c(5,10,15,19)/100, #true probabilities
  prior = c(5,10,15,19)/100, #initial guesses 
  target = 1/3, #target DLT rate, standard is 1/3
  n = 35, #trial sample size
  x0 = 1, #initial dose
  nsim = 100,
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

#now do the seamless oWo
#STEIN design, reminescent of MPTI
#adapted from Lin, Yin 2017

# pj tox probability @ level j -> estimated via xj/nj
# qj efficacy probability @ level j -> estimated via yj/nj
# tox prob monotonic assumption, no constraint
# nj #of patients treated at level j
# xj # of tox outcomes
# yj # of efficacy outcomes

#GOAL: find Optimal Biological Dose
#integrate toxicity and efficacy finding in one go

# 0 < phi1 < phi0 < phi2 < 1 
#Design parameters:
#phi0 = 0.3 (default) -> target tox probability / acceptable tox probability
#phi1: highest tox probability that is overly safe (if dose j pj < phi1 next level still ok)
#phi2: lowest tox prob that is excessiely toxic (if dose j has pj > phi2 then de-escalate)

#psi1,psi2: efficacy monitoring parameters
#
#impelement boundary values xd
#default reccomendation, phi1 = 0.75*phi0 , phi2 = 1.25*phi0
#psi1 & psi2 from clinician, (drug higher response, then increase)
# psi1 = 0.3 , psi2 = 0.8 default values

#early stopping Pr(pj > phi0 | xj) > 0.95) [unif(0,1) prior for each pj)]
#futility: Pr (qj < psi1 | yj ) > 0.98, dose level j eliminated for futility
#if all suck (futile) then terminate early

#ALGORITHM CODE
STEIN_sim <- function(nsims = 100, npatients, dose_levels, 
                phi_0 = 0.3, phi_1 = 0.75*phi_0, phi_2 = 1.25*phi_0,
                psi_1 = 0.3, psi_2 = 0.8,
                ncohort = 1) {
#loop over simulations

#variable setup
phi_L = log( (1-phi_1) / (1-phi_0) ) / log( (phi_0*(1-phi_1))/ (phi_1*(1-phi_0)) ) 
phi_U = log( (1-phi_0) / (1-phi_2) ) / log( (phi_2*(1-phi_0))/ (phi_0*(1-phi_2)) ) 
psi_opt = log( (1-psi_1) / (1-psi_2) ) / log( (psi_2*(1-psi_1))/ (psi_1*(1-psi_2)) )
#loop variable setup
ndose <- length(dose_levels)
ntox <- rep(0,ndose)
neff <- rep(0,ndose)
n_lvl <- rep(0,ndose)
#step 1: treat first cohort at lowest dose level
j <- 1 #current dose level
#step 2: follow rules: 
#   a.) p_estj >= phi_U, de-escalate  "inadmissible"
#   b.) p_estj < phi_u and q_estj >= psi, retain dose "promising"
#   c.) o.w select next dose based on: "exploratory"
#     j_next = argmax( Pr(qj > psi | yj')) , largest posterior probability in doses
# in admissible set, Unif(0,1) prior to all qj's (ties go to lower dose)
for (i in 1:npatients%/%ncohort+1) {
  #simulate cohort on dose level
  np <- min(npatients-ncohort*i,ncohort)
  if (np > 0) {
  ntox[j] <- sim_dose(j,ntox,tox_levels,np)
  neff[j] <- sim_dose(j,neff,eff_levels,np)
  n_lvl[j] <- n_lvl[j] + np
  }
  #calculate p_estj and q_estj
  p_estj <- ntox[j]/n_lvl[j]
  q_estj <- neff[j]/n_lvl[j]
  
  if (p_estj >= phi_U) { #a
    j <- max(1,j-1)
  }
  else if(p_estj < phi_U & q_estj >= psi_opt) { #b
    j <- j #not even sure this line is necessary
  }
  else { #c
    #admissible set:
    if(p_estj <= phi_L) {
      set <- c(j-1,j,j+1) #make sure this is not greater or less than set 
    }
    else {
      set <- c(j-1,j) #same as above
    }
    #dunno how to do Pr(q_estj' > psi_opt | yj') yet
    
  }
  #early stopping Pr(pj > phi0 | xj) > 0.95) [unif(0,1) prior for each pj)]
  #futility: Pr (qj < psi1 | yj ) > 0.98, dose level j eliminated for futility
  #if all suck (futile) then terminate early
  
  
  #utility function at the end to calculate optimal dose
  #in remark 3, U(pj,qj) = qj - w1pj - w2 pjI(pj > phi_0)
}
#step 3: repeat step 2 until max sample size

#returns a matrix of selected doses, ntox, and idk what else
}

#helper function to sim at the dose level maybe....
sim_dose <- function(j,nevent,est_curv,people) {
  #MOVE TO FUNCTION, eff probs have separate fxn
  cohort <- rbinom(np,1,tox_levels[j])
  nevent[j] <- nevent[j] + sum(cohort)
  return(nevent)
}



#actual simulation
#35 people recruited, cohort size 3
#dose levels: above
#target tox 1/3

