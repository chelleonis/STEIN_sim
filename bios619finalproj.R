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
ndose <- length(dose_levels)

  
#variable setup
phi_L = log( (1-phi_1) / (1-phi_0) ) / log( (phi_0*(1-phi_1))/ (phi_1*(1-phi_0)) ) 
phi_U = log( (1-phi_0) / (1-phi_2) ) / log( (phi_2*(1-phi_0))/ (phi_0*(1-phi_2)) ) 
psi_opt = log( (1-psi_1) / (1-psi_2) ) / log( (psi_2*(1-psi_1))/ (psi_1*(1-psi_2)) )
#loop variable setup
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
  
  
}
#step 3: repeat step 2 until max sample size

#utility function at the end to calculate optimal dose
#in remark 3, U(pj,qj) = qj - w1pj - w2 pjI(pj > phi_0)


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

