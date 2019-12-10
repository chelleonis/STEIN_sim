#STEIN design, reminescent of MPTI
#adapted from Lin, Yin 2017
#Author: Allen Li

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
set.seed(619)
#ALGORITHM CODE
STEIN_sim <- function(nsims, npatients, tox_levels, eff_levels, 
                phi_0 = 0.3, phi_1 = 0.75*phi_0, phi_2 = 1.25*phi_0,
                psi_1 = 0.3, psi_2 = 0.8,
                ncohort = 1, w_tox = 1, w_eff = 1 ) {
#loop over simulations
if (npatients%%ncohort == 0) {cohorts = npatients/ncohort} 
else{ cohorts <- npatients%/%ncohort+1 }
ndose <- length(tox_levels)
OBD <- rep(0,ndose)
sim_tox <- matrix(0,nsims,ncol = ndose)
sim_eff <- matrix(0,nsims,ncol = ndose)
sim_level <- matrix(0,nsims,ncol = ndose)
  for(r in 1:nsims) {
#variable setup
phi_L = log( (1-phi_1) / (1-phi_0) ) / log( (phi_0*(1-phi_1))/ (phi_1*(1-phi_0)) ) 
phi_U = log( (1-phi_0) / (1-phi_2) ) / log( (phi_2*(1-phi_0))/ (phi_0*(1-phi_2)) ) 
psi_opt = log( (1-psi_1) / (1-psi_2) ) / log( (psi_2*(1-psi_1))/ (psi_1*(1-psi_2)) )
#loop variable setup
ntox <- rep(0,ndose)
neff <- rep(0,ndose)
n_lvl <- rep(0,ndose)
futile <- rep(0,ndose)
unsafe <- rep(0,ndose)
#step 1: treat first cohort at lowest dose level
j <- 1 #current dose level
#step 2: follow rules: 
#   a.) p_estj >= phi_U, de-escalate  "inadmissible"
#   b.) p_estj < phi_u and q_estj >= psi, retain dose "promising"
#   c.) o.w select next dose based on: "exploratory"
#     j_next = argmax( Pr(qj > psi | yj')) , largest posterior probability in doses
# in admissible set, Unif(0,1) prior to all qj's (ties go to lower dose)

for (i in 1:cohorts) {
  #simulate cohort on dose level
  np <- min(npatients-ncohort*(i-1),ncohort)
  ntox[j] <- sim_dose(j,ntox,tox_levels,np)
  neff[j] <- sim_dose(j,neff,eff_levels,np)
  n_lvl[j] <- n_lvl[j] + np
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
      good_set <- set[set > 0 & set <= ndose]
      j <- c_doserule(good_set,neff[good_set],n_lvl[good_set],psi_opt) 
    }
    else {
      set <- c(j-1,j) #same as above
      good_set <- set[set > 0 & set <= ndose]
      j <- c_doserule(good_set,neff[good_set],n_lvl[good_set],psi_opt) 
    }
    #dunno how to do Pr(q_estj' > psi_opt | yj') yet
  }
  #Safety Rule: Pr(pj > phi0 | xj) > 0.95) [unif(0,1) prior for each pj)]
  if (calc_dl(ntox[j],n_lvl[j],phi_0) > 0.95) {
    unsafe[j] <- 1
  }
  if(all(unsafe == 1) == TRUE) {
    #add params or w/e
    break #stop trial
  }
 
  #futility: Pr (qj < psi1 | yj ) > 0.98, dose level j eliminated for futility
  if (calc_dl(neff[j],n_lvl[j],psi_1) < 1-0.98) {
    futile[j] <- 1
  }
  if (all(futile == 1) == TRUE) {
    print(paste0("Trial stopped early for futility at cohort"), i)
    #add params or w/e
    break #stop trial
  }
  #if all futile then terminate early, TBD
}
#step 3: repeat step 2 until max sample size
#utility function at the end to calculate optimal dose
#in remark 3, U(pj,qj) = qj - w1pj - w2 pjI(pj > phi_0)
OBD[r] <- OBD_select(ndose,n_lvl,ntox,neff, phi_0, w_tox, w_eff)
sim_tox[r,] <- ntox
sim_eff[r,] <- neff
sim_level[r,] <- n_lvl

}
OBD_sims <- rep(0,ndose)
for (i in 1:ndose) {
  OBD_sims[i] <- sum(OBD == i)/nsims
}
#returns a matrix of OBD selection (%), average ntox, and average neffs
info_mat <- matrix(c(OBD_sims,colMeans(sim_tox),colMeans(sim_eff),colMeans(sim_level)),ncol = ndose, byrow = TRUE)
rownames(info_mat) <- c("OBD Selection (%)", "Average Ntox", "Average Neff", "Average Nlvl")
return(info_mat)
}

#helper function to sim at the dose level
sim_dose <- function(j,nevent,est_curv,people) {
  #MOVE TO FUNCTION, eff probs have separate fxn
  cohort <- rbinom(people,1,est_curv[j])
  nevent[j] <- nevent[j] + sum(cohort)
  return(nevent[j])
}

c_doserule <- function(doseset,neff,nlvl,psi) {
  m <- length(doseset) 
  if(m == 1) {
    return(doseset[1])
  }
  else{ #argmax Pr(qj > psi | yj)
    probs <- rep(0,m)
    #print(doseset)
    for (i in 1:m) {
      probs[i] <- calc_dl(neff[i],nlvl[i],psi)
      #print(paste("entry:",i, ":" , probs[i] ))
    }
    max_dose <- which.max(probs)
    return(doseset[max_dose])
  }
}
#unif 0,1 prior default -> beta(alpha = 1, beta = 1)
#doing beta for now
calc_dl <- function(nevents,ntotal,prob) { 
  #posterior is then a+sum xi, beta + n - sum xi
  posterior <- 1-pbeta(prob,1+ nevents, 1 + ntotal - nevents) #checkcheckechk
}

OBD_select <- function(j,nlvl,ntox,neff,phi_0,w1,w2) {
  utility <- rep(0,j)
  tox_est <- est_calc(nlvl,ntox)
  eff_est <- est_calc(nlvl,neff)
  for (i in 1:j) {
    utility[i] <- eff_est[i] - w1*tox_est[i] - w2*tox_est[i]*indicator(tox_est[i],phi_0)
  }
  which.max(utility)
}

est_calc <- function(nlvl,nevent) {
  nlvl[nlvl == 0] <- 1  
  return(nevent/nlvl)
}
indicator <- function(tox,phi) {
  if (tox > phi) {return (1)}
  else { return(0)}
}

#actual simulation
#35 people recruited, cohort size 5
#Target Tox Level: 0.2

prior_tox <- c(5,10,15,19)/100 #given
prior_eff_mono <- c(10,20,30,40)/100
prior_eff_level <- c(10,25,42,50)/100
prior_eff_quad <- c(15,40,35,30)/100

dose_labels <- c(1,2,3,4)

set_mono <- STEIN_sim(nsims = 1000, npatients = 57, 
                      tox_levels = prior_tox, eff_levels = prior_eff_mono, 
                      phi_0 = 0.2,
                      psi_1 = 0.25, psi_2 = 0.45,
                      ncohort = 5, w_tox = 1, w_eff = 1 )
colnames(set_mono) <- dose_labels
set_mono

set_level1 <- STEIN_sim(nsims = 1000, npatients = 57, 
                      tox_levels = prior_tox, eff_levels = prior_eff_level, 
                      phi_0 = 0.2,
                      psi_1 = 0.25, psi_2 = 0.45,
                      ncohort = 5, w_tox = 1, w_eff = 1 )
colnames(set_level1) <- dose_labels
set_level1

set_quad <- STEIN_sim(nsims = 1000, npatients = 57, 
                      tox_levels = prior_tox, eff_levels = prior_eff_quad, 
                      phi_0 = 0.2,
                      psi_1 = 0.25, psi_2 = 0.45,
                      ncohort = 5, w_tox = 1, w_eff = 1 )
colnames(set_quad) <- dose_labels
set_quad
