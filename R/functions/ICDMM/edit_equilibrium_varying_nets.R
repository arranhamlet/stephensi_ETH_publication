edit_equilibrium_varying_nets <- function(wh, population_split){
  # Need to add population split and use to edit initial conditions
  wh$state$pop_split <- population_split # Population split
  cov <- wh$mpl$cov
  
  # This code here only works when IRS is turned off
  # Divide the no intervention compartment by cov[1] to get original eqm soln
  tmp_init_S <- wh$state$init_S[,,1] / cov[1]
  wh$state$init_S[,,1] <- tmp_init_S * wh$state$pop_split[1]
  wh$state$init_S[,,2] <- tmp_init_S * wh$state$pop_split[2]
  tmp_init_T <- wh$state$init_T[,,1] / cov[1]
  wh$state$init_T[,,1] <- tmp_init_T * wh$state$pop_split[1]
  wh$state$init_T[,,2] <- tmp_init_T * wh$state$pop_split[2]
  tmp_init_D <- wh$state$init_D[,,1] / cov[1]
  wh$state$init_D[,,1] <- tmp_init_D * wh$state$pop_split[1]
  wh$state$init_D[,,2] <- tmp_init_D * wh$state$pop_split[2]
  tmp_init_A <- wh$state$init_A[,,1] / cov[1]
  wh$state$init_A[,,1] <- tmp_init_A * wh$state$pop_split[1]
  wh$state$init_A[,,2] <- tmp_init_A * wh$state$pop_split[2]
  tmp_init_U <- wh$state$init_U[,,1] / cov[1]
  wh$state$init_U[,,1] <- tmp_init_U * wh$state$pop_split[1]
  wh$state$init_U[,,2] <- tmp_init_U * wh$state$pop_split[2]
  tmp_init_P <- wh$state$init_P[,,1] / cov[1]
  wh$state$init_P[,,1] <- tmp_init_P * wh$state$pop_split[1]
  wh$state$init_P[,,2] <- tmp_init_P * wh$state$pop_split[2]
  
  return(wh)
}