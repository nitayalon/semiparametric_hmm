#' @title Apply Grid Search algorithm
#' @description Apply GS algorithm for profile MLE
#' @param observations - Vector, vector of observations
#' @param config_file - List, configuration file
#' @param symmetric - Bool, is the transition matrix symmetric
#' @param n_states - Numeric, number of hidden states
#' @param stationary - Bool, is the process stationary
#' @return Grid Search class
gridSearch <- function(observations, config_file, symmetric = T,
                       n_states = 3, stationary = T){
  #1) Validate input
  #2) Create 10*10 grid of theta p
    #2.1) Compute b for each entry
  #3) Apply Constrained BW per entry - return transition matrix, forward, llk
  #4) Find minimum -llk and the corresponding tuple (theta,p,b,Gamma)
  #5) create 10*10 thin grid centered at the above
  #6) Apply 2.1-4 on the thin grid
  #7) Fit parabolid to the mle 3*3 grid, compute hessian
  #8) Compute CI for p,theta
  #9) Return GS object
}
