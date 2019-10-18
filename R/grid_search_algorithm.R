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
  stopifnot(length(observations) > 0)
  if(is.null(config_file))
  {
    stop("Must pass a config file")
  }
  #2) Create 10*10 grid of theta p
  outer_grid <- createOuterGrid()
  #3) Apply Constrained BW per entry - return transition matrix, forward, llk
  outer_grid_search_results <- pbapply(outer_grid, 1, function(x)
  {constrainedBaumWelch(observations, x[1], x[2],
                        constrained = symmetric,
                        number_of_iterations = config_file$em_iteration,
                        config_file = config_file)})
  #4) Find minimum -llk and the corresponding tuple (theta,p,b,Gamma)
  browser()
  outer_mle_tuple <- findMleFromGridSearchResults(outer_grid_search_results,outer_grid)
  #5) create 10*10 thin grid centered at the above
  inner_grid <- createInnerGrid()
  #6) Apply 2.1-4 on the thin grid
  inner_grid_search_results <-pblapply(inner_grid, function(x)
  {constrainedBaumWelch(observations, x$theta, x$p,constrained = symmetric,
                        number_of_iterations = config_file$em_iteration,
                        config_file = config_file)})
  #4) Find minimum -llk and the corresponding tuple (theta,p,b,Gamma)
  inner_mle_tuple <- findMleFromGridSearchResults(inner_grid_search_results,inner_grid)
  #7) Fit parabolid to the mle 3*3 grid, compute hessian
  #8) Return GS object
  gs_results <- new("GridSearchResults",
                    outer_grid = outer_grid,
                    outer_grid_results = outer_grid_search_results,
                    outer_grid_mle = outer_mle_tuple,
                    inner_grid = inner_grid,
                    inner_grid_results = inner_grid_search_results,
                    inner_grid_mle = inner_mle_tuple,
                    paraboloid_mle = NULL,
                    fisher_matrix = NULL)
  return(gs_results)
}

createOuterGrid <- function(){
  theta <- seq(0,1,0.1)
  p <- c(0.01,seq(0.1,0.9,0.1),0.99)
  outer_grid <- expand.grid(theta,p)
  names(outer_grid) <- c("theta","p")
  return(outer_grid)
}

createInnerGrid <- function(theta_center,p_center){
  theta <- seq(max(theta_center - 0.05,0),
               min(theta_center + 0.05,1),
               0.01)
  p <- seq(max(p_center - 0.05,0),
           min(p_center + 0.05,1),
           0.01)
  inner_grid <- expand.grid(theta,p)
  names(inner_grid) <- c("theta","p")
  return(inner_grid)
}

findMleFromGridSearchResults <- function(grid_search_results,grid) {
  llk_list <- sapply(grid_search_results, function(x){x$llk[length(x$llk)]})
  mle_set <- grid[which.min(-llk_list),]
  return(mle_set)
}
