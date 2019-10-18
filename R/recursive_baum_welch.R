#' @title Apply Baum-Welch algorithm recursively
#' @description Apply BW single iteration
#' @param current_transition_matrix - Matrix, the current estiamtion of transition matrix
#' @param observations - Vector, vector of observations
#' @param theta - Numeric, exponential load parameter
#' @param p - Numeric, the stationary distribution parameter
#' @param b - Numeric, partition function size
#' @return MLE for transition matrix

baumWelchRecursion = function(current_transition_matrix,
                              observations,
                              theta,
                              p,
                              b,
                              constrained = FALSE,
                              config_file = NULL)
{
  bw_transition_matrix <- rcpp_RecursiveBaumWelch(observations,
                                                  current_transition_matrix,
                                                  theta,p,b)

  constrained_bw_results <- constrainedBaumWelchSolver(bw_transition_matrix,observations,
                                                      theta,p,b,config_file$pseudoCount)

  updated_forward <- rcpp_forwardVector(observations,
                                        constrained_bw_results$bw_transition_matrix,
                                        theta, p, b)
  non_normalized_llk <- updated_forward[length(observations),]
  average_llk <- mean(non_normalized_llk)
  delta_llk <- non_normalized_llk - average_llk
  llk <- average_llk + log(sum(exp(delta_llk)))
  return(
    list(updated_transition_matrix = constrained_bw_results$bw_transition_matrix,
         llk = llk,
         s_value = constrained_bw_results$s,
         t_value = constrained_bw_results$t,
         forward = updated_forward)
  )
}


