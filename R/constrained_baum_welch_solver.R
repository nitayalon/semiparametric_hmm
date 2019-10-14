#' @title Solve constrained Baum-Welch equations
#' @description Solves the constrained BW equations
#' @param bw_transition_matrix - Matrix, the result of BW
#' @param observations - Vector, vector of observations
#' @param theta - Numeric, exponential load parameter
#' @param p - Numeric, the stationary distribution parameter
#' @param b - Numeric, partition function size
#' @return MLE for transition matrix


constranedBaumWelchSolver <- function(bw_transition_matrix,observations,
                                      theta,p,b,pseudo_count)
  {
  polynom_coefficients <- computePolynomCoefficients(bw_transition_matrix,p)

  s <-
    tryCatch(
      {solveForS(polynom_coefficients)},
      error=function(e){
        return(0)
      }
    )

  t <- solveForT(s,polynom_coefficients)

  if(length(s) > 1)
  {
    s_min <- s[1]
    t_min <- t[1]
    q_min <- p/(1-p) * s_min
    if((q_min) > 1){q_min = 1}
    first_row <- c(1-s_min-t_min, s_min, t_min)
    second_row <- c(q_min / 2, 1 - q_min, q_min / 2)
    third_row <- c(t_min, s_min, 1-s_min-t_min)
    transition_matrix_min <- matrix(c(first_row, second_row, third_row),
                                    ncol = 3, nrow = 3, byrow = T)
    f_min = rcpp_forwardVector(observations,transition_matrix_min,
                               theta,p,b)
    non_normalized_llk_min = f_min[,length(observation)]
    average_min_llk <- mean(non_normalized_llk_min)
    delta_min_llk <- non_normalized_llk_min - average_min_llk
    llk_min <- average_min_llk + log(sum(exp(delta_min_llk)))
    s_max <- s[2]
    t_max <- t[2]
    q_max <- p/(1-p) * s_max
    if((q_max) > 1){q_max = 1}
    first_row <- c(1-s_max-t_max, s_max, t_max)
    second_row <- c(q_max / 2, 1 - q_max, q_max / 2)
    third_row <- c(t_max, s_max, 1-s_max-t_max)
    transition_matrix_max <- matrix(c(first_row, second_row, third_row),
                                    ncol = 3, nrow = 3, byrow = T)
    f_max = rcpp_forwardVector(observations,transition_matrix_max,
                               theta,p,b)
    non_normalized_llk_max = f_max[,length(observation)]
    average_max_llk <- mean(non_normalized_llk_max)
    delta_max_llk <- non_normalized_llk_max - average_max_llk
    llk_max <- average_max_llk + log(sum(exp(delta_max_llk)))
    if(llk_max > llk_min)
    {
      transition_matrix = ((transition_matrix_max + pseudo_count) / apply(transition_matrix_max + pseudo_count,1,sum))
      return(list(bw_transition_matrix=transition_matrix,
                  s_value = s_max,
                  t_value = t_max))
    }
    else
    {
      transition_matrix = ((transition_matrix_min + pseudo_count) / apply(transition_matrix_min + pseudo_count,1,sum))
      return(list(bw_transition_matrix=transition_matrix,
                  s_value = s_max,
                  t_value = t_max))
    }
  }
  q <- p/(1-p) * s
  if((q) > 1){q = 1}
  first_row <- c(1-s-t, s, t)
  second_row <- c(q / 2, 1 - q, q / 2)
  third_row <- c(t, s, 1-s-t)
  TransitionMatrix <- matrix(c(first_row, second_row, third_row),
                             ncol = 3, nrow = 3, byrow = T)
  transition_matrix = ((TransitionMatrix + pseudo_count) / apply(TransitionMatrix + pseudo_count,1,sum))
  return(
    list(bw_transition_matrix=transition_matrix,
         s_value = s,
         t_value = t)
  )
}
