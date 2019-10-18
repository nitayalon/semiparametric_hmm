#' @title Apply Baum-Welch algorithm
#' @description Computes the MLE for the transition matrix using constrained BW algorithm
#' @param observations - Vector, vector of observations
#' @param theta - Numeric, exponential load parameter
#' @param p - Numeric, the stationary distribution parameter
#' @param constrained - Boolean, solve for general case or constrained?
#' @param number_of_iterations - Numeric, number of iteration
#' @param config_file - List, configuration file
#' @return MLE for transition matrix

constrainedBaumWelch = function(observations,
                                theta,
                                p,
                                constrained = FALSE,
                                number_of_iterations = NULL,
                                config_file= NULL){
  if(is.null(number_of_iterations))
  {
    if(is.null(config_file$em_iteration))
    {
      number_of_iterations <- 5000
    }else{
      number_of_iterations <- config_file$em_iteration
    }
  }
  b <- NewtonRaphson(observations,theta,p,config_file$NRiteration)
  current_transition_matrix <- createStartingTransitionMatrixForEM()
  llk <- c()
  s_ratio <- c()
  current_s <- c()
  current_t <- c()
  for(i in 1:number_of_iterations)
  {
    # Expectation Step (Calculate expected transition_matrixransitions and Emissions)
    bw = baumWelchRecursion(current_transition_matrix,
                            observations,
                            theta,
                            p,
                            b,
                            constrained = constrained,
                            config_file = config_file)
    llk <- c(llk, bw$llk)
    current_transition_matrix  = bw$updated_transition_matrix
    # Maximization Step (Maximise Log-Likelihood for transition_matrixransitions and Emissions-Probabilities)
    current_s <- c(current_s,bw$s_value)
    current_t <- c(current_t,bw$t_value)
    if(i > 3)
    {
      s_ratio <- c(s_ratio, (current_s[i] - current_s[i - 1]) /
                     (current_s[i - 1] - current_s[i - 2]))
    }

    # if(i > 4 && round(llk[i],5) ==
    #    round(llk[i-1],5))
    # {
    #   cat(sprintf("Number of effective iterations %d",i))
    #   break
    # }
  }
  return(list(mle_transition_matrix = current_transition_matrix,
              llk = llk,
              s_ratio = s_ratio,
              current_s = current_s,
              current_t = current_t,
              forward = bw$forward))
}

createStartingTransitionMatrixForEM <- function(n_states = 3)
{
  first_row <- c(1-(7/400 + 1/400), 1/400, 7/400)
  third_row <- rev(first_row)
  middle_row <- c(1/2 * (2/3 / 1/3) * 1/400, 1 - (2/3 / 1/3) * 1/400, 1/2 * (2/3 / 1/3) * 1/400)
  return(matrix(c(first_row,middle_row,third_row), nrow = n_states, byrow = T))
}


