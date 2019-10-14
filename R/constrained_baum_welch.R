#' @title Apply Baum-Welch algorithm
#' @description Computes the MLE for the transition matrix using constrained BW algorithm
#' @param hmm - HiddenMarkovModel object
#' @param constrained - Boolean, solve for general case or constrained?
#' @param number_of_iterations - Numeric, number of iteration
#' @param config_file - List, configuration file
#' @return MLE for transition matrix

setGeneric(name = "constrainedBaumWelch",
           def = function(hidden_markov_model_object,
                          constrained = FALSE,
                          number_of_iterations = NULL,
                          config_file = NULL)
           {
             standardGeneric("constrainedBaumWelch")
           }
)

setMethod(f = "constrainedBaumWelch",
          signature = "HiddenMarkovModel",
          definition = function(hidden_markov_model_object,
                                constrained = FALSE,
                                number_of_iterations = NULL,
                                config_file= NULL)
          {
            if(is.null(number_of_iterations))
            {
              if(is.null(config_file$em_iteration))
              {
                number_of_iterations <- 5000
              }else{
                number_of_iterations <- config_file$em_iteration
              }
            }
            observations <- getObservationsFromHMM(hidden_markov_model_object)
            if(is.null(observations)){observations <- generateSampleFromHMM(hidden_markov_model_object)}
            theta <- getTheta(hidden_markov_model_object)
            p <- hidden_markov_model_object@p
            b <- NewtonRaphson(observations,theta,p,config_file$NRiteration)
            starting_transition_matrix <- createStartingTransitionMatrixForEM()
            llk <- c()
            s_ratio <- c()
            current_s <- c()
            current_t <- c()
            for(i in 1:number_of_iterations)
            {
              # Expectation Step (Calculate expected transition_matrixransitions and Emissions)
              bw = baumWelchRecursion(starting_transition_matrix,
                                      observations,
                                      theta,
                                      p,
                                      b,
                                      constrained = constrained,
                                      config_file = config_file)
              llk <- c(llk, bw$llk)
              transition_matrix  = bw$updated_transition_matrix
              # Maximization Step (Maximise Log-Likelihood for transition_matrixransitions and Emissions-Probabilities)
              current_s <- c(current_s,bw$s_value)
              current_t <- c(current_t,bw$t_value)
              if(i > 3)
              {
                s_ratio <- c(s_ratio, (current_s[i] - current_s[i - 1]) /
                               (current_s[i - 1] - current_s[i - 2]))
              }

              if(i > 1 && round(current_s[i],5) ==
                 round(current_s[i-1],5))
              {
                break(cat(sprintf("Number of effective iterations %d",i), "\n"))
              }
            }
            return(list(mle_transition_matrix = transition_matrix,
                        llk = llk,
                        s_ratio = s_ratio,
                        current_s = current_s,
                        current_t = current_t,
                        forward = bw$forward))
          })

createStartingTransitionMatrixForEM <- function(n_states = 3)
{
  first_row <- c(1-(7/400 + 1/400), 1/400, 7/400)
  third_row <- rev(first_row)
  middle_row <- c(1/2 * (2/3 / 1/3) * 1/400, 1 - (2/3 / 1/3) * 1/400, 1/2 * (2/3 / 1/3) * 1/400)
  return(matrix(c(first_row,middle_row,third_row), nrow = n_states, byrow = T))
}


