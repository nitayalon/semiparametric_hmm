#' @title Apply Baum-Welch algorithm
#' @description Computes the MLE for the transition matrix using BW algorithm and Meilijson EM boosting
#' @param hmm - HiddenMarkovModel object
#' @return MLE for transition matrix

setGeneric(name = "estimateTransitionMatrix",
           def = function(theObject, transition_matrix = NULL, 
                          constraint = FALSE,
                          density_matrix = NULL, 
                          number_of_iterations = NULL,
                          config_file)
           {
             standardGeneric("estimateTransitionMatrix")
           }
)

setMethod(f = "estimateTransitionMatrix", 
          signature = "HiddenMarkovModel", 
          definition = function(theObject, transition_matrix = NULL, constraint = FALSE,
                                density_matrix = NULL, number_of_iterations = NULL,
                                config_file)
          {
            # browser()
            if(is.null(number_of_iterations))
            {
              number_of_iterations <- config_file$em_iteration
            }
            theObject <- computeB(theObject)
            llk <- c()
            s_ratio <- c()
            current_s <- c()
            current_t <- c()
            temphmm_object = theObject
            if(!is.null(transition_matrix))
            {
              temphmm_object = setTransitionMatrixForHMM(temphmm_object, transition_matrix)
            }
            
            for(i in 1:number_of_iterations)
            {
              
              # Expectation Step (Calculate expected transition_matrixransitions and Emissions)
              bw = baumWelchRecursion(temphmm_object, constraint = constraint,
                                      config_file = config_file)
              llk <- c(llk, bw$llk)
              transition_matrix  = bw$TransitionMatrix
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
              temphmm_object = bw$updated_hmm
              if(is.null(temphmm_object)){browser()}
              temphmm_object@s = transition_matrix[1,2]
              temphmm_object@t = transition_matrix[1,3]
            }
            return(list(hmm_object = temphmm_object,
                        llk = llk,
                        s_ratio = s_ratio,
                        current_s = current_s,
                        current_t = current_t,
                        forward = bw$forward))
          })




