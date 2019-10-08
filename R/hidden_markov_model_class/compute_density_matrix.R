#' @title Compute semi parametric density
#' @description Computes density mtrix for all the observation, stors the results in a matrix
#' @param hmm - HiddenMarkovModel object
#' @param b - Normalizing size
#' @return Matrix with observation density based on the semi parametric model

# Set

setGeneric(name = "computeFullDensityMatrix",
           def = function(theObject)
           {
             standardGeneric("computeFullDensityMatrix")
           }
)

setMethod(f = "computeFullDensityMatrix", 
          signature = "HiddenMarkovModel", 
          definition = function(theObject)
          {
            observation <- getObservationsFromHMM(theObject, "Vector")
            n_obs = length(observation)
            p <- theObject@p
            theta <- getTheta(theObject)
            b0 <- theObject@b
            theObject@semi_parametric_distribution <- sapply(observation, 
                  function(x){computeDensityMatrix(n_obs, x, b0, theta, p, send_diag = T)})
            return(theObject)
          }
)

# Get

setGeneric(name = "getDensityMatrix",
           def = function(theObject,b)
           {
             standardGeneric("getDensityMatrix")
           }
)

setMethod(f = "getDensityMatrix", 
          signature = "HiddenMarkovModel", 
          definition = function(theObject,b)
          {
            return(theObject@semi_parametric_distribution)
          }
)
