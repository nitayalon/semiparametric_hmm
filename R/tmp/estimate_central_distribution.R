#' Estimating the central distribtuion for each pair of (p,theta)
#' This method plots a histogram for the empirical central distribution 

# Set
setGeneric(name = "estimateCentralDistribution",
           def = function(theObject)
           {
             standardGeneric("estimateCentralDistribution")
           }
)

setMethod(f = "estimateCentralDistribution", 
          signature = "HiddenMarkovModel", 
          definition = function(theObject)
          {
            conf_file <- config::get(file = "~/MasterThesis/Main/Code/Mock_data_classes/hidden_markov_model_class/hmm_config.yml")
            observation <- getObservationsFromHMM(theObject)
            n <- length(observation)
            theta <- getTheta(theObject)
            p <- theObject@p
            b <- NewtonRaphson(observation, theta, p)
            central_distribution <- 1 / n * 1 / exponentialShift(theta,p,b,observation) 
            return(central_distribution)
          }
)