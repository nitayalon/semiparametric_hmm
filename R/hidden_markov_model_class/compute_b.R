#' Computing the normalization size b for tuple (theta,p)

# Set
setGeneric(name = "computeB",
           def = function(theObject)
           {
             standardGeneric("computeB")
           }
)

setMethod(f = "computeB", 
          signature = "HiddenMarkovModel", 
          definition = function(theObject)
          {
            conf_file <- config::get(file = "~/MastersThesis/Main/Code/Mock_data_classes/hidden_markov_model_class/hmm_config.yml")
            observation <- getObservationsFromHMM(theObject)
            n <- length(observation)
            theta <- getTheta(theObject)
            p <- theObject@p
            b <- NewtonRaphson(observation, theta, p)
            theObject@b <- b
            return(theObject)
          }
)
