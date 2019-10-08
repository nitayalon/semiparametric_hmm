# Generate sequence from HMM object

sampleFromHMM <- function(n_obs,p,s,t,theta,sigma,type, stationarity = T)
{
  init_probs <- c(p/2,1-p,p/2)
  probs_matrix_params = c(s,t)
  parameters_for_conditional_distribution = 
    data.frame(mu = c(-theta,0,theta), sigma = rep(sigma,3))
  
  obs <- createObservations(n_obs, 
                            init_probs,
                            p,
                            stationarity,
                            probs_matrix_params = probs_matrix_params, 
                            parameters_for_distribution = 
                              parameters_for_conditional_distribution,
                            conditional_distribution_type = type)
  return(obs)
}

# Set

setGeneric(name = "sampleFromHMM",
           def = function(theObject)
           {
             standardGeneric("sampleFromHMM")
           }
)

setMethod(f = "sampleFromHMM", 
          signature = "HiddenMarkovModel", 
          definition = function(theObject)
          {
            init_probs <- computeStartingVector(theObject)
            probs_matrix_params = c(theObject@s,theObject@t)
            parameters_for_conditional_distribution = 
              data.frame(mu = c(-theObject@epsilon * theObject@sigma
                                ,0,theObject@epsilon * theObject@sigma), sigma = rep(theObject@sigma,3))
            obs <- createObservations(theObject@n_obs, 
                                      probs_matrix_params = probs_matrix_params, 
                                      theObject@p,
                                      T,
                                      init_probs,
                                      parameters_for_distribution = 
                                        parameters_for_conditional_distribution,
                                      conditional_distribution_type = theObject@type)
            theObject@observations <- obs
            return(theObject)
          })

# Get 
setGeneric(name = "getObservationsFromHMM",
           def = function(theObject,type = "Vector")
           {
             standardGeneric("getObservationsFromHMM")
           }
)

setMethod(f = "getObservationsFromHMM", 
          signature = "HiddenMarkovModel", 
          definition = function(theObject,type = "Vector")
          {
            obs <- switch (type,
                           "Data.Frame" = data.frame(observation = theObject@observations$observations),
                           "Vector" = theObject@observations$observations
            )
            return(obs)
          })
