#' Generate samples of hidden states from HMM
#' @param n_obs - Numeric
#' @param transition_matrix - Matrix, transition matrix
#' @param delta - Numeric, stationary vector
#' @return vector of length \code{n_obs} of hidden states

generateHiddenSequance <- function(n_obs, transition_matrix, delta){
  if(dim(transition_matrix)[1] != dim(transition_matrix)[2]){
    stop("The transition matrix is not square")
  }
  if(sum(delta) != 1){
    stop("Stationary vector must sum to 1")
  }
  if(!all ( apply( transition_matrix,1,sum) == 1 ) ){
    stop("The transition matrix rows don't sum to 1")
  }
  m <- dim(transition_matrix)[1]
  states <- numeric(n_obs)
  states[1] <- sample(m, 1, prob = delta)
  for (i in 2:n_obs){
    states[i] <- sample(m, 1 , prob = transition_matrix[states[i-1], ])
  }
  return(states)
}

#' Generate observation from HMM
#' @param hidden_states - Numeric, vector of hidden states
#' @param parameters_for_distribution - List, emission distribution parameters
#' @param type - Character, type of emission distribution
#' @return list of hidden states and observations

generateObservedSequance <- function(hidden_states,
                                     parameters_for_distribution,
                                     type = c("Normal", "Laplace"))
{
  if( any(parameters_for_distribution[,2] < 0) ){
    stop("Scale params can't be negative")
  }
  cond_distribution <- match.arg(type)
  observations <- switch(cond_distribution,
                         Normal = sapply(hidden_states, function(x) sampleDataFromShiftedNormalDistribution(1,
                                                                                                            parameters_for_distribution[x,1],
                                                                                                            0,
                                                                                                            parameters_for_distribution[x,2])),
                         Laplace = sapply(hidden_states, function(x) sampleDataFromDoubleExponentialDistribution(1,
                                                                                                                 parameters_for_distribution[x,1]))
  )
  return(list(
    observations = observations,
    hidden_states = hidden_states
  )
  )
}

#' Sample from Normal distribution
#' @param n - Numeric, number of samples
#' @param theta - Numeric, perturbation parameter
#' @param mu - Numeric, location parameter
#' @param simga2 - Numeric, variance parameter
#' @return list of hidden states and observations

sampleDataFromShiftedNormalDistribution <- function(n, theta, mu = 0, sigma2 = 1){
  if(sigma2 <= 0){
    warning("Negative variance")
    sigma2 <- 1
    }
  samples <- rnorm(n,mu + theta, sqrt(sigma2))
  return(samples)
}

#' Sampling from A-Symmetric double Laplace
#' @param n - int, number of observations
#' @param mu - numeric, drift parameter
#' @return x - n samples from A-symmetric laplace distribution

sampleDataFromDoubleExponentialDistribution <- function(n, mu = 0){
  stopifnot(mu < 1)
  ind <- runif(n)
  direction <- ifelse(mu < 0 , "-1", "1")
  theta <- abs(mu)
  y <- rexp(n,1)
  h_x <- 1 / (1 - theta) * y
  g_x <- 1 / (1 + theta) * y
  p_positive <- (1-theta) / 2
  positive_ind <- ind <= p_positive
  negative_sign <- switch (direction,
                           "1" = c(-1,1),
                           "-1" = c(1,-1))
  x <- negative_sign[1] * g_x * positive_ind + negative_sign[2] * h_x * (1 - positive_ind)
  return(x)
}

createObservations <- function(n_obs,
                               probs_matrix_params,
                               p,
                               is_stationarity = T,
                               init_vector,
                               parameters_for_distribution,
                               conditional_distribution_type){
  if(any(probs_matrix_params > 1)){stop("The transition matrix params are not valid")}
  if(n_obs <= 1){n_obs = 1000}
  return(generateObservedSequance(
    generateHiddenSequance(n_obs,
                           createTransitionMatrix(probs_matrix_params[1], probs_matrix_params[2], p, stationary_matrix = T), init_vector),
    parameters_for_distribution,conditional_distribution_type))
}

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

setGeneric(name = "generateSampleFromHMM",
           def = function(theObject)
           {
             standardGeneric("generateSampleFromHMM")
           }
)

setMethod(f = "generateSampleFromHMM",
          signature = "HiddenMarkovModel",
          definition = function(theObject)
          {
            p <- theObject@p
            init_probs <- c(p/2, 1-p, p/2)
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
            return(list(hidden_markov_model_with_observations = theObject,
                        observations = obs))
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
          definition = function(theObject,type = c("Vector","Data.Frame"))
          {
            type_of_data = match.arg(type)
            obs <- switch (type_of_data,
                           Data.frame = data.frame(observation = theObject@observations$observations),
                           Vector = theObject@observations$observations
            )
            return(obs)
          })
