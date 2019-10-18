#' Create transition matrix for Semi-Parametric model
#' @param s - Numeric
#' @param t - Numeric
#' @param p - Numeric, default null, the value of the stationary vector
#' @param check_input - Logical, indicate if to check the caonstraints
#' @param stationary_matrix - Logical, indicate if the transition matrix is assumed to by stationary
#' @return transition matrix

createTransitionMatrix <- function(s, t, p = NULL, check_input = T,
                                   stationary_matrix = T){

  #Input validation:
  if(any(!is.numeric(c(s,t)))){stop("The input must be numeric")}

  if(check_input){
    if(any(c(s,t) > 1)){stop("The input must be less than 1")}
    if(s + t > 1){
      warning("s and t must sum to less than 1")
      s <- s / (s + t)
      t <- t / (s + t)
    }
  }

  if(stationary_matrix)
  {
    stopifnot(p > 0 & p < 1)
    s <- (1-p)/p * 0.01
    q <- p/(1-p) * s
    t = 0.02 - s
  }

  first_row <- c(max(1 - s - t,0), s ,t)
  second_row <- c(q / 2 , 1 - q ,q / 2)
  third_row <- rev(first_row)
  transition_matrix <- matrix(c(first_row, second_row, third_row), nrow = 3, ncol = 3, byrow = T)

  return(transition_matrix)
}

#' An S4 class to represent a HMM with exponential perturbations.
#'
#' @slot sigma numeric, the data standard deviation
#' @slot epsilon numeric, the proportion of the mean wrt to std
#' @slot theta numeric, exponential perturbation size
#' @slot p numeric, stationary vector proportion
#' @slot s numeric, transition probability to non-perturbed state
#' @slot t numeric, transition probability between perturbed states
#' @slot n_obs numeric, length of sample
#' @slot type character, type of conditional distribution
#' @slot transition_matrix matrix, HMM transition matrix
#' @slot observations list, Observed data and Hidden data


hiddenMarkovModel <- setClass(

  # Set class name
  "HiddenMarkovModel",

  # Set class properties
  slots = c(
    sigma = "numeric",
    epsilon = "numeric",
    theta = "numeric",
    p = "numeric",
    s = "numeric",
    t = "numeric",
    n_obs = "numeric",
    type = "character",
    transition_matrix = "matrix",
    observations = "list"
  ),

  # Constructor
  prototype=list(
    sigma = 1,
    epsilon = 0.1,
    theta = NULL,
    p = 2/3,
    s = 1/400,
    t = 7/400,
    n_obs = 10000,
    type = "Normal",
    transition_matrix = createTransitionMatrix(1/400,7/400,2/3,stationary_matrix = T),
    observations = NULL
  ),

  validity=function(object)
  {
    errors <- character()

    if((object@p < 0) || (object@p > 1))
    {
      msg <- "A mixture proportion must by between 0 and 1"
      errors <- c(errors,msg)
    }
    if(object@sigma < 0)
    {
      msg <- "Scale parameter must be positive"
      errors <- c(errors,msg)
    }
    if(!(object@type %in% c("Normal","Laplace")))
    {
      msg <- "Only ADL and normal distributions are supported"
      errors <- c(errors,msg)
    }
    if(object@n_obs <= 0)
    {
      msg <- "Number of observations must be positive and natural"
      errors <- c(errors,msg)
    }

    if(any(c(object@s,object@t) < 0) || any(c(object@s,object@t) > 1))
    {
      msg <- "All the transition matrix entries must be between 0 and 1"
      errors <- c(errors,msg)
    }

    if(object@s + object@t > 1)
    {
      msg <- "t and s must be less than 1!"
      errors <- c(errors,msg)
    }

    if(length(errors) == 0) TRUE else errors
  }

)

#' @title Create transition matrix function
#' @description Create a 3*3 transition matrix for the semi parametric mode
#' @param s - float in range(0,1) indicating the probability of transition from 1/3->2
#' @param t - float in range(0,1) indicating the probability of transition from 1/3->1/3
#' @return transition matrix

### Set

setGeneric(name = "createTransitionMatrixForHMM",
           def = function(theObject)
           {
             standardGeneric("createTransitionMatrixForHMM")
           }
)

setMethod(f = "createTransitionMatrixForHMM",
          signature = "HiddenMarkovModel",
          definition = function(theObject)
          {
            if(any(c(theObject@t,theObject@s,theObject@p) < 0))
            {
              cat("Transition probalities cannot be negative","\n")
              return(theObject)
            }
            if(any(c(theObject@t,theObject@s,theObject@p) > 1))
            {
              cat("Transition probalities cannot exceed 1","\n")
              return(theObject)
            }
            tryCatch(
              {
                theObject@transition_matrix <- createTransitionMatrix(theObject@s,theObject@t,theObject@p,stationary_matrix = T)
              },
              error=function(cond)
              {
                message("There's a problem with creating transition matrix")
                message(cond)
                return(NULL)
              }
            )
            return(theObject)
          }
)


setGeneric(name = "setTransitionMatrixForHMM",
           def = function(theObject,transition_matrix)
           {
             standardGeneric("setTransitionMatrixForHMM")
           }
)

setMethod(f = "setTransitionMatrixForHMM",
          signature = "HiddenMarkovModel",
          definition = function(theObject,transition_matrix)
          {
            if(any(transition_matrix < 0))
            {
              cat("Transition probalities cannot be negative","\n")
              return(theObject)
            }
            if(any(transition_matrix > 1))
            {
              cat("Transition probalities cannot exceed 1","\n")
              return(theObject)
            }
            tryCatch(
              {
                theObject@transition_matrix <- transition_matrix
              },
              error=function(cond)
              {
                message("There's a problem with creating transition matrix")
                message(cond)
                return(NULL)
              }
            )
            return(theObject)
          }
)
### Get

setGeneric(name = "getTransitionMatrix",
           def = function(theObject)
           {
             standardGeneric("getTransitionMatrix")
           }
)

setMethod(f = "getTransitionMatrix",
          signature = "HiddenMarkovModel",
          definition = function(theObject)
          {
            return(theObject@transition_matrix)
          }
)

setGeneric(name="plotHmmData",
           def=function(theObject)
           {
             standardGeneric("plotHmmData")
           }
)

setMethod(f="plotHmmData",
          signature="HiddenMarkovModel",
          definition=function(theObject)
          {
            return(NULL)
          }
)



