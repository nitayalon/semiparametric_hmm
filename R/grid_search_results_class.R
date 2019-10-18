#' An S4 class for Grdi Search Results.
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


GridSearchResults <- setClass(

  # Set class name
  "GridSearchResults",

  # Set class properties
  slots = c(
    outer_grid = "list",
    outer_grid_results = "list",
    outer_grid_mle = "list",
    inner_grid = "list",
    inner_grid_results = "list",
    inner_grid_mle = "list",
    paraboloid_mle = "list",
    fisher_matrix = "list"
  ),

  # Constructor
  prototype=list(
    outer_grid = NULL,
    outer_grid_results = NULL,
    outer_grid_mle = NULL,
    inner_grid = NULL,
    inner_grid_results = NULL,
    inner_grid_mle = NULL,
    paraboloid_mle = NULL,
    fisher_matrix = NULL
  )
)

setGeneric(name="getOuterGridResults",
           def=function(theObject)
           {
             standardGeneric("getOuterGridResults")
           }
)

setMethod(f="getOuterGridResults",
          signature="GridSearchResults",
          definition=function(theObject)
          {
            return(theObject@outer_grid_results)
          }
)

setGeneric(name="getInnerGridResults",
           def=function(theObject)
           {
             standardGeneric("getInnerGridResults")
           }
)

setMethod(f="getInnerGridResults",
          signature="GridSearchResults",
          definition=function(theObject)
          {
            return(theObject@inner_grid_results)
          }
)

setGeneric(name="plotOuterGridResults",
           def=function(theObject)
           {
             standardGeneric("plotOuterGridResults")
           }
)

setMethod(f="plotOuterGridResults",
          signature="GridSearchResults",
          definition=function(theObject)
          {
            return(NULL)
          }
)

setGeneric(name="plotInnerGridResults",
           def=function(theObject)
           {
             standardGeneric("plotInnerGridResults")
           }
)

setMethod(f="plotInnerGridResults",
          signature="GridSearchResults",
          definition=function(theObject)
          {
            return(NULL)
          }
)

setGeneric(name="computeCI",
           def=function(theObject)
           {
             standardGeneric("computeCI")
           }
)

setMethod(f="computeCI",
          signature="GridSearchResults",
          definition=function(theObject)
          {
            return(NULL)
          }
)
