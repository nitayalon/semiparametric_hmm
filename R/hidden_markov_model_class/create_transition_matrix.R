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
                return("NONONONONO")
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

