computeTheta <- function(sigma, epsilon){sigma * epsilon}

# Set

setGeneric("computeTheta")
setGeneric("computeTheta", function(object){
  standardGeneric("computeTheta")
})
setMethod("computeTheta", signature(object = "HiddenMarkovModel"), function(object){
  object@theta <- object@sigma * object@epsilon
  return(object)
})


# Get
setGeneric(name = "getTheta",
           def = function(theObject)
           {
             standardGeneric("getTheta")
           }
)

setMethod(f = "getTheta",
          signature = "HiddenMarkovModel",
          definition = function(theObject)
          {
            return(theObject@sigma * theObject@epsilon)
          }
)
