#' @title compute gradient embedded HMM
#' @description compute gradient for embedded HMM
#' @return gradient for MLE

# Set

setGeneric(name = "computeGradientEmbeddedHMM",
           def = function(theObject)
           {
             standardGeneric("computeGradientEmbeddedHMM")
           }
)

setMethod(f = "computeGradientEmbeddedHMM", 
          signature = "embeddedHMMModel", 
          definition = function(theObject)
          {
            gradient <- computeGradient(getHMM(theObject))
            theObject@Gradient <- gradient
            return(theObject)          
          }
)

# Get
setGeneric(name = "getGradient",
           def = function(theObject)
           {
             standardGeneric("getGradient")
           }
)

setMethod(f = "getGradient", 
          signature = "embeddedHMMModel", 
          definition = function(theObject)
          {
            return(theObject@Gradient)
          }
)

