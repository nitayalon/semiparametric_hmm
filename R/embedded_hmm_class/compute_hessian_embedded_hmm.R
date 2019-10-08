#' @title compute hessian embedded HMM
#' @description compute hessian for embedded HMM
#' @return hessian for MLE

# Set
setGeneric(name = "computeHessianEmbeddedHMM",
           def = function(theObject)
           {
             standardGeneric("computeHessianEmbeddedHMM")
           }
)

setMethod(f = "computeHessianEmbeddedHMM", 
          signature = "embeddedHMMModel", 
          definition = function(theObject)
          {
            mle_hmm <- createParametersGridForHessian(getHMM(theObject))
            Hessian <- computeHessian(mle_hmm)
            theObject@Hessian <- Hessian
            return(theObject)          
          }
)

# Get
setGeneric(name = "getHessian",
           def = function(theObject)
           {
             standardGeneric("getHessian")
           }
)

setMethod(f = "getHessian", 
          signature = "embeddedHMMModel", 
          definition = function(theObject)
          {
            return(theObject@Hessian)
          }
)

