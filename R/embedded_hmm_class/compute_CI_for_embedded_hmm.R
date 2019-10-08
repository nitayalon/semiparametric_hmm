#' @title compute CI for embedded HMM
#' @description compute confidance interval for embedded HMM parameters
#' @return CI for MLE 

# Set
setGeneric(name = "computeCIForParameters",
           def = function(theObject,n)
           {
             standardGeneric("computeCIForParameters")
           }
)

setMethod(f = "computeCIForParameters", 
          signature = "embeddedHMMModel", 
          definition = function(theObject,n)
          {
            Hessian <- getHessian(theObject)
            if(length(Hessian) == 0)
            {
              cat("You must compute the Hessian first")
              return(theObject)
            }
            MLE <- getMLE(theObject)
            fisher_matrix <- Hessian$fisher_information
            variance <- sqrt(1/diag(fisher_matrix))
            upper <- MLE + qnorm(0.975) * variance
            lower <- MLE - qnorm(0.975) * variance
            CI <- data.frame(MLE = MLE, lower = lower, upper = upper)
            theObject@CI = CI  
            return(theObject)          
          }
)

# Get
setGeneric(name = "getCI",
           def = function(theObject)
           {
             standardGeneric("getCI")
           }
)

setMethod(f = "getCI", 
          signature = "embeddedHMMModel", 
          definition = function(theObject)
          {
            return(theObject@CI)
          }
)

setGeneric(name = "getMLE",
           def = function(theObject)
           {
             standardGeneric("getMLE")
           }
)

setMethod(f = "getMLE", 
          signature = "embeddedHMMModel", 
          definition = function(theObject)
          {
            HMM_model <- getHMM(theObject)
            # MLE's
            mle_theta <- getTheta(HMM_model)
            mle_p <- HMM_model@p
            mle_t <- HMM_model@transition_matrix[1,1]
            mle_s <- HMM_model@transition_matrix[1,2]
            mle_q <- HMM_model@transition_matrix[2,2]
            return(c(mle_theta,mle_p,mle_t,mle_s,mle_q))
          }
)
