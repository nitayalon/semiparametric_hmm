#' @title Create parameter grid for hessian computation
#' @description Creates parameters grid (cartesian product)
#' @param hmm - HiddenMarkovModel object
#' @return Grid parameter for Hessian matrix

# Set

setGeneric(name = "createParametersGridForHessian",
           def = function(theObject)
           {
             standardGeneric("createParametersGridForHessian")
           }
)

setMethod(f = "createParametersGridForHessian", 
          signature = "HiddenMarkovModel", 
          definition = function(theObject)
          {
            conf_file <- config::get(file = "~/mastersdegree/Thesis/Main/Code/Mock_data_classes/hidden_markov_model_class/hmm_config.yml")
            observation <- getObservationsFromHMM(theObject)
            hessian_step <- as.numeric(conf_file$gridStepForHessian)
            # MLE's
            mle_theta <- getTheta(theObject)
            mle_p <- theObject@p
            mle_t <- theObject@transition_matrix[1,1]
            mle_s <- theObject@transition_matrix[1,2]
            mle_q <- theObject@transition_matrix[2,2]
            # Create grid 
            mle_grid <- data.frame(mle = c(mle_theta,mle_p,mle_t,mle_s,mle_q))
            mle_grid$down <- mle_grid$mle - hessian_step
            mle_grid$up <- mle_grid$mle + hessian_step
            theObject@hessian_grid_parameters = mle_grid
            return(theObject)
          })


# Get
setGeneric(name = "getGridParametersForHessian",
           def = function(theObject)
           {
             standardGeneric("getGridParametersForHessian")
           }
)

setMethod(f = "getGridParametersForHessian", 
          signature = "HiddenMarkovModel", 
          definition = function(theObject)
          {
            return(theObject@hessian_grid_parameters)
          })

