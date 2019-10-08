#' @title Create parameter grid for hessian computation
#' @description Creates parameters grid (cartesian product)
#' @param hmm - HiddenMarkovModel object
#' @return Grid parameter for Hessian matrix

# Set

setGeneric(name = "createParametersGridForHMM",
           def = function(theObject, seq_out = T, conf_file)
           {
             standardGeneric("createParametersGridForHMM")
           }
)

setMethod(f = "createParametersGridForHMM", 
          signature = "HiddenMarkovModel", 
          definition = function(theObject, seq_out = T, conf_file)
          {
            observation <- getObservationsFromHMM(theObject)
            # Create Grid
            if(seq_out)
            {
              theta <- getTheta(theObject)
              theta_step_size = conf_file$thetaStepSize
              theta_grid <- seq(max(1e-2,theta - theta_step_size), min(theta + theta_step_size,1), by = 0.05)
              p <- theObject@p
              p_step_size <- conf_file$pStepSize
              p_grid <- seq(max(0,p - p_step_size), min(p + p_step_size,1), by = 0.05)
              parameters_grid <- expand.grid(theta_grid,p_grid)
            }
            else
              # For numeric gradient computation
            {
              mle_theta <- getTheta(theObject)
              gradient_step <- as.numeric(conf_file$gridStepForGradient)
              max_theta <- mle_theta + gradient_step
              min_theta <- max(mle_theta - gradient_step, 0)
              theta_grid <- c(min_theta,mle_theta,max_theta)
              mle_p <- theObject@p
              max_p <- min(mle_p + gradient_step,1)
              min_p <- max(mle_p - gradient_step, 1e-3)
              p_grid <- c(min_p,mle_p,max_p)
              parameters_grid <- expand.grid(theta_grid,p_grid)  
            }
            # Apply NR to compute b         
            # b <-mapply(function(x,y)NewtonRaphson(observation, x, y), parameters_grid[,1], parameters_grid[,2])
            data_for_central_ditribution <- data.frame(theta = parameters_grid[,1], p = parameters_grid[,2])
            theObject@grid_search_parameters = data_for_central_ditribution
            return(theObject)
          })

# Get

setGeneric(name = "getGridSearchParameters",
           def = function(theObject)
           {
             standardGeneric("getGridSearchParameters")
           }
)

setMethod(f = "getGridSearchParameters", 
          signature = "HiddenMarkovModel", 
          definition = function(theObject)
          {
            return(theObject@grid_search_parameters)
          })

