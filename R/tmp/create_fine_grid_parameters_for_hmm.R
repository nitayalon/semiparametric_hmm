#' @title Create parameter grid for hessian computation
#' @description Creates parameters grid (cartesian product)
#' @param hmm - HiddenMarkovModel object
#' @return Grid parameter for Hessian matrix

# Set

setGeneric(name = "createFineParametersGridForHMM",
           def = function(theObject, seq_out = T,conf_file)
           {
             standardGeneric("createFineParametersGridForHMM")
           }
)

setMethod(f = "createFineParametersGridForHMM", 
          signature = "HiddenMarkovModel", 
          definition = function(theObject, seq_out = T,conf_file)
          {
            observation <- getObservationsFromHMM(theObject)
            grid_width <- conf_file$fineGridWidth
            step_size <- conf_file$fineGridStepSize
            # Create Grid
            if(seq_out)
            {
              theta <- getTheta(theObject)
              theta_grid <- seq(max(0,theta - grid_width), min(theta + grid_width,1), by = step_size)
              p <- theObject@p
              p_grid <- seq(max(0,p - grid_width), min(p + grid_width,1), by = step_size)
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

setGeneric(name = "getFineGridSearchParameters",
           def = function(theObject)
           {
             standardGeneric("getFineGridSearchParameters")
           }
)

setMethod(f = "getFineGridSearchParameters", 
          signature = "HiddenMarkovModel", 
          definition = function(theObject)
          {
            return(theObject@grid_search_parameters)
          })

