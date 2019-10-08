#' @title Mock data grid search
#' @description Apply grid search method over embedded hmm
#' @return results of grid search computation

# Set

setGeneric(name = "MockDataGridSearchMethod",
           def = function(theObject,sigma,epsilon,p,type, constraint, fine_grid = FALSE,config_file)
           {
             standardGeneric("MockDataGridSearchMethod")
           }
)

setMethod(f = "MockDataGridSearchMethod", 
          signature = "embeddedHMMModel", 
          definition = function(theObject,sigma,epsilon,p,type, constraint,fine_grid = FALSE,config_file)
          {
            grid_search_results <- 
              MockDataGridSearch(sigma, epsilon, p, type, constraint = constraint,fine_grid = fine_grid,
                                 config_file = config_file)
            theObject@Grid_search_result <- grid_search_results
            return(theObject)          
          }
)

# Get

# HMM
setGeneric(name = "getHMM",
           def = function(theObject)
           {
             standardGeneric("getHMM")
           }
)

setMethod(f = "getHMM", 
          signature = "embeddedHMMModel", 
          definition = function(theObject)
          {
            return(theObject@Grid_search_result$HMM$model)
          }
)

# Grid
setGeneric(name = "getGrid",
           def = function(theObject)
           {
             standardGeneric("getGrid")
           }
)

setMethod(f = "getGrid", 
          signature = "embeddedHMMModel", 
          definition = function(theObject)
          {
            return(theObject@Grid_search_result$grid)
          }
)

# llk
setGeneric(name = "getLlk",
           def = function(theObject)
           {
             standardGeneric("getLlk")
           }
)

setMethod(f = "getLlk", 
          signature = "embeddedHMMModel", 
          definition = function(theObject)
          {
            return(theObject@Grid_search_result$llk)
          }
)


