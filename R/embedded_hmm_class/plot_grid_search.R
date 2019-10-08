#' @title plot grid search
#' @description Plot the manifold of the grid search
#' @return A 3D plot of the p-theta-llk surface

# Set

setGeneric(name = "PlotLlkSurface",
           def = function(theObject)
           {
             standardGeneric("PlotLlkSurface")
           }
)

setMethod(f = "PlotLlkSurface", 
          signature = "embeddedHMMModel", 
          definition = function(theObject)
          {
            
            dat = getGrid(theObject)
            dat$llk <- unlist(getLlk(theObject))
            
            # need to pivot the grid data frame
            p <- unique(dat$p)
            theta <- unique(dat$theta)
            llk_matrix = matrix(dat$llk, ncol = length(theta),
                                nrow = length(P), byrow = T)
            surface_plot <- plot_ly(x = theta, y = p , z = llk_matrix) %>% 
              add_surface()
            surface_plot %>% 
              layout(
                title = "Grid parameters surface",
                scene = list(
                  xaxis = list(title = "Theta"),
                  yaxis = list(title = "P"),
                  zaxis = list(title = "Log-Likelihood")
                ))
          }
)
