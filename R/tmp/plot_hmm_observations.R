# Set
library(ggplot2)
setGeneric(name = "plotHMMObservations",
           def = function(theObject, conditional = TRUE)
           {
             standardGeneric("plotHMMObservations")
           }
)

setMethod(f = "plotHMMObservations", 
          signature = "HiddenMarkovModel", 
          definition = function(theObject, conditional = TRUE)
          {
            if(conditional)
            {
              tryCatch(
                {
                  plot_data <- data.frame(obs = theObject@observations$observations,
                                          hidden = theObject@observations$hidden_states)
                  histogram <- ggplot(plot_data, aes(obs)) + 
                    facet_wrap(hidden ~ .) +
                    geom_histogram(binwidth = 0.1)  
                  return(histogram)
                },
                error=function(cond)
                {
                  message("There's a problem with plotting the histgram")
                  message(cond)
                  return(theObject)
                })  
            }
            else
            {
            tryCatch(
              {
                plot_data <- data.frame(obs = theObject@observations$observations,
                                        hidden = theObject@observations$hidden_states)
                histogram <- ggplot(plot_data, aes(obs, fill = factor(hidden))) + 
                  geom_histogram(bins = 50)
                return(histogram)
              },
              error=function(cond)
              {
                message("There's a problem with plotting the histgram")
                message(cond)
                return(theObject)
              }
            )
            }
          })

