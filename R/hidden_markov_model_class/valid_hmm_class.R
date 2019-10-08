check_hidden_markov_model <- function(object)
{
  errors <- character()
  
  if((object@p < 0) || (object@p > 1))
  {
    msg <- "A mixture proportion must by between 0 and 1"
    errors <- c(errors,msg)
  }
  if(object@sigma < 0)
  {
    msg <- "Scale parameter must be positive"
    errors <- c(errors,msg)
  }
  if(!(object@type %in% c("normal","laplace")))
  {
    msg <- "Only ADL and normal distributions are supported"
    errors <- c(errors,msg)
  }
  if(object@n_obs <= 0)
  {
    msg <- "Number of observations must be positive and natural"
    errors <- c(errors,msg)
  }
  
  if(any(c(object@s,object@t) < 0) || any(c(object@s,object@t) > 1))
  {
    msg <- "All the transition matrix entries must be between 0 and 1"
    errors <- c(errors,msg)
  }
  
  if(object@s + object@t > 1)
  {
    msg <- "t and s must be less than 1!"
    errors <- c(errors,msg)
  }
  
  if(length(errors) == 0) TRUE else errors
}