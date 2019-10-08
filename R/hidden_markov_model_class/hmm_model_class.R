hiddenMarkovModel <- setClass(
  
  # Set class name
  "HiddenMarkovModel",
  
  # Set class properties
  slots = c(
    sigma = "numeric",
    epsilon = "numeric",
    theta = "numeric",
    p = "numeric",
    s = "numeric",
    t = "numeric",
    n_obs = "numeric",
    type = "character",
    transition_matrix = "matrix",
    observations = "list",
    hidden_states = "list",
    b = "numeric",
    semi_parametric_distribution = "matrix",
    llk = "vector",
    grid_search_parameters = "data.frame",
    hessian_grid_parameters = "data.frame"
  ),
  
  # Constructor
  prototype=list(
    sigma = 1,
    epsilon = 0.1,
    theta = NULL,
    p = 2/3,
    s = 0.0025,
    t = 0.0175,
    n_obs = 10000,
    type = "normal",
    transition_matrix = NULL,
    observations = NULL,
    hidden_states = NULL,
    b = NULL,
    semi_parametric_distribution = NULL,
    llk = c(),
    grid_search_parameters = data.frame(),
    hessian_grid_parameters = data.frame()
  ),
  validity = check_hidden_markov_model
)



