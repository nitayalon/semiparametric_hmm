#' An S4 class to represent a HMM with exponential perturbations.
#'
#' @slot sigma numeric, the data standard deviation
#' @slot epsilon numeric, the proportion of the mean wrt to std
#' @slot theta numeric, exponential perturbation size
#' @slot p numeric, stationary vector proportion
#' @slot s numeric, transition probability to non-perturbed state
#' @slot t numeric, transition probability between perturbed states
#' @slot n_obs numeric, length of sample
#' @slot type character, type of conditional distribution
#' @slot transition_matrix matrix, HMM transition matrix
#' @slot observations list, Observed data and Hidden data


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
    observations = "list"
  ),

  # Constructor
  prototype=list(
    sigma = 1,
    epsilon = 0.1,
    theta = NULL,
    p = 2/3,
    s = 1/400,
    t = 7/400,
    n_obs = 10000,
    type = "Normal",
    transition_matrix = NULL,
    observations = NULL
  ),
  validity = check_hidden_markov_model
)



