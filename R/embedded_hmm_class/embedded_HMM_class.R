#' An S4 class to represent a HMM with exponential perturbations.
#'
#' @slot balance A length-one numeric vector

embeddedHMM <- setClass(

  # Set class name
  "embeddedHMMModel",

  # Set class properties
  slots = c(
    Grid_search_result = "list",
    MLE = "list",
    Gradient = "vector",
    Hessian = "list",
    CI = "data.frame"
  ),

  # Constructor
  prototype=list(
    Grid_search_result = list(),
    MLE = list(),
    Gradient = list(),
    Hessian = list(),
    CI = data.frame()
  ),

  validity = validate_embedded_hmm
)



