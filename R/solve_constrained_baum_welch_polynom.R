#' @title Computes polynom coefficients for S and T
#' @description Converts BW counts to polynom coefficients
#' @param non_normalized_transition_matrix - Matrix, count matrix from BW iteration
#' @param p - Numeric, the stationary distribution parameter
#' @return List of polynom coefficients

computePolynomCoefficients <- function(non_normalized_transition_matrix, p)
{
  stopifnot(class(non_normalized_transition_matrix) == "matrix")
  if(any(non_normalized_transition_matrix < 0))
  {
    warning("One of the entries is zero")
    non_normalized_transition_matrix[non_normalized_transition_matrix < 0] <- 1
  }
  alpha = p / (1-p)

  V = non_normalized_transition_matrix
  A1 = V[1,2] + V[3,2] + V[2,1] + V[2,3]
  A2 = V[1,1] + V[3,3]
  A3 = V[2,2]
  B1 = V[3,1] + V[1,3]
  B2 = V[1,1] + V[3,3]

  return(list(A1 = A1,
              A2 = A2,
              A3 = A3,
              B1 = B1,
              B2 = B2,
              alpha = alpha))
}

#' @title Solves polynom for S
#' @description Compute solutions to quadratic equation for S
#' @param polynom_coefficients - List, polynom coefficients
#' @return vector of solution for s

solveForS <- function(polynom_coefficients)
{
  A1 = polynom_coefficients$A1
  A2 = polynom_coefficients$A2
  A3 = polynom_coefficients$A3
  B1 = polynom_coefficients$B1
  B2 = polynom_coefficients$B2
  alpha = polynom_coefficients$alpha

  upper_limit <- min(1,1/alpha)
  P1 = A1 + A2 + A3 - (A1*B1 + A3*B1) / (B1 + B2)
  P2 = A1/alpha - A1 + A1*B1 / (alpha * (B1+B2)) + A1*B1 / (B1+B2) -
    A2/alpha - A3 + A1*B1 / (B1+B2)
  P3 = A1/alpha - A1*B1*(1/(alpha*(B1+B2)) + 1/((B1+B2)))

  discreminante <- P2^2 - 4 * P1 * P3
  if(discreminante < 0)
  {
    return(c(0))
  }

  solutions_1 <- (-P2 - sqrt(P2^2 - 4 * P1 * P3)) / (2 * P1)
  solutions_2 <- (-P2 + sqrt(P2^2 - 4 * P1 * P3)) / (2 * P1)

  solution_vector <- c(solutions_1,solutions_2)

  if(!all(is.numeric(solution_vector)))
  {
    return(0)
  }
  if(all(solution_vector < 0))
  {
    return(0)
  }
  if(all(solution_vector > upper_limit))
  {
    return(upper_limit)
  }
  if(any(solution_vector < upper_limit) &
     any(solution_vector > 0))
  {
    return(findFeasibleSolution(solution_vector,upper_limit))
  }
  return(0.5)
}

findFeasibleSolution <- function(solution_vector, upper_limit)
{
  if(all(solution_vector > 0) & all(solution_vector < upper_limit))
  {
    return(solution_vector)
  }
  feasible_value <- solution_vector[which(solution_vector > 0 & solution_vector < upper_limit)]
  if(length(feasible_value) == 0){return(upper_limit)}
  return(feasible_value)
}

#' @title Solves polynom for T
#' @description Compute solutions to linear equation for T
#' @param feasible_s_value - List, feasible solutions for S
#' @param polynom_coefficients - List, polynom coefficients
#' @return vector of solution for t

solveForT <- function(feasible_s_value, polynom_coefficients)
{
  B1 <- polynom_coefficients$B1
  B2 <- polynom_coefficients$B2
  epsilon = 0.01
  B1 <- B1 + epsilon
  P4 <- B1 / (B1 + B2)
  current_value <- (1-feasible_s_value) * P4
  if(all(current_value <= 0))
  {
    return(0)
  }
  if(all(current_value >= 1))
  {
    return(1)
  }
  return(current_value[current_value > 0 && current_value < 1])
}
