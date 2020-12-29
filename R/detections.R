######################################
######################################
#### detection_pr()

#' @title A detection probability function based on distance
#' @description This function calculates detection probability (e.g., of an acoustic detection) at specified distances from the sampling device (e.g., a passive acoustic telemetry receiver) using user-defined parameters (a model intercept, an coefficient for the effect of distance and an inverse link function). The function returns a plot of detection probability with distance and/or a vector of detection probabilities.
#'
#' @param distance A numeric vector of distances at which to calculate detection probability.
#' @param beta_0,beta_1 Single numbers that define the model coefficients (i.e., the intercept and gradient on the scale of the link function).
#' @param inv_link A function that defines the inverse link function.The default function is the logistic (inverse logit) function.
#' @param output An integer (\code{1L}, \code{2L} or \code{3L}) that defines the output type. \code{1L} returns a plot of detection probability against distance; \code{2L} returns a numeric vector of detection probabilities; and \code{3L} returns both of the above.
#' @param ... Additional arguments passed to \code{\link[prettyGraphics]{pretty_plot}} to customise the plot. These are only implemented if \code{output = 1L} or \code{output = 3L}.
#'
#' @return The function calculates detection probability at each specified distance and returns a plot, a vector of detection probabilities, or both, depending on the value of the \code{output} argument. If a vector of detection probabilities is returned, this contains the following attributes: X', the model matrix; 'beta', the regression coefficients; and 'inv_link', the inverse link function.
#'
#' @examples
#' #### Example (1): Implement the function using the default parameters
#' # The function returns a graph and a vector of detection probabilities
#' det_pr <- detection_pr()
#' utils::head(det_pr)
#' # The vector has attributes:
#' # ... 'X' (the model matrix),
#' # ... 'beta' (the regression coefficient)
#' # ... 'inv_link' (the inverse link function)
#' utils::str(det_pr)
#'
#' #### Example (2): Adjust model parameters
#' # Change regression coefficients
#' det_pr <- detection_pr(beta_0 = 2.5, beta_1 = -0.006)
#' # Use inverse probit
#' det_pr <- detection_pr(beta_0 = 2.5, beta_1 = -0.006, inv_link = stats::pnorm)
#'
#' #### Example (3): Modify graphical properties
#' det_pr <- detection_pr(beta_0 = 2.5,
#'                           beta_1 = -0.006,
#'                           type = "l",
#'                           xlab = "Distance (m)",
#'                           ylab = "Detection Probability")
#'
#' #### Example (4): Modify return options
#' # Only graph
#' detection_pr(output = 1L)
#' # Only values
#' detection_pr(output = 2L)
#' # Both graph and values (the default)
#' detection_pr(output = 3L)
#'
#' @author Edward Lavender
#' @export
#'

detection_pr <- function(distance = 1:1000,
                         beta_0 = 2.5,
                         beta_1 = -0.01,
                         inv_link = stats::plogis,
                         output = 3L,...){
  #### Checks
  stopifnot(length(beta_0) == 1 & length(beta_1) == 1)
  output <- check_value(input = output, supp = 1:3, warn = TRUE, default = 3L)

  #### Calculate detection probabilities
  X <- matrix(c(rep(1, length(distance)), distance), ncol = 2, byrow = FALSE)
  beta <- matrix(c(beta_0, beta_1), nrow = 2)
  rownames(beta) <- c("beta_0", "beta_1")
  Y <- inv_link(X %*% beta)
  Y <- as.numeric(Y)
  # Add attributes
  attributes(Y)$X        <- X
  attributes(Y)$beta     <- beta
  attributes(Y)$inv_link <- inv_link

  #### Visualise detection probabilities
  if(output %in% c(1, 3)){
    prettyGraphics::pretty_plot(X[, 2], Y,...)
  }

  #### Return detection probabilities
  if(output %in% c(2, 3)){
    return(Y)
  } else return(invisible())
}


