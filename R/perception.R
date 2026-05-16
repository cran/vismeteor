#' @title Perception Probabilities of Visual Meteor Magnitudes
#' @description
#' Provides the perception probability of visual meteor magnitudes.
#' @param dm numeric; difference between the limiting magnitude and the meteor magnitude.
#' @details
#' The perception probabilities of _Koschack R., Rendtel J., 1990b_
#' are estimated with the formula
#' \deqn{
#'     p(dm) = \begin{cases}
#'         1.0 - \exp\left(-z(dm + 0.5)\right)\  & \text{ if } dm > -0.5,\\
#'         0.0 \  & \text{ otherwise,}
#'     \end{cases}
#' }
#' where
#' \deqn{
#' z(x) = 0.0037 \, x + 0.0019 \, x^2 + 0.00271 \, x^3 + 0.0009 \, x^4
#' }
#' and `dm` is the difference between the limiting magnitude and the meteor magnitude.
#' @return This function returns the visual perception probabilities.
#' @references
#' Koschack R., Rendtel J., 1990b
#' _Determination of spatial number density and mass index from visual meteor observations (II)._
#' WGN 18, 119–140.
#' @examples
#' # Perception probability of visually estimated meteor of magnitude 3.0
#' # with a limiting magnitude of 5.6.
#' vmperception(5.6 - 3.0)
#'
#' # plot
#' old_par <- par(mfrow = c(1, 1))
#' plot(
#'     vmperception,
#'     -0.5, 8,
#'     main = paste(
#'         "perception probability of",
#'         "visual meteor magnitudes"
#'     ),
#'     col = "blue",
#'     xlab = "dm",
#'     ylab = "p"
#' )
#'
#' par(old_par)
#' @export
vmperception <- function(dm) {
    poly_coef <- c(0.0, 0.0037, 0.0019, 0.00271, 0.0009)
    names(poly_coef) <- seq(along = poly_coef) - 1 # exponents

    dm <- dm + 0.5
    p <- rep(0.0, length(dm))
    idx <- dm > .Machine$double.eps
    if (any(idx)) {
        f0 <- .f_polynomial(dm[idx], poly_coef)
        p[idx] <- 1.0 - exp(-f0)
    }

    p
}

#' build polynomial sum
#'
#' @noRd
.f_polynomial <- function(m, poly_coef) {
    exponents <- as.numeric(names(poly_coef))
    margin.table(poly_coef * t(outer(m, exponents, "^")), 2)
}

#' returns polynomial coefficients
#'
#' @noRd
.f_polynomial_coef <- function(poly_coef, deriv_degree = 1L) {
    if (0L == deriv_degree) {
        return(poly_coef)
    }

    if (1 == length(poly_coef)) {
        return(0)
    }

    exponents <- as.numeric(names(poly_coef))
    poly_coef <- poly_coef * exponents
    if (0L %in% exponents) {
        intercept_idx <- 0L == exponents
        poly_coef <- poly_coef[!intercept_idx]
        exponents <- exponents[!intercept_idx]
    }
    exponents <- exponents - 1L
    names(poly_coef) <- exponents

    .f_polynomial_coef(poly_coef, deriv_degree - 1L)
}
