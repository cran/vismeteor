#' @title Perception Probabilities of Visual Meteor Magnitudes
#' @description
#' Provides the perception probability of visual meteor magnitudes.
#' @param m numeric; difference between the limiting magnitude and the meteor magnitude.
#' @details
#' The perception probabilities of _Koschack R., Rendtel J., 1990b_
#' are estimated with the formula
#' \deqn{
#'     p(m) = \begin{cases}
#'         1.0 - \exp\left(-z(m + 0.5)\right)\  & \text{ if } m > -0.5,\\
#'         0.0 \  & \text{ otherwise,}
#'     \end{cases}
#' }
#' where
#' \deqn{
#' z(x) = 0.0037 \, x + 0.0019 \, x^2 + 0.00271 \, x^3 + 0.0009 \, x^4
#' }
#' and `m` is the difference between the limiting magnitude and the meteor magnitude.
#' @return This function returns the visual perception probabilities.
#' @references Koschack R., Rendtel J., 1990b _Determination of spatial number density and mass index from visual meteor observations (II)._ WGN 18, 119–140.
#' @examples
#' # Perception probability of visually estimated meteor of magnitude 3.0
#' # with a limiting magnitude of 5.6.
#' vmperception(5.6 - 3.0)
#'
#' # plot
#' old_par <- par(mfrow = c(1,1))
#' plot(
#'     vmperception,
#'     -0.5, 8,
#'     main = paste(
#'         'perception probability of',
#'         'visual meteor magnitudes'
#'     ),
#'     col = "blue",
#'     xlab = 'm',
#'     ylab = 'p'
#' )
#'
#' par(old_par)
#' @export
vmperception <- function(m) {
    poly.coef <- c(0.0, 0.0037, 0.0019, 0.00271, 0.0009)
    names(poly.coef) <- seq(along = poly.coef) - 1 # exponents

    m <- m + 0.5
    p <- rep(0.0, length(m))
    idx <- m > .Machine$double.eps
    if (any(idx)) {
        f0 <- f.polynomial(m[idx], poly.coef)
        p[idx] <- 1.0 - exp(-f0)
    }

    p
}

#' build polynomial sum
#'
#' @noRd
f.polynomial <- function(m, poly.coef) {
    exponents <- as.numeric(names(poly.coef))
    margin.table(poly.coef * t(outer(m, exponents, "^")), 2)
}

#' returns polynomial coefficients
#'
#' @noRd
f.polynomial.coef <- function(poly.coef, deriv.degree = 1L) {
    if (0L == deriv.degree)
        return(poly.coef)

    if (1 == length(poly.coef))
        return(0)

    exponents <- as.numeric(names(poly.coef))
    poly.coef <- poly.coef * exponents
    if (0L %in% exponents) {
        intercept.idx <- 0L == exponents
        poly.coef <- poly.coef[! intercept.idx]
        exponents <- exponents[! intercept.idx]
    }
    exponents <- exponents - 1L
    names(poly.coef) <- exponents

    f.polynomial.coef(poly.coef, deriv.degree - 1L)
}
