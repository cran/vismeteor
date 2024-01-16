#' @title Perception Probabilities of Visual Meteor Magnitudes
#' @description
#' Provides the perception probability of visual meteor magnitudes and its first derivative.
#' @param m numerical; difference between the limiting magnitude and the meteor magnitude.
#' @param deriv.degree integer; degree of derivative of the perception probability.
#' Currently, valid values of `deriv.degree` are `0`, `1` and `2`.
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
#' z(x) = 0.003 \, x + 0.0056 \, x^2 + 0.0014 \, x^4
#' }
#' and `m` is the difference between the limiting magnitude and the meteor magnitude.
#' @return This function returns the visual perception probabilities.
#' If `deriv.degree` is specified, it will return the `deriv.degree`-th order derivative
#' of the perception probability.
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
#' plot(
#'     function(m) {
#'         vmperception(m, deriv.degree=1L)/vmperception(m)
#'     },
#'     -0.3, 8,
#'     main = paste(
#'         'q-values of',
#'         'visual meteor magnitudes'
#'     ),
#'     col = "blue",
#'     log = 'y',
#'     xlab = 'm',
#'     ylab = 'q'
#' )
#'
#' par(old_par)
#' @export
vmperception <- function(m, deriv.degree = 0L) {
    poly.coef <- c(0.0, 0.003, 0.0056, 0, 0.0014)
    names(poly.coef) <- seq(along = poly.coef) - 1 # exponents

    m <- m + 0.5
    p <- rep(0.0, length(m))
    idx <- m > .Machine$double.eps
    if (any(idx)) {
        f0 <- f.polynomial(m[idx], poly.coef)
        if (0L == deriv.degree) {
            # 1 - exp(-f(x))
            p[idx] <- 1.0 - exp(-f0)
        } else if (1L == deriv.degree) {
            # f'(x) * exp(-f(x))
            poly.coef1 <- f.polynomial.coef(poly.coef, deriv.degree = 1L)
            f1 <- f.polynomial(m[idx], poly.coef1)
            p[idx] <- exp(-f0) * f1
        } else if (2L == deriv.degree) {
            # (f''(x) -f'(x)^2) * exp(-f(x))
            poly.coef1 <- f.polynomial.coef(poly.coef, deriv.degree = 1L)
            f1 <- f.polynomial(m[idx], poly.coef1)
            poly.coef2 <- f.polynomial.coef(poly.coef, deriv.degree = 2L)
            f2 <- f.polynomial(m[idx], poly.coef2)
            p[idx] <- exp(-f0) * (f2 - f1^2)
        } else {
            stop(paste('deriv.degree', deriv.degree, 'not implemented!'))
        }
    }

    p
}

#' @title Laplace-Transformed Perception Probabilities of Visual Meteor Magnitudes
#' @description
#' Provides the Laplace-transformed perception probability of visual meteor magnitudes
#' and its first derivative.
#' @param s numerical; Real (non-complex) parameter for the Laplace transformation.
#' @param deriv.degree integer; degree of derivative of the transformation.
#' Currently, valid values of `deriv.degree` are `0`, `1` and `2`.
#' @details
#' The Laplace-transformed [perception probabilities][vismeteor::vmperception] `F(s)`, given as
#' \deqn{
#' F(s) = \mathcal{L} \left\{p\right\}(s)
#' = \int_{-0.5}^{\infty} \, f(m) \, \mathrm e^{-s \, m} \,\mathrm{d}m \,,
#' }
#' are approximately
#' \deqn{
#'     P(s) = \begin{cases}
#'         s^{-1} \, \exp\left(-4.11 \, s + 1.32 \, s^2 - 0.15 \, s^3\right)\ & \text{ if } s >= 0.0,\\
#'         \text{undefined} \  & \text{ otherwise.}
#'     \end{cases}
#' }
#' Here, `m` is the difference between the limiting magnitude and the meteor magnitude,
#' and `f(m)` denotes the perception probabilities as a function of `m`.
#' The \eqn{\mathcal{L}} recalls here the one-sided Laplace transform.
#'
#' The Laplace transform is notably effective for determining the mean and variance
#' of observed meteor magnitudes, which are measured relative to the limiting magnitude.
#' This is just one example of its application.
#' This approach is valid only when the actual magnitude distribution adheres
#' to \eqn{p(m) \sim r^{-m}}, where \eqn{s = \log(r)}.
#' In this scenario, the mean of the observable meteor magnitudes is given by
#' \eqn{-\mathcal{L}'/\mathcal{L}}, and their variance is calculated as
#' \eqn{\mathcal{L}''/\mathcal{L} - (\mathcal{L}'/\mathcal{L})^2}.
#' @return returns the Laplace-transformed perception probabilities.
#' If `deriv.degree` is specified, it will return the `deriv.degree`-th order derivative
#' of these Laplace-transformed values.
#' @seealso
#'   [vismeteor::vmperception]
#'   [vismeteor::vmgeom]
#' @examples
#' r <- 2.0
#' s <- log(r)
#' F0 <- vmperception.l(s)
#' F1 <- vmperception.l(s, deriv.degree=1L)
#' # magnitude mean
#' -F1/F0
#' F2 <- vmperception.l(s, deriv.degree=2L)
#' # magnitude variance
#' F2/F0 - (F1/F0)^2
#' # plot the Laplace-transformed perception probabilities
#' old_par <- par(mfrow = c(1,1))
#' plot(
#'     vmperception.l,
#'     0.2, 1.1,
#'     main = paste(
#'         'Laplace-transformed perception',
#'         'probability of visual meteor magnitudes'
#'     ),
#'     col = "blue",
#'     log = 'y',
#'     xlab = 's',
#'     ylab = 'L'
#' )
#' par(old_par)
#' @export
vmperception.l <- function(s, deriv.degree = 0L) {
    poly.coef <- c(0.0, -4.11, 1.32, -0.15)
    names(poly.coef) <- seq(along = poly.coef) - 1 # exponents

    f0 <- f.polynomial(s, poly.coef)
    if (0L == deriv.degree) {
        # exp(f(s))/s
        exp(f0)/s
    } else if (1L == deriv.degree) {
        # e^f(s) ((f'(s))/s - 1/s^2)
        poly.coef1 <- f.polynomial.coef(poly.coef, deriv.degree = 1L)
        f1 <- f.polynomial(s, poly.coef1)
        exp(f0) * (f1/s - s^-2)
    } else if (2L == deriv.degree) {
        # e^f(s) ((f''(s))/s - (2 f'(s))/s^2 + f'(s)^2/s + 2/s^3)
        poly.coef1 <- f.polynomial.coef(poly.coef, deriv.degree = 1L)
        f1 <- f.polynomial(s, poly.coef1)
        poly.coef2 <- f.polynomial.coef(poly.coef, deriv.degree = 2L)
        f2 <- f.polynomial(s, poly.coef2)
        exp(f0) * ( f2/s - 2*f1*s^-2 + (f1^2)/s + 2*(s^-3))
    } else {
        stop(paste('deriv.degree', deriv.degree, 'not implemented!'))
    }
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
