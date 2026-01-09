#' @name vmgeomVst
#' @aliases vmgeomVstFromMagn
#' @aliases vmgeomVstToR
#' @title Variance-Stabilizing Transformation for Geometric Visual Meteor Magnitudes
#'
#' @description
#' Applies a variance-stabilizing transformation to visual meteor magnitudes
#' under the geometric model.
#'
#' @param m integer; meteor magnitude.
#' @param lm numeric; limiting magnitude.
#' @param tm numeric; transformed magnitude.
#' @param log logical; if `TRUE`, the logarithm of the population index `r` is returned.
#' @param deriv.degree integer; the order of the derivative at `tm` to return
#'   instead of `r` or `log(r)`. Must be `0`, `1`, or `2`.
#'
#' @details
#' Many linear models require the variance of visual meteor magnitudes to be
#' homoscedastic. The function `vmgeomVstFromMagn` applies a transformation
#' that produces homoscedastic distributions of visual meteor magnitudes if the
#' underlying distribution follows a geometric model.
#'
#' The geometric model of visual meteor magnitudes
#' depends on the [population index][vismeteor::vmgeom] `r` and the limiting magnitude `lm`,
#' resulting in a two-parameter distribution. Without detection probabilities,
#' the magnitude distribution is purely geometric, and for integer limiting
#' magnitudes the variance depends only on the population index `r`. Since the
#' limiting magnitude `lm` is a fixed parameter and never estimated
#' statistically, the magnitudes can be transformed such that, for example,
#' the mean of the transformed magnitudes directly provides an estimate of `r`
#' using the function `vmgeomVstToR`.
#'
#' A key advantage of this transformation is that the limiting magnitude `lm`
#' is already incorporated into subsequent analyses. In this sense, the
#' transformation acts as a normalization of meteor magnitudes and yields a
#' variance close to `1.0`.
#'
#' This transformation is valid for \eqn{1.4 \le r \le 3.5}.
#' The numerical form of the transformation is version-specific and may change
#' substantially in future releases. Do not rely on equality of transformed
#' values across package versions.
#'
#' @return
#' * `vmgeomVstFromMagn`: numeric value, the transformed meteor magnitude.
#' * `vmgeomVstToR`: numeric value of the population index `r`, derived from
#'   the mean of `tm`.
#'
#' The argument `deriv.degree` can be used to apply the delta method.
#' If `log = TRUE`, the logarithm of `r` is returned.
#'
#' @note
#' The internal approximations used here are derived from the perception
#' probabilities produced by [vismeteor::vmperception].
#' For details on the derivation, see the script `inst/derivation/vmgeom_vst.R` in the
#' package's source code.
#'
#' @seealso [vismeteor::vmgeom] [vismeteor::vmperception]
#' @examples
#' N <- 100
#' r <- 2.0
#' limmag <- 6.3
#'
#' # Simulate magnitudes
#' m <- rvmgeom(N, limmag, r)
#'
#' # Variance-stabilizing transformation
#' tm <- vmgeomVstFromMagn(m, limmag)
#' tm.mean <- mean(tm)
#' tm.var  <- var(tm)
#'
#' # Estimator for r from the transformed mean
#' r.hat  <- vmgeomVstToR(tm.mean)
#'
#' # Derivative dr/d(tm) at tm.mean (needed for the delta method)
#' dr_dtm <- vmgeomVstToR(tm.mean, deriv.degree = 1L)
#'
#' # Variance of the sample mean of tm
#' var_tm.mean <- tm.var / N
#'
#' # Delta method: variance and standard error of r.hat
#' var_r.hat <- (dr_dtm^2) * var_tm.mean
#' se_r.hat  <- sqrt(var_r.hat)
#'
#' # Results
#' print(r.hat)
#' print(se_r.hat)

#' @rdname vmgeomVst
#' @export
vmgeomVstFromMagn <- function(m, lm) {
    offset <- lm - round(lm)
    if (1L == length(lm)) {
        limmag <- rep(lm, length(m))
        offset <- rep(offset, length(m))
    }

    a <- .vmgeomVstFromMagn.params$pa.fun(offset)
    b <- .vmgeomVstFromMagn.params$pb.fun(offset)
    c <- .vmgeomVstFromMagn.params$pc.fun(offset)
    d <- .vmgeomVstFromMagn.params$pd.fun(offset)

    x <- lm - m
    a - exp(b + c * (x + 0.5)^d)
}

#' @rdname vmgeomVst
#' @export
vmgeomVstToR <- function(tm, log = FALSE, deriv.degree = 0L) {
    poly.coef0 <- c(10.31637299, -5.51610811, 1.38003791, -0.18014183, 0.00946975)
    names(poly.coef0) <- seq_along(poly.coef0) - 1 # exponents

    # tm min 3.96 (r approx 3.5)
    # tm max 5.74 (r approx 1.4)
    tm[tm < 3.96 | tm > 5.74] <- NA

    if(deriv.degree > 0L) {
        poly.coef1 <- f.polynomial.coef(poly.coef0, deriv.degree = 1L)
    }
    if(deriv.degree > 1L) {
        poly.coef2 <- f.polynomial.coef(poly.coef1, deriv.degree = 1L)
    }
    if(deriv.degree > 2L) {
       stop(paste('deriv.degree', deriv.degree, 'not implemented!'))
    }

    if(log) {
        if (2L == deriv.degree) {
            f.polynomial(tm, poly.coef2)
        } else if (1L == deriv.degree) {
            f.polynomial(tm, poly.coef1)
        } else {
            f.polynomial(tm, poly.coef0)
        }
    } else {
        if(2L == deriv.degree) {
            exp(f.polynomial(tm, poly.coef0)) * (
               f.polynomial(tm, poly.coef1)^2 + f.polynomial(tm, poly.coef2)
            )
        } else if(1L == deriv.degree) {
            f.polynomial(tm, poly.coef1) * exp(f.polynomial(tm, poly.coef0))
        } else {
            exp(f.polynomial(tm, poly.coef0))
        }
    }
}

#' @keywords internal
.vmgeomVstFromMagn.params <- (function() {
    param.df <- data.frame(
        offset = c(-0.5, -0.48, -0.45, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.48, 0.5),
        a = c(8.810993, 8.496702, 8.29026, 8.140408, 8.031893, 8.024524, 8.046554, 8.093328, 8.165077, 8.265932, 8.384126, 8.564082, 8.678373, 8.757368, 8.811488),
        b = c(2.176, 2.077, 2.004, 1.945, 1.897, 1.885, 1.889, 1.906, 1.933, 1.97, 2.02, 2.087, 2.128, 2.156, 2.176),
        c = c(-0.317, -0.278, -0.249, -0.225, -0.206, -0.2, -0.201, -0.207, -0.217, -0.231, -0.252, -0.279, -0.296, -0.308, -0.317),
        d = c(0.669, 0.748, 0.813, 0.872, 0.922, 0.934, 0.928, 0.91, 0.881, 0.843, 0.796, 0.739, 0.706, 0.684, 0.669)
    )

    list(
        pa.fun = stats::approxfun(param.df$offset, param.df$a),
        pb.fun = stats::approxfun(param.df$offset, param.df$b),
        pc.fun = stats::approxfun(param.df$offset, param.df$c),
        pd.fun = stats::approxfun(param.df$offset, param.df$d)
    )
})()
