#' @name vmgeom
#' @aliases dvmgeom
#' @aliases pvmgeom
#' @aliases qvmgeom
#' @aliases rvmgeom
#' @title Geometric Model of Visual Meteor Magnitudes
#' @description
#' Density, distribution function, quantile function, and random generation
#' for the geometric model of visual meteor magnitudes.
#' @param m integer; the meteor magnitude.
#' @param p numeric; probability.
#' @param lm numeric; limiting magnitude.
#' @param r numeric; the population index.
#' @param n numeric; count of meteor magnitudes.
#' @param log logical; if `TRUE`, probabilities p are given as `log(p)`.
#' @param lower.tail logical; if `TRUE` (default) probabilities are
#'     \eqn{P[M < m]}, otherwise, \eqn{P[M \ge m]}.
#' @param perception_fun function; perception probability function (optional).
#'     Default is [vismeteor::vmperception].
#' @details
#' In visual meteor observations, magnitudes are estimated as integer values.
#' Consequently, the distribution of observed magnitudes is discrete, and its
#' probability mass function is given by
#'
#' \deqn{
#' P[M = m] \sim
#' \begin{cases}
#'   f(m_{\mathrm{lim}} - m)\, r^m, & \text{if } m_{\mathrm{lim}} - m > -0.5,\\[5pt]
#'   0 & \text{otherwise,}
#' \end{cases}
#' }
#'
#' where \eqn{m_{\mathrm{lim}}} denotes the limiting (non-integer) magnitude of the observation,
#' and \eqn{m} the integer meteor magnitude.
#' The function \eqn{f(\cdot)} denotes the [perception probability function][vismeteor::vmperception].
#'
#' Thus, the distribution is the product of the perception probabilities and the
#' underlying [geometric distribution][stats::Geometric] of meteor magnitudes.
#' Therefore, the parameter `p` of the geometric distribution is given by `p = 1 - 1/r`.
#'
#' The parameter `lm` specifies the limiting magnitude for the meteor magnitude `m`.
#' `m` must be an integer meteor magnitude.
#' The length of the vector `lm` must either equal the length of the vector `m`,
#' or `lm` must be a scalar value.
#' In the case of `rvmgeom`, the length of the vector `lm` must equal `n`,
#' or `lm` must be a scalar value.
#'
#' If a different perception probability function `perception_fun` is provided,
#' it must have the signature `function(x)` and return the perception probability of
#' the difference `x` between the limiting magnitude and the meteor magnitude.
#' If `x >= 15.0`, the function `perception_fun` should return a perception probability of `1.0`.
#' The argument `perception_fun` is resolved using [match.fun].
#'
#' @return
#' - `dvmgeom`: density
#' - `pvmgeom`: distribution function
#' - `qvmgeom`: quantile function
#' - `rvmgeom`: random generation
#'
#' The length of the result is determined by `n` for `rvmgeom`, and by the maximum
#' of the lengths of the numeric vector arguments for the other functions.
#' All arguments are vectorized; standard R recycling rules apply.
#'
#' Since the distribution is discrete, `qvmgeom` and `rvmgeom` always return integer values.
#' `qvmgeom` may return `NaN` with a warning.
#' @seealso [vismeteor::vmperception]
#'   [stats::Geometric]
#' @examples
#' N <- 100
#' r <- 2.0
#' limmag <- 6.5
#' (m <- seq(6, -7))
#'
#' # discrete density of `N` meteor magnitudes
#' (freq <- round(N * dvmgeom(m, limmag, r)))
#'
#' # log likelihood function
#' lld <- function(r) {
#'     -sum(freq * dvmgeom(m, limmag, r, log = TRUE))
#' }
#'
#' # maximum likelihood estimation (MLE) of r
#' est <- optim(2, lld, method = "Brent", lower = 1.1, upper = 4)
#'
#' # estimations
#' est$par # mean of r
#'
#' # generate random meteor magnitudes
#' m <- rvmgeom(N, r, lm = limmag)
#'
#' # log likelihood function
#' llr <- function(r) {
#'     -sum(dvmgeom(m, limmag, r, log = TRUE))
#' }
#'
#' # maximum likelihood estimation (MLE) of r
#' est <- optim(2, llr, method = "Brent", lower = 1.1, upper = 4, hessian = TRUE)
#'
#' # estimations
#' est$par # mean of r
#' sqrt(1 / est$hessian[1][1]) # standard deviation of r
#'
#' m <- seq(6, -4, -1)
#' p <- vismeteor::dvmgeom(m, limmag, r)
#' barplot(
#'     p,
#'     names.arg = m,
#'     main = paste0("Density (r = ", r, ", limmag = ", limmag, ")"),
#'     col = "blue",
#'     xlab = "m",
#'     ylab = "p",
#'     border = "blue",
#'     space = 0.5
#' )
#' axis(side = 2, at = pretty(p))

#' @rdname vmgeom
#' @export
dvmgeom <- function(m, lm, r, log = FALSE, perception_fun = vmperception) {
    if (anyNA(m) || anyNA(lm) || anyNA(r)) {
        stop("NA's are not allowed!")
    }

    if (any(is.infinite(lm))) {
        stop("Infinite limiting magnitudes are not allowed!")
    }

    if (any(is.infinite(r))) {
        stop("Infinite r values are not allowed!")
    }

    if (any(1.0 > r)) {
        stop(paste0('r must be greater than 1.0 instead of "', r, '"!'))
    }

    is_wholenumber <- function(x, tol = .Machine$double.eps^0.5) all(is.infinite(x) | abs(x - round(x)) < tol)
    if (!is_wholenumber(m)) {
        stop("magnitudes must be integer values!")
    }

    perception_fun <- match.fun(perception_fun)

    p_geom <- 1.0 - 1.0 / r

    if (1 == length(r)) {
        r <- rep(r, length(m))
        p_geom <- rep(p_geom, length(m))
    }

    std_res <- .vmgeom_std(m, lm)
    m <- std_res$m
    offset <- std_res$offset

    f_density <- function(m, offset, p_geom) {
        m_max <- 15L
        idx <- m <= m_max
        d <- stats::dgeom(m, p_geom, log = TRUE)
        if (any(idx)) {
            d[idx] <- d[idx] + base::log(perception_fun(m[idx] + offset))
        }

        d - base::log(.vmgeom_norm(offset, p_geom, m_max, perception_fun))
    }

    arg_data <- data.frame(
        m = m,
        offset = offset,
        p_geom = p_geom
    )

    data_f <- as.factor(paste0(offset, "/", p_geom))
    data_s <- split(arg_data, data_f)
    d <- lapply(data_s, \(data) {
        m <- data$m
        offset <- data$offset[1]
        p_geom <- data$p_geom[1]
        d <- rep(-Inf, length(m))

        idx <- m > -1
        if (any(idx)) {
            d[idx] <- f_density(m[idx], offset, p_geom)
        }

        if (!log) {
            d[idx] <- exp(d[idx])
            d[!idx] <- 0.0
        }

        d
    })

    unsplit(d, data_f)
}

#' @rdname vmgeom
#' @export
pvmgeom <- function(m, lm, r, lower.tail = TRUE, log = FALSE, perception_fun = vmperception) {
    if (anyNA(m) || anyNA(lm) || anyNA(r)) {
        stop("NA's are not allowed!")
    }

    if (any(is.infinite(lm))) {
        stop("Infinite limiting magnitudes are not allowed!")
    }

    if (any(is.infinite(r))) {
        stop("Infinite r values are not allowed!")
    }

    if (any(1.0 > r)) {
        stop(paste0('r must be greater than 1.0 instead of "', r, '"!'))
    }

    is_wholenumber <- function(x, tol = .Machine$double.eps^0.5) all(is.infinite(x) | abs(x - round(x)) < tol)
    if (!is_wholenumber(m)) {
        stop("magnitudes must be integer values!")
    }

    perception_fun <- match.fun(perception_fun)

    p_geom <- 1.0 - 1.0 / r

    if (1 == length(r)) {
        r <- rep(r, length(m))
        p_geom <- rep(p_geom, length(m))
    }

    std_res <- .vmgeom_std(m, lm)
    m <- std_res$m
    offset <- std_res$offset

    f_density <- function(m, offset, p_geom) {
        stats::dgeom(m, p_geom) * perception_fun(m + offset)
    }

    f_sum <- Vectorize(\(m, offset, p_geom) {
        m <- as.integer(seq(0, m))
        sum(f_density(m, offset, p_geom))
    })

    f_prob <- function(m, offset, p_geom) {
        m_max <- 15L
        norm <- .vmgeom_norm(offset, p_geom, m_max, perception_fun)
        p <- rep(0.0, length(m))

        if (lower.tail) {
            idx <- m <= m_max
            if (any(idx)) {
                p[idx] <- 1.0 - f_sum(m[idx], offset, p_geom) / norm
            }

            idx <- m > m_max
            if (any(idx)) {
                p[idx] <- stats::pgeom(m[idx], p_geom, lower.tail = FALSE) / norm
            }
        } else {
            idx <- m <= m_max
            if (any(idx)) {
                p[idx] <- f_sum(m[idx], offset, p_geom) / norm
            }

            idx <- m > m_max
            if (any(idx)) {
                p[idx] <- 1.0 - stats::pgeom(m[idx], p_geom, lower.tail = FALSE) / norm
            }
        }

        p
    }

    arg_data <- data.frame(
        m = m,
        offset = offset,
        p_geom = p_geom
    )

    data_f <- as.factor(paste0(offset, "/", p_geom))
    data_s <- split(arg_data, data_f)
    p <- lapply(data_s, \(data) {
        m <- data$m
        offset <- data$offset[1]
        p_geom <- data$p_geom[1]

        if (lower.tail) {
            p <- rep(1.0, length(m))
        } else {
            p <- rep(0.0, length(m))
        }

        idx <- m > -1
        if (any(idx)) {
            p[idx] <- f_prob(m[idx], offset, p_geom)
        }

        if (log) {
            p[idx] <- base::log(p[idx])
            if (lower.tail) {
                p[!idx] <- 0.0
            } else {
                p[!idx] <- -Inf
            }
        }

        p
    })

    unsplit(p, data_f)
}

#' @rdname vmgeom
#' @export
qvmgeom <- function(p, lm, r, lower.tail = TRUE, perception_fun = vmperception) {
    if (anyNA(p) || anyNA(lm) || anyNA(r)) {
        stop("NA's are not allowed!")
    }

    if (any(is.infinite(lm))) {
        stop("Infinite limiting magnitudes are not allowed!")
    }

    if (any(is.infinite(r))) {
        stop("Infinite r values are not allowed!")
    }

    if (any(1.0 > r)) {
        stop(paste0('r must be greater than 1.0 instead of "', r, '"!'))
    }

    perception_fun <- match.fun(perception_fun)

    p_geom <- 1.0 - 1.0 / r

    if (1 == length(r)) {
        r <- rep(r, length(p))
        p_geom <- rep(p_geom, length(p))
    }

    lm_round <- round(lm)
    offset <- lm - lm_round
    lm <- lm_round

    idx <- -0.5 == offset
    lm[idx] <- lm[idx] - 1L
    offset[idx] <- offset[idx] + 1.0

    if (1 == length(offset)) {
        offset <- rep(offset, length(p))
    }

    arg_data <- data.frame(
        p = p,
        offset = offset,
        p_geom = p_geom,
        r = r
    )

    data_f <- as.factor(paste0(offset, "/", p_geom))
    data_s <- split(arg_data, data_f)
    m <- lapply(data_s, \(data) {
        m_max <- 15L
        p <- data$p
        offset <- data$offset[1]
        p_geom <- data$p_geom[1]
        r <- data$r[1]
        m <- rep(NA, length(p))

        if (lower.tail) {
            m[1.0 == p] <- 0L
            p_max <- 1.0 - vismeteor::pvmgeom(0, m_max + offset, r, lower.tail = FALSE, perception_fun = perception_fun)
            idx <- p >= 0.0 & p < p_max
            if (any(idx)) {
                m[idx] <- m_max + 1 + stats::qgeom(p[idx] / p_max, p_geom, lower.tail = FALSE)
            }

            idx <- p >= p_max & p < 1.0
            if (any(idx)) {
                m0 <- seq(-m_max, 0, 1)
                p0 <- c(vismeteor::pvmgeom(m0, offset, r, lower.tail = TRUE, perception_fun = perception_fun), 1.0)
                m0 <- c(m0, 1L)
                p_idx <- findInterval(p[idx], p0, left.open = FALSE)
                m[idx] <- -m0[p_idx]
            }
        } else {
            m[0.0 == p] <- 0L
            p_max <- vismeteor::pvmgeom(0, m_max + offset, r, lower.tail = FALSE, perception_fun = perception_fun)
            idx <- p > p_max & p <= 1.0
            if (any(idx)) {
                m[idx] <- m_max + 1L + stats::qgeom((p[idx] - p_max) / (1.0 - p_max), p_geom)
            }

            idx <- p > 0.0 & p <= p_max
            if (any(idx)) {
                m0 <- seq(1, -m_max, -1)
                p0 <- vismeteor::pvmgeom(m0, offset, r, lower.tail = FALSE, perception_fun = perception_fun)
                p_idx <- findInterval(p[idx], p0, left.open = TRUE) + 1
                m[idx] <- -m0[p_idx]
                m[m < 0] <- NA
            }
        }

        m
    })

    m <- lm - unsplit(m, data_f)

    if (anyNA(m)) {
        warning("NaNs produced")
    }

    as.numeric(m)
}

#' @rdname vmgeom
#' @export
rvmgeom <- function(n, lm, r, perception_fun = vmperception) {
    if (anyNA(lm) || anyNA(r)) {
        stop("NA's are not allowed!")
    }

    if (any(is.infinite(lm))) {
        stop("Infinite limiting magnitudes are not allowed!")
    }

    if (any(is.infinite(r))) {
        stop("Infinite r values are not allowed!")
    }

    if (any(1.0 > r)) {
        stop(paste0('r must be greater than 1.0 instead of "', r, '"!'))
    }

    perception_fun <- match.fun(perception_fun)

    p <- stats::runif(n)
    m <- rep(NA, n)
    idx <- p < 0.5

    if (any(idx)) {
        m[idx] <- vismeteor::qvmgeom(p[idx], lm, r, lower.tail = TRUE, perception_fun = perception_fun)
    }

    if (any(!idx)) {
        m[!idx] <- vismeteor::qvmgeom(1.0 - p[!idx], lm, r, lower.tail = FALSE, perception_fun = perception_fun)
    }

    m
}

#' standardization
#'
#' @noRd
.vmgeom_std <- function(m, lm) {
    m <- lm - m
    m_round <- round(m)
    offset <- rep(0.0, length(m))
    idx <- is.infinite(m)
    offset[!idx] <- m[!idx] - m_round[!idx]
    m <- m_round

    idx <- -0.5 == offset
    m[idx] <- m[idx] - 1L
    offset[idx] <- offset[idx] + 1.0

    list(
        m = m,
        offset = offset
    )
}

#' normalization
#'
#' @noRd
.vmgeom_norm <- function(offset, p_geom, m_max, perception_fun) {
    m <- as.integer(seq(0, m_max))
    sum(stats::dgeom(m, p_geom) * perception_fun(m + offset)) +
        stats::pgeom(m_max, p_geom, lower.tail = FALSE)
}
