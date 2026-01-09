#' @name vmideal
#' @aliases dvmideal
#' @aliases pvmideal
#' @aliases qvmideal
#' @aliases rvmideal
#' @aliases cvmideal
#' @title Ideal Distribution of Visual Meteor Magnitudes
#' @description
#' Density, distribution function, quantile function, and random generation
#' for the ideal distribution of visual meteor magnitudes.
#' @param psi numeric; the location parameter of the probability distribution.
#' @param m integer; visual meteor magnitude.
#' @param lm numeric; limiting magnitude.
#' @param p numeric; probability.
#' @param n numeric; count of meteor magnitudes.
#' @param log logical; if `TRUE`, probabilities are returned as `log(p)`.
#' @param lower.tail logical; if `TRUE` (default), probabilities are
#'     \eqn{P[M < m]}; otherwise, \eqn{P[M \ge m]}.
#' @param perception.fun function; optional perception probability function.
#'     The default is [vismeteor::vmperception].
#'
#' @details
#' The density of the [ideal distribution of meteor magnitudes][vismeteor::mideal] is
#' \deqn{
#'     {\displaystyle f(m) = \frac{\mathrm{d}p}{\mathrm{d}m} = \frac{3}{2} \, \log(r) \sqrt{\frac{r^{3 \, \psi + 2 \, m}}{(r^\psi + r^m)^5}}}
#' }
#' where \eqn{m} denotes the continuous (real-valued) meteor magnitude,
#' \eqn{r = 10^{0.4} \approx 2.51189 \dots} is a constant, and
#' \eqn{\psi} is the only parameter of this magnitude distribution.
#'
#' In visual meteor observations, magnitudes are usually estimated as integer values.
#' Hence, this distribution is discrete and its probability mass function is given by
#' \deqn{
#' P[M = m] \sim
#' \begin{cases}
#'   g(m_{\mathrm{lim}} - m) \displaystyle \int\limits_{m-0.5}^{m+0.5} f(u) \, \mathrm{d}u, & \text{if } m_{\mathrm{lim}} - m > -0.5,\\[5pt]
#'   0 & \text{otherwise,}
#' \end{cases}
#' }
#' where \eqn{m_{\mathrm{lim}}} denotes the limiting (non-integer) magnitude of the observation,
#' and \eqn{m} the integer meteor magnitude.
#' The function \eqn{f(\cdot)} is the continuous density of the ideal magnitude distribution,
#' and \eqn{g(\cdot)} denotes the [perception probability function][vismeteor::vmperception].
#'
#' If a different perception probability function `perception.fun` is supplied,
#' it must have the signature `function(x)` and return the perception probabilities of
#' the difference `x` between the limiting magnitude and the meteor magnitude.
#' If `x >= 15.0`, the `perception.fun` function should return a perception probability of `1.0`.
#' The argument `perception.fun` is resolved using [match.fun].
#'
#' @return
#' - `dvmideal`: density
#' - `pvmideal`: distribution function
#' - `qvmideal`: quantile function
#' - `rvmideal`: random generation
#' - `cvmideal`: partial convolution of the ideal distribution of meteor magnitudes
#'   with the perception probabilities.
#'
#' The length of the result is determined by `n` for `rvmideal`, and is the maximum
#' of the lengths of the numeric vector arguments for the other functions.
#' All arguments are vectorized; standard R recycling rules apply.
#'
#' Since the distribution is discrete, `qvmideal` and `rvmideal` always return integer values.
#' `qvmideal` may return `NaN` with a warning.
#' @seealso [vismeteor::mideal]
#'   [vismeteor::vmperception]
#'
#' @references Richter, J. (2018) \emph{About the mass and magnitude distributions of meteor showers}.
#'   WGN, Journal of the International Meteor Organization, vol. 46, no. 1, p. 34-38
#' @examples
#' N <- 100
#' psi <- 5.0
#' limmag <- 6.5
#' (m <- seq(6, -4))
#'
#' # discrete density of `N` meteor magnitudes
#' (freq <- round(N * dvmideal(m, limmag, psi)))
#'
#' # log likelihood function
#' lld <- function(psi) {
#'     -sum(freq * dvmideal(m, limmag, psi, log=TRUE))
#' }
#'
#' # maximum likelihood estimation (MLE) of psi
#' est <- optim(2, lld, method='Brent', lower=0, upper=8, hessian=TRUE)
#'
#' # estimations
#' est$par # mean of psi
#'
#' # generate random meteor magnitudes
#' m <- rvmideal(N, limmag, psi)
#'
#' # log likelihood function
#' llr <- function(psi) {
#'     -sum(dvmideal(m, limmag, psi, log=TRUE))
#' }
#'
#' # maximum likelihood estimation (MLE) of psi
#' est <- optim(2, llr, method='Brent', lower=0, upper=8, hessian=TRUE)
#'
#' # estimations
#' est$par # mean of psi
#' sqrt(1/est$hessian[1][1]) # standard deviation of psi
#'
#' m <- seq(6, -4, -1)
#' p <- vismeteor::dvmideal(m, limmag, psi)
#' barplot(
#'     p,
#'     names.arg = m,
#'     main = paste0('Density (psi = ', psi, ', limmag = ', limmag, ')'),
#'     col = "blue",
#'     xlab = 'm',
#'     ylab = 'p',
#'     border = "blue",
#'     space = 0.5
#' )
#' axis(side = 2, at = pretty(p))
#'
#' plot(
#'     function(lm) vismeteor::cvmideal(lm, psi, log = TRUE),
#'     -5, 10,
#'     main = paste0(
#'         'Partial convolution of the ideal meteor magnitude distribution\n',
#'         'with the perception probabilities (psi = ', psi, ')'
#'     ),
#'     col = "blue",
#'     xlab = 'lm',
#'     ylab = 'log(rate)'
#' )


#' @rdname vmideal
#' @export
dvmideal <- function(m, lm, psi, log = FALSE, perception.fun = NULL) {
    if (anyNA(m) | anyNA(lm) | anyNA(psi)) {
        stop("NA's are not allowed!")
    }

    if (any(is.infinite(lm))) {
        stop("Infinite limiting magnitudes are not allowed!")
    }

    if (any(is.infinite(psi))) {
        stop("Infinite psi values are not allowed!")
    }

    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) all(is.infinite(x) | abs(x - round(x)) < tol)
    if (! is.wholenumber(m)) {
        stop("magnitudes must be integer values!")
    }

    if (is.null(perception.fun)) {
        perception.fun <- vismeteor::vmperception
    } else {
        perception.fun <- match.fun(perception.fun)
    }

    if (1 == length(lm)) {
        lm <- rep(lm, length(m))
    }

    if (1 == length(psi)) {
        psi <- rep(psi, length(m))
    }

    # density function
    f.density <- function(m, lm, psi, log) {
        psi.exp <- 10.0
        if (lm + psi.exp < psi) {
            return(vismeteor::dvmgeom(m, lm, 10^0.4, log = log))
        }

        norm.res <- vmideal.norm(lm, psi, perception.fun)
        p <- norm.res$p
        m.lower <- norm.res$m.lower
        m.upper <- norm.res$m.upper

        if (log) {
            d <- rep(-Inf, length(m))
        } else {
            d <- rep(0.0, length(m))
        }

        idx <- m >= m.lower & m <= m.upper
        if (any(idx)) {
            d.tmp <- p[as.character(m[idx])]

            if (log) {
                d.tmp[0.0 == d.tmp] <- -Inf

                log.idx <- -Inf != d.tmp
                if (any(log.idx)) {
                    d.tmp[log.idx] <- base::log(d.tmp[log.idx])
                }
            }

            d[idx] <- d.tmp

        }

        idx <- m > -Inf & m < m.lower
        if (any(idx)) {
            if (log) {
                d[idx] <- dmideal.int(m[idx], psi, log = TRUE) - base::log(norm.res$norm)
            } else {
                d[idx] <- dmideal.int(m[idx], psi)/norm.res$norm
            }
        }

        d
    }

    arg.data <- data.frame(
        m = m,
        lm = lm,
        psi = psi
    )

    data.f <- as.factor(paste0(lm, '/', psi))
    data.s <- split(arg.data, data.f)
    d <- lapply(data.s, function(data) {
        m <- data$m
        lm <- data$lm[1]
        psi <- data$psi[1]
        f.density(m, lm, psi, log = log)
    })

    unsplit(d, data.f)
}

#' @rdname vmideal
#' @export
pvmideal <- function(m, lm, psi, lower.tail = TRUE, log = FALSE, perception.fun = NULL) {
    if (anyNA(m) | anyNA(lm) | anyNA(psi)) {
        stop("NA's are not allowed!")
    }

    if (any(is.infinite(lm))) {
        stop("Infinite limiting magnitudes are not allowed!")
    }

    if (any(is.infinite(psi))) {
        stop("Infinite psi values are not allowed!")
    }

    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) all(is.infinite(x) | abs(x - round(x)) < tol)
    if (! is.wholenumber(m)) {
        stop("magnitudes must be integer values!")
    }

    if (is.null(perception.fun)) {
        perception.fun <- vismeteor::vmperception
    } else {
        perception.fun <- match.fun(perception.fun)
    }

    if (1 == length(lm)) {
        lm <- rep(lm, length(m))
    }

    if (1 == length(psi)) {
        psi <- rep(psi, length(m))
    }

    # probability function
    f.prob <- function(m, lm, psi, log) {
        psi.exp <- 10.0
        if (lm + psi.exp < psi) {
            return(vismeteor::pvmgeom(m, lm, 10^0.4, lower.tail = lower.tail, log = log))
        }

        norm.res <- vmideal.norm(lm, psi, perception.fun)
        m.lower <- norm.res$m.lower
        m.upper <- norm.res$m.upper

        if (lower.tail) {
            if (log) {
                p <- rep(0.0, length(m))
                p[-Inf == m] <- -Inf
            } else {
                p <- rep(1.0, length(m))
                p[-Inf == m] <- 0.0
            }

            idx <- m > -Inf & m <= m.lower
            if (any(idx)) {
                if (log) {
                    p[idx] <- vismeteor::pmideal(m[idx] - 0.5, psi, lower.tail = TRUE, log = TRUE) -
                        base::log(norm.res$norm)
                } else {
                    p[idx] <- vismeteor::pmideal(m[idx] - 0.5, psi, lower.tail = TRUE, log = FALSE) / norm.res$norm
                }
            }

            idx <- m > m.lower & m <= m.upper
            if (any(idx)) {
                p.gen <- cumsum(norm.res$p)
                p[idx] <- norm.res$p.lower.tail + p.gen[as.character(m[idx] - 1)]
                p[p>1.0] <- 1.0
                if (log) {
                    log.idx <- p > 0.0
                    if (any(log.idx)) {
                        p[log.idx] <- base::log(p[log.idx])
                    }
                }
            }
        } else {
            p <- rep(0.0, length(m))
            p[-Inf == m] <- 1.0

            idx <- m > -Inf & m < m.lower
            if (any(idx)) {
                p[idx] <- 1.0 - vismeteor::pmideal(m[idx] - 0.5, psi, lower.tail = TRUE) / norm.res$norm
            }

            idx <- m >= m.lower & m <= m.upper
            if (any(idx)) {
                p.gen <- base::rev(cumsum(base::rev(norm.res$p)))
                p[idx] <- p.gen[as.character(m[idx])]
            }

            p[m > 0.0 & is.infinite(m)] <- 0.0
            p[p>1.0] <- 1.0

            if (log) {
                log.idx <- p > 0.0
                if (any(log.idx)) {
                    p[log.idx] <- base::log(p[log.idx])
                }
                p[!log.idx] <- -Inf
            }
        }

        p
    }

    arg.data <- data.frame(
        m = m,
        lm = lm,
        psi = psi
    )

    data.f <- as.factor(paste0(lm, '/', psi))
    data.s <- split(arg.data, data.f)
    p <- lapply(data.s, function(data) {
        m <- data$m
        lm <- data$lm[1]
        psi <- data$psi[1]
        f.prob(m, lm, psi, log = log)
    })

    unsplit(p, data.f)
}

#' @rdname vmideal
#' @export
qvmideal <- function(p, lm, psi, lower.tail = TRUE, perception.fun = NULL) {
    if (anyNA(p) | anyNA(lm) | anyNA(psi)) {
        stop("NA's are not allowed!")
    }

    if (any(is.infinite(lm))) {
        stop("Infinite limiting magnitudes are not allowed!")
    }

    if (any(is.infinite(psi))) {
        stop("Infinite psi values are not allowed!")
    }

    if (is.null(perception.fun)) {
        perception.fun <- vismeteor::vmperception
    } else {
        perception.fun <- match.fun(perception.fun)
    }

    if (1 == length(lm)) {
        lm <- rep(lm, length(p))
    }

    if (1 == length(psi)) {
        psi <- rep(psi, length(p))
    }

    # quantile function
    f.q <- function(p, lm, psi) {
        m.max <- 15L
        psi.exp <- 10.0

        if (lm + psi.exp < psi) {
            r.lower <- 10^0.4
            return(vismeteor::qvmgeom(p, lm, r.lower, lower.tail = lower.tail))
        }

        m.upper <- vmideal.upper.lm(lm)
        m.lower <- m.upper - m.max
        m <- rep(NA, length(p))

        if(lower.tail) {
            m[0.0 == p] <- -Inf
            m[1.0 == p] <- m.upper
            p.max <- vismeteor::pvmideal(m.lower, lm, psi, lower.tail = TRUE, perception.fun = perception.fun)

            idx <- p > 0.0 & p < p.max
            if (any(idx)) {
                p.max.ideal <- vismeteor::pmideal(m.lower - 0.5, psi, lower.tail = TRUE)
                m[idx] <- floor(0.5 + vismeteor::qmideal(
                    p.max.ideal * p[idx] / p.max,
                    psi,
                    lower.tail = TRUE
                ))
            }
            idx <- p>=p.max & p<1.0
            if (any(idx)) {
                m.gen <- seq(m.lower, m.upper + 1)
                p.gen <- vismeteor::pvmideal(m.gen, lm, psi, lower.tail = lower.tail, perception.fun = perception.fun)
                p.idx <- findInterval(p[idx], p.gen, left.open = FALSE)
                m[idx] <- m.gen[p.idx]
            }
        } else {
            m[0.0 == p] <- m.upper
            m[1.0 == p] <- -Inf
            p.max <- vismeteor::pvmideal(m.lower, lm, psi, lower.tail = FALSE, perception.fun = perception.fun)

            idx <- p > p.max & p < 1.0
            if (any(idx)) {
                p.max.ideal <- vismeteor::pmideal(m.lower - 0.5, psi, lower.tail = FALSE)
                m[idx] <- floor(0.5 + vismeteor::qmideal(
                    1.0 - (1.0 - p[idx]) * ((1.0 - p.max.ideal) / (1.0 - p.max)),
                    psi,
                    lower.tail = FALSE
                ))
            }

            idx <- p>0.0 & p<=p.max
            if (any(idx)) {
                m.gen <- seq(m.upper + 1, m.lower)
                p.gen <- vismeteor::pvmideal(m.gen, lm, psi, lower.tail = FALSE, perception.fun = perception.fun)
                p.idx <- findInterval(p[idx], p.gen, left.open = TRUE) + 1
                m[idx] <- m.gen[p.idx]
            }
        }

        m
    }

    arg.data <- data.frame(
        p = p,
        lm = lm,
        psi = psi
    )

    data.f <- as.factor(paste0(lm, '/', psi))
    data.s <- split(arg.data, data.f)
    m <- lapply(data.s, function(data) {
        p <- data$p
        lm <- data$lm[1]
        psi <- data$psi[1]
        f.q(p, lm, psi)
    })

    m <- unsplit(m, data.f)

    if (anyNA(m)) {
        warning('NaNs produced')
    }

    m
}

#' @rdname vmideal
#' @export
rvmideal <- function(n, lm, psi, perception.fun = NULL) {
    if (anyNA(lm) | anyNA(psi)) {
        stop("NA's are not allowed!")
    }

    if (any(is.infinite(lm))) {
        stop("Infinite limiting magnitudes are not allowed!")
    }

    if (any(is.infinite(psi))) {
        stop("Infinite psi values are not allowed!")
    }

    if (is.null(perception.fun)) {
        perception.fun <- vismeteor::vmperception
    } else {
        perception.fun <- match.fun(perception.fun)
    }

    p <- stats::runif(n)
    m <- rep(NA, n)

    idx <- p < 0.5
    if (any(idx)) {
        m[idx] <- vismeteor::qvmideal(p[idx], lm, psi, lower.tail = TRUE, perception.fun = perception.fun)
    }

    if (any(!idx)) {
        m[!idx] <- vismeteor::qvmideal(1.0 - p[!idx], lm, psi, lower.tail = FALSE, perception.fun = perception.fun)
    }

    m
}

#' @rdname vmideal
#' @export
cvmideal <- function(lm, psi, log = FALSE, perception.fun = NULL) {
    if (anyNA(lm) | anyNA(psi)) {
        stop("NA's are not allowed!")
    }

    if (is.null(perception.fun)) {
        perception.fun <- vismeteor::vmperception
    } else {
        perception.fun <- match.fun(perception.fun)
    }

    if (1 == length(psi)) {
        psi <- rep(psi, length(lm))
    }

    # Integration - similar to vmideal.norm()
    f.integrate <- function(lm, psi) {
        m.max <- 15L
        m.upper <- vmideal.upper.lm(lm)
        m.lower <- m.upper - m.max
        m <- as.integer(seq(m.lower, m.upper))
        p <- dmideal.int(m, psi) * perception.fun(lm - m)
        p.lower.tail <- vismeteor::pmideal(m.lower - 0.5, psi, lower.tail = TRUE)
        sum(p) + p.lower.tail
    }

    p <- mapply(function(lm, psi){
        if (Inf == lm & Inf == psi) return(NA)
        if (-Inf == lm & -Inf == psi) return(NA)
        if (Inf == lm & -Inf == psi) return(if(log) 0.0 else 1.0)
        if (Inf == psi & -Inf == lm) return(if(log) -Inf else 0.0)

        if (log) {
            base::log(f.integrate(lm, psi))
        } else {
            f.integrate(lm, psi)
        }
    }, lm, psi, SIMPLIFY = TRUE)

    if (anyNA(p)) {
        warning('NaNs produced')
    }

    p
}

#' upper available magnitude
#'
#' @noRd
vmideal.upper.lm <- function(lm) {
    lm.round <- round(lm)
    offset <- lm - lm.round
    if (-0.5 == offset) {
        lm.round <- lm.round - 1L
    }

    lm.round
}

#' density of ideal integer magnitude distribution
#'
#' @noRd
dmideal.int <- function(m, psi, log = FALSE) {
    psi.exp <- 10.0
    r.lower <- 10^0.4

    a <- -base::log(r.lower)
    d <- rep(NA, length(m))

    idx <- m > (psi + psi.exp)
    if (any(idx)) {
        if (log) {
            d[idx] <- base::log(1 - base::exp(1.5 * a)) + a * (1.5 * (m[idx] - psi) - 0.75)
        } else {
            d[idx] <- base::exp(a * (1.5 * (m[idx] - psi) - 0.75)) -
                base::exp(a * (1.5 * (m[idx] - psi) + 0.75))
        }
    }

    idx <- m < (psi - psi.exp)
    if (any(idx)) {
        if (log) {
            d[idx] <- base::log(1.5) + a * (psi - m[idx] - 0.5) + base::log(1 - base::exp(a))
        } else {
            d[idx] <- 1.5 * base::exp(a * (psi - m[idx] - 0.5)) -
                1.5 * base::exp(a * (psi - m[idx] + 0.5))
        }
    }

    idx <- is.na(d)
    if (any(idx)) {
        d[idx] <- vismeteor::pmideal(m[idx] + 0.5, psi, lower.tail = TRUE) -
            vismeteor::pmideal(m[idx] - 0.5, psi, lower.tail = TRUE)

        if (log) {
            d[idx] <- base::log(d[idx])
        }
    }

    d
}

#' normalization
#'
#' @noRd
vmideal.norm <- function(lm, psi, perception.fun) {
    m.max <- 15L
    m.upper <- vmideal.upper.lm(lm)
    m.lower <- m.upper - m.max
    m <- as.integer(seq(m.lower, m.upper))
    p <- dmideal.int(m, psi) * perception.fun(lm - m)
    names(p) <- as.character(m)
    p.lower.tail <- vismeteor::pmideal(m.lower - 0.5, psi, lower.tail = TRUE)
    norm <- sum(p) + p.lower.tail

    list(
        m.lower = m.lower,
        m.upper = m.upper,
        norm = norm,
        p = p/norm,
        p.lower.tail = p.lower.tail/norm
    )
}
