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
#' @param perception_fun function; optional perception probability function.
#'     The default is [vismeteor::vmperception].
#'
#' @details
#' The density of the [ideal distribution of meteor magnitudes][vismeteor::mideal] is
#' \deqn{
#'     {\displaystyle f(m) = \frac{\mathrm{d}p}{\mathrm{d}m} =
#'      \frac{3}{2} \, \log(r) \sqrt{\frac{r^{3 \, \psi + 2 \, m}}{(r^\psi + r^m)^5}}}
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
#'   g(m_{\mathrm{lim}} - m) \displaystyle \int\limits_{m-0.5}^{m+0.5} f(u) \, \mathrm{d}u,
#'   & \text{if } m_{\mathrm{lim}} - m > -0.5,\\[5pt]
#'   0 & \text{otherwise,}
#' \end{cases}
#' }
#' where \eqn{m_{\mathrm{lim}}} denotes the limiting (non-integer) magnitude of the observation,
#' and \eqn{m} the integer meteor magnitude.
#' The function \eqn{f(\cdot)} is the continuous density of the ideal magnitude distribution,
#' and \eqn{g(\cdot)} denotes the [perception probability function][vismeteor::vmperception].
#'
#' If a different perception probability function `perception_fun` is supplied,
#' it must have the signature `function(x)` and return the perception probabilities of
#' the difference `x` between the limiting magnitude and the meteor magnitude.
#' If `x >= 15.0`, the `perception_fun` function should return a perception probability of `1.0`.
#' The argument `perception_fun` is resolved using [match.fun].
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
#'     -sum(freq * dvmideal(m, limmag, psi, log = TRUE))
#' }
#'
#' # maximum likelihood estimation (MLE) of psi
#' est <- optim(2, lld, method = "Brent", lower = 0, upper = 8, hessian = TRUE)
#'
#' # estimations
#' est$par # mean of psi
#'
#' # generate random meteor magnitudes
#' m <- rvmideal(N, limmag, psi)
#'
#' # log likelihood function
#' llr <- function(psi) {
#'     -sum(dvmideal(m, limmag, psi, log = TRUE))
#' }
#'
#' # maximum likelihood estimation (MLE) of psi
#' est <- optim(2, llr, method = "Brent", lower = 0, upper = 8, hessian = TRUE)
#'
#' # estimations
#' est$par # mean of psi
#' sqrt(1 / est$hessian[1][1]) # standard deviation of psi
#'
#' m <- seq(6, -4, -1)
#' p <- vismeteor::dvmideal(m, limmag, psi)
#' barplot(
#'     p,
#'     names.arg = m,
#'     main = paste0("Density (psi = ", psi, ", limmag = ", limmag, ")"),
#'     col = "blue",
#'     xlab = "m",
#'     ylab = "p",
#'     border = "blue",
#'     space = 0.5
#' )
#' axis(side = 2, at = pretty(p))
#'
#' plot(
#'     \(lm) vismeteor::cvmideal(lm, psi, log = TRUE),
#'     -5, 10,
#'     main = paste0(
#'         "Partial convolution of the ideal meteor magnitude distribution\n",
#'         "with the perception probabilities (psi = ", psi, ")"
#'     ),
#'     col = "blue",
#'     xlab = "lm",
#'     ylab = "log(rate)"
#' )

#' @rdname vmideal
#' @export
dvmideal <- function(m, lm, psi, log = FALSE, perception_fun = vmperception) {
    if (anyNA(m) || anyNA(lm) || anyNA(psi)) {
        stop("NA's are not allowed!")
    }

    if (any(is.infinite(lm))) {
        stop("Infinite limiting magnitudes are not allowed!")
    }

    if (any(is.infinite(psi))) {
        stop("Infinite psi values are not allowed!")
    }

    is_wholenumber <- function(x, tol = .Machine$double.eps^0.5) all(is.infinite(x) | abs(x - round(x)) < tol)
    if (!is_wholenumber(m)) {
        stop("magnitudes must be integer values!")
    }

    perception_fun <- match.fun(perception_fun)

    if (1 == length(lm)) {
        lm <- rep(lm, length(m))
    }

    if (1 == length(psi)) {
        psi <- rep(psi, length(m))
    }

    # density function
    f_density <- function(m, lm, psi, log) {
        psi_exp <- 10.0
        if (lm + psi_exp < psi) {
            return(vismeteor::dvmgeom(m, lm, 10^0.4, log = log))
        }

        norm_res <- .vmideal_norm(lm, psi, perception_fun)
        p <- norm_res$p
        m_lower <- norm_res$m_lower
        m_upper <- norm_res$m_upper

        if (log) {
            d <- rep(-Inf, length(m))
        } else {
            d <- rep(0.0, length(m))
        }

        idx <- m >= m_lower & m <= m_upper
        if (any(idx)) {
            d_tmp <- p[as.character(m[idx])]

            if (log) {
                d_tmp[0.0 == d_tmp] <- -Inf

                log_idx <- -Inf != d_tmp
                if (any(log_idx)) {
                    d_tmp[log_idx] <- base::log(d_tmp[log_idx])
                }
            }

            d[idx] <- d_tmp
        }

        idx <- m > -Inf & m < m_lower
        if (any(idx)) {
            if (log) {
                d[idx] <- .dmideal_int(m[idx], psi, log = TRUE) - base::log(norm_res$norm)
            } else {
                d[idx] <- .dmideal_int(m[idx], psi) / norm_res$norm
            }
        }

        d
    }

    arg_data <- data.frame(
        m = m,
        lm = lm,
        psi = psi
    )

    data_f <- as.factor(paste0(lm, "/", psi))
    data_s <- split(arg_data, data_f)
    d <- lapply(data_s, \(data) {
        m <- data$m
        lm <- data$lm[1]
        psi <- data$psi[1]
        f_density(m, lm, psi, log = log)
    })

    unsplit(d, data_f)
}

#' @rdname vmideal
#' @export
pvmideal <- function(m, lm, psi, lower.tail = TRUE, log = FALSE, perception_fun = vmperception) {
    if (anyNA(m) || anyNA(lm) || anyNA(psi)) {
        stop("NA's are not allowed!")
    }

    if (any(is.infinite(lm))) {
        stop("Infinite limiting magnitudes are not allowed!")
    }

    if (any(is.infinite(psi))) {
        stop("Infinite psi values are not allowed!")
    }

    is_wholenumber <- function(x, tol = .Machine$double.eps^0.5) all(is.infinite(x) | abs(x - round(x)) < tol)
    if (!is_wholenumber(m)) {
        stop("magnitudes must be integer values!")
    }

    perception_fun <- match.fun(perception_fun)

    if (1 == length(lm)) {
        lm <- rep(lm, length(m))
    }

    if (1 == length(psi)) {
        psi <- rep(psi, length(m))
    }

    # probability function
    f_prob <- function(m, lm, psi, log) {
        psi_exp <- 10.0
        if (lm + psi_exp < psi) {
            return(vismeteor::pvmgeom(m, lm, 10^0.4, lower.tail = lower.tail, log = log))
        }

        norm_res <- .vmideal_norm(lm, psi, perception_fun)
        m_lower <- norm_res$m_lower
        m_upper <- norm_res$m_upper

        if (lower.tail) {
            if (log) {
                p <- rep(0.0, length(m))
                p[-Inf == m] <- -Inf
            } else {
                p <- rep(1.0, length(m))
                p[-Inf == m] <- 0.0
            }

            idx <- m > -Inf & m <= m_lower
            if (any(idx)) {
                if (log) {
                    p[idx] <- vismeteor::pmideal(m[idx] - 0.5, psi, lower.tail = TRUE, log = TRUE) -
                        base::log(norm_res$norm)
                } else {
                    p[idx] <- vismeteor::pmideal(m[idx] - 0.5, psi, lower.tail = TRUE, log = FALSE) / norm_res$norm
                }
            }

            idx <- m > m_lower & m <= m_upper
            if (any(idx)) {
                p_gen <- cumsum(norm_res$p)
                p[idx] <- norm_res$p_lower_tail + p_gen[as.character(m[idx] - 1)]
                p[p > 1.0] <- 1.0
                if (log) {
                    log_idx <- p > 0.0
                    if (any(log_idx)) {
                        p[log_idx] <- base::log(p[log_idx])
                    }
                }
            }
        } else {
            p <- rep(0.0, length(m))
            p[-Inf == m] <- 1.0

            idx <- m > -Inf & m < m_lower
            if (any(idx)) {
                p[idx] <- 1.0 - vismeteor::pmideal(m[idx] - 0.5, psi, lower.tail = TRUE) / norm_res$norm
            }

            idx <- m >= m_lower & m <= m_upper
            if (any(idx)) {
                p_gen <- base::rev(cumsum(base::rev(norm_res$p)))
                p[idx] <- p_gen[as.character(m[idx])]
            }

            p[m > 0.0 & is.infinite(m)] <- 0.0
            p[p > 1.0] <- 1.0

            if (log) {
                log_idx <- p > 0.0
                if (any(log_idx)) {
                    p[log_idx] <- base::log(p[log_idx])
                }
                p[!log_idx] <- -Inf
            }
        }

        p
    }

    arg_data <- data.frame(
        m = m,
        lm = lm,
        psi = psi
    )

    data_f <- as.factor(paste0(lm, "/", psi))
    data_s <- split(arg_data, data_f)
    p <- lapply(data_s, \(data) {
        m <- data$m
        lm <- data$lm[1]
        psi <- data$psi[1]
        f_prob(m, lm, psi, log = log)
    })

    unsplit(p, data_f)
}

#' @rdname vmideal
#' @export
qvmideal <- function(p, lm, psi, lower.tail = TRUE, perception_fun = vmperception) {
    if (anyNA(p) || anyNA(lm) || anyNA(psi)) {
        stop("NA's are not allowed!")
    }

    if (any(is.infinite(lm))) {
        stop("Infinite limiting magnitudes are not allowed!")
    }

    if (any(is.infinite(psi))) {
        stop("Infinite psi values are not allowed!")
    }

    perception_fun <- match.fun(perception_fun)

    if (1 == length(lm)) {
        lm <- rep(lm, length(p))
    }

    if (1 == length(psi)) {
        psi <- rep(psi, length(p))
    }

    # quantile function
    f_q <- function(p, lm, psi) {
        m_max <- 15L
        psi_exp <- 10.0

        if (lm + psi_exp < psi) {
            r_lower <- 10^0.4
            return(vismeteor::qvmgeom(p, lm, r_lower, lower.tail = lower.tail))
        }

        m_upper <- .vmideal_upper_lm(lm)
        m_lower <- m_upper - m_max
        m <- rep(NA, length(p))

        if (lower.tail) {
            m[0.0 == p] <- -Inf
            m[1.0 == p] <- m_upper
            p_max <- vismeteor::pvmideal(m_lower, lm, psi, lower.tail = TRUE, perception_fun = perception_fun)

            idx <- p > 0.0 & p < p_max
            if (any(idx)) {
                p_max_ideal <- vismeteor::pmideal(m_lower - 0.5, psi, lower.tail = TRUE)
                m[idx] <- floor(0.5 + vismeteor::qmideal(
                    p_max_ideal * p[idx] / p_max,
                    psi,
                    lower.tail = TRUE
                ))
            }
            idx <- p >= p_max & p < 1.0
            if (any(idx)) {
                m_gen <- seq(m_lower, m_upper + 1)
                p_gen <- vismeteor::pvmideal(m_gen, lm, psi, lower.tail = lower.tail, perception_fun = perception_fun)
                p_idx <- findInterval(p[idx], p_gen, left.open = FALSE)
                m[idx] <- m_gen[p_idx]
            }
        } else {
            m[0.0 == p] <- m_upper
            m[1.0 == p] <- -Inf
            p_max <- vismeteor::pvmideal(m_lower, lm, psi, lower.tail = FALSE, perception_fun = perception_fun)

            idx <- p > p_max & p < 1.0
            if (any(idx)) {
                p_max_ideal <- vismeteor::pmideal(m_lower - 0.5, psi, lower.tail = FALSE)
                m[idx] <- floor(0.5 + vismeteor::qmideal(
                    1.0 - (1.0 - p[idx]) * ((1.0 - p_max_ideal) / (1.0 - p_max)),
                    psi,
                    lower.tail = FALSE
                ))
            }

            idx <- p > 0.0 & p <= p_max
            if (any(idx)) {
                m_gen <- seq(m_upper + 1, m_lower)
                p_gen <- vismeteor::pvmideal(m_gen, lm, psi, lower.tail = FALSE, perception_fun = perception_fun)
                p_idx <- findInterval(p[idx], p_gen, left.open = TRUE) + 1
                m[idx] <- m_gen[p_idx]
            }
        }

        m
    }

    arg_data <- data.frame(
        p = p,
        lm = lm,
        psi = psi
    )

    data_f <- as.factor(paste0(lm, "/", psi))
    data_s <- split(arg_data, data_f)
    m <- lapply(data_s, \(data) {
        p <- data$p
        lm <- data$lm[1]
        psi <- data$psi[1]
        f_q(p, lm, psi)
    })

    m <- unsplit(m, data_f)

    if (anyNA(m)) {
        warning("NaNs produced")
    }

    m
}

#' @rdname vmideal
#' @export
rvmideal <- function(n, lm, psi, perception_fun = vmperception) {
    if (anyNA(lm) || anyNA(psi)) {
        stop("NA's are not allowed!")
    }

    if (any(is.infinite(lm))) {
        stop("Infinite limiting magnitudes are not allowed!")
    }

    if (any(is.infinite(psi))) {
        stop("Infinite psi values are not allowed!")
    }

    perception_fun <- match.fun(perception_fun)

    p <- stats::runif(n)
    m <- rep(NA, n)

    idx <- p < 0.5
    if (any(idx)) {
        m[idx] <- vismeteor::qvmideal(p[idx], lm, psi, lower.tail = TRUE, perception_fun = perception_fun)
    }

    if (any(!idx)) {
        m[!idx] <- vismeteor::qvmideal(1.0 - p[!idx], lm, psi, lower.tail = FALSE, perception_fun = perception_fun)
    }

    m
}

#' @rdname vmideal
#' @export
cvmideal <- function(lm, psi, log = FALSE, perception_fun = vmperception) {
    if (anyNA(lm) || anyNA(psi)) {
        stop("NA's are not allowed!")
    }

    perception_fun <- match.fun(perception_fun)

    if (1 == length(psi)) {
        psi <- rep(psi, length(lm))
    }

    # Integration - similar to .vmideal_norm()
    f_integrate <- function(lm, psi) {
        m_max <- 15L
        m_upper <- .vmideal_upper_lm(lm)
        m_lower <- m_upper - m_max
        m <- as.integer(seq(m_lower, m_upper))
        p <- .dmideal_int(m, psi) * perception_fun(lm - m)
        p_lower_tail <- vismeteor::pmideal(m_lower - 0.5, psi, lower.tail = TRUE)
        sum(p) + p_lower_tail
    }

    p <- mapply(\(lm, psi) {
        if (Inf == lm & Inf == psi) {
            return(NA)
        }
        if (-Inf == lm & -Inf == psi) {
            return(NA)
        }
        if (Inf == lm & -Inf == psi) {
            return(if (log) 0.0 else 1.0)
        }
        if (Inf == psi & -Inf == lm) {
            return(if (log) -Inf else 0.0)
        }

        if (log) {
            base::log(f_integrate(lm, psi))
        } else {
            f_integrate(lm, psi)
        }
    }, lm, psi, SIMPLIFY = TRUE)

    if (anyNA(p)) {
        warning("NaNs produced")
    }

    p
}

#' upper available magnitude
#'
#' @noRd
.vmideal_upper_lm <- function(lm) {
    lm_round <- round(lm)
    offset <- lm - lm_round
    if (-0.5 == offset) {
        lm_round <- lm_round - 1L
    }

    lm_round
}

#' density of ideal integer magnitude distribution
#'
#' @noRd
.dmideal_int <- function(m, psi, log = FALSE) {
    psi_exp <- 10.0
    r_lower <- 10^0.4

    a <- -base::log(r_lower)
    d <- rep(NA, length(m))

    idx <- m > (psi + psi_exp)
    if (any(idx)) {
        if (log) {
            d[idx] <- base::log(1 - base::exp(1.5 * a)) + a * (1.5 * (m[idx] - psi) - 0.75)
        } else {
            d[idx] <- base::exp(a * (1.5 * (m[idx] - psi) - 0.75)) -
                base::exp(a * (1.5 * (m[idx] - psi) + 0.75))
        }
    }

    idx <- m < (psi - psi_exp)
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
.vmideal_norm <- function(lm, psi, perception_fun) {
    m_max <- 15L
    m_upper <- .vmideal_upper_lm(lm)
    m_lower <- m_upper - m_max
    m <- as.integer(seq(m_lower, m_upper))
    p <- .dmideal_int(m, psi) * perception_fun(lm - m)
    names(p) <- as.character(m)
    p_lower_tail <- vismeteor::pmideal(m_lower - 0.5, psi, lower.tail = TRUE)
    norm <- sum(p) + p_lower_tail

    list(
        m_lower = m_lower,
        m_upper = m_upper,
        norm = norm,
        p = p / norm,
        p_lower_tail = p_lower_tail / norm
    )
}
