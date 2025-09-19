#' @name vmgeom
#' @aliases dvmgeom
#' @aliases pvmgeom
#' @aliases qvmgeom
#' @aliases rvmgeom
#' @title Geometric Model of Visual Meteor Magnitudes
#' @description
#' Density, distribution function, quantile function, and random generation
#' for the geometric model of visual meteor magnitudes.
#' @param m numeric; the meteor magnitude.
#' @param p numeric; probability.
#' @param lm numeric; limiting magnitude.
#' @param r numeric; the population index.
#' @param n numeric; count of meteor magnitudes.
#' @param log logical; if `TRUE`, probabilities p are given as `log(p)`.
#' @param lower.tail logical; if `TRUE` (default) probabilities are
#'     \eqn{P[M < m]}, otherwise, \eqn{P[M \ge m]}.
#' @param perception.fun function; perception probability function (optional).
#'     Default is [vismeteor::vmperception].
#' @details
#' In visual meteor observations, magnitudes are usually estimated as integer values.
#' Hence, this distribution is discrete and its probability mass function is
#' \deqn{
#'     {\displaystyle P[X = x] \sim f(x) \, \mathrm r^{-x}} \,\mathrm{,}
#' }
#' where \eqn{x \ge -0.5} denotes the difference between the limiting magnitude `lm`
#' and the meteor magnitude `m`, and \eqn{f(x)} is the perception probability function.
#' Thus, the distribution is the product of the
#' [perception probabilities][vismeteor::vmperception] and the
#' underlying [geometric distribution][stats::Geometric] of meteor magnitudes.
#' Therefore, the parameter `p` of the geometric distribution is given by `p = 1 - 1/r`.
#'
#' The parameter `lm` specifies the reference for the meteor magnitude `m`.
#' `m` must be an integer meteor magnitude.  
#' The length of the vector `lm` must either equal the length of the vector `m`,
#' or `lm` must be a scalar value.  
#' In the case of `rvmgeom`, the length of the vector `lm` must equal `n`,
#' or `lm` must be a scalar value.
#'
#' If a perception probability function `perception.fun` is provided,
#' it must have the signature `function(x)` and return the perception probability of
#' the difference `x` between the limiting magnitude and the meteor magnitude.
#' If `x >= 15.0`, the function `perception.fun` should return a perception probability of `1.0`.
#' The argument `perception.fun` is resolved using [match.fun].
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
#'     -sum(freq * dvmgeom(m, limmag, r, log=TRUE))
#' }
#'
#' # maximum likelihood estimation (MLE) of r
#' est <- optim(2, lld, method='Brent', lower=1.1, upper=4)
#'
#' # estimations
#' est$par # mean of r
#'
#' # generate random meteor magnitudes
#' m <- rvmgeom(N, r, lm=limmag)
#'
#' # log likelihood function
#' llr <- function(r) {
#'     -sum(dvmgeom(m, limmag, r, log=TRUE))
#' }
#'
#' # maximum likelihood estimation (MLE) of r
#' est <- optim(2, llr, method='Brent', lower=1.1, upper=4, hessian=TRUE)
#'
#' # estimations
#' est$par # mean of r
#' sqrt(1/est$hessian[1][1]) # standard deviation of r
#'
#' m <- seq(6, -4, -1)
#' p <- vismeteor::dvmgeom(m, limmag, r)
#' barplot(
#'     p,
#'     names.arg = m,
#'     main = paste0('Density (r = ', r, ', limmag = ', limmag, ')'),
#'     col = "blue",
#'     xlab = 'm',
#'     ylab = 'p',
#'     border = "blue",
#'     space = 0.5
#' )
#' axis(side = 2, at = pretty(p))


#' @rdname vmgeom
#' @export
dvmgeom <- function(m, lm, r, log = FALSE, perception.fun = NULL) {
    if (anyNA(m) | anyNA(lm) | anyNA(r)) {
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

    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) all(is.infinite(x) | abs(x - round(x)) < tol)
    if (! is.wholenumber(m)) {
        stop("magnitudes must be integer values!")
    }

    if (is.null(perception.fun)) {
        perception.fun <- vismeteor::vmperception
    } else {
        perception.fun <- match.fun(perception.fun)
    }

    p.geom <- 1.0 - 1.0/r

    if (1 == length(r)) {
        r <- rep(r, length(m))
        p.geom <- rep(p.geom, length(m))
    }

    std_res <- vmgeom.std(m, lm)
    m <- std_res$m
    offset <- std_res$offset

    f.density <- function(m, offset, p.geom) {
        m.max <- 15L
        idx <- m <= m.max
        d <- stats::dgeom(m, p.geom, log = TRUE)
        if (any(idx)) {
            d[idx] <- d[idx] + base::log(perception.fun(m[idx] + offset))
        }

        d - base::log(vmgeom.norm(offset, p.geom, m.max, perception.fun))
    }

    arg.data <- data.frame(
        m = m,
        offset = offset,
        p.geom = p.geom
    )

    data.f <- as.factor(paste0(offset, '/', p.geom))
    data.s <- split(arg.data, data.f)
    d <- lapply(data.s, function(data) {
        m <- data$m
        offset <- data$offset[1]
        p.geom <- data$p.geom[1]
        d <- rep(-Inf, length(m))

        idx <- m > -1
        if (any(idx)) {
            d[idx] <- f.density(m[idx], offset, p.geom)
        }

        if (! log) {
            d[idx] <- exp(d[idx])
            d[!idx] <- 0.0
        }

        d
    })

    unsplit(d, data.f)
}

#' @rdname vmgeom
#' @export
pvmgeom <- function(m, lm, r, lower.tail = TRUE, log = FALSE, perception.fun = NULL) {
    if (anyNA(m) | anyNA(lm) | anyNA(r)) {
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

    is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) all(is.infinite(x) | abs(x - round(x)) < tol)
    if (! is.wholenumber(m)) {
        stop("magnitudes must be integer values!")
    }

    if (is.null(perception.fun)) {
        perception.fun <- vismeteor::vmperception
    } else {
        perception.fun <- match.fun(perception.fun)
    }

    p.geom <- 1.0 - 1.0/r

    if (1 == length(r)) {
        r <- rep(r, length(m))
        p.geom <- rep(p.geom, length(m))
    }

    std_res <- vmgeom.std(m, lm)
    m <- std_res$m
    offset <- std_res$offset

    f.density <- function(m, offset, p.geom) {
        stats::dgeom(m, p.geom) * perception.fun(m + offset)
    }

    f.sum <- Vectorize(function(m, offset, p.geom) {
        m <- as.integer(seq(0, m))
        sum(f.density(m, offset, p.geom))
    })

    f.prob <- function(m, offset, p.geom) {
        m.max <- 15L
        norm <- vmgeom.norm(offset, p.geom, m.max, perception.fun)
        p <- rep(0.0, length(m))

        if (lower.tail) {
            idx <- m <= m.max
            if(any(idx)) {
                p[idx] <- 1.0 - f.sum(m[idx], offset, p.geom)/norm
            }

            idx <- m > m.max
            if(any(idx)) {
                p[idx] <- stats::pgeom(m[idx], p.geom, lower.tail = FALSE)/norm
            }
        } else {
            idx <- m <= m.max
            if(any(idx)) {
                p[idx] <- f.sum(m[idx], offset, p.geom)/norm
            }

            idx <- m > m.max
            if(any(idx)) {
                p[idx] <- 1.0 - stats::pgeom(m[idx], p.geom, lower.tail = FALSE)/norm
            }
        }

        p
    }

    arg.data <- data.frame(
        m = m,
        offset = offset,
        p.geom = p.geom
    )

    data.f <- as.factor(paste0(offset, '/', p.geom))
    data.s <- split(arg.data, data.f)
    p <- lapply(data.s, function(data) {
        m <- data$m
        offset <- data$offset[1]
        p.geom <- data$p.geom[1]

        if (lower.tail) {
            p <- rep(1.0, length(m))
        } else {
            p <- rep(0.0, length(m))
        }

        idx <- m > -1
        if(any(idx)) {
            p[idx] <- f.prob(m[idx], offset, p.geom)
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

    unsplit(p, data.f)
}

#' @rdname vmgeom
#' @export
qvmgeom <- function(p, lm, r, lower.tail = TRUE, perception.fun = NULL) {
    if (anyNA(p) | anyNA(lm) | anyNA(r)) {
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

    if (is.null(perception.fun)) {
        perception.fun <- vismeteor::vmperception
    } else {
        perception.fun <- match.fun(perception.fun)
    }

    p.geom <- 1.0 - 1.0/r

    if (1 == length(r)) {
        r <- rep(r, length(p))
        p.geom <- rep(p.geom, length(p))
    }

    lm.round <- round(lm)
    offset <- lm - lm.round
    lm <- lm.round

    idx <- -0.5 == offset
    lm[idx] <- lm[idx] - 1L
    offset[idx] <- offset[idx] + 1.0

    if (1 == length(offset)) {
        offset <- rep(offset, length(p))
    }

    arg.data <- data.frame(
        p = p,
        offset = offset,
        p.geom = p.geom,
        r = r
    )

    data.f <- as.factor(paste0(offset, '/', p.geom))
    data.s <- split(arg.data, data.f)
    m <- lapply(data.s, function(data) {
        m.max <- 15L
        p <- data$p
        offset <- data$offset[1]
        p.geom <- data$p.geom[1]
        r <- data$r[1]
        m <- rep(NA, length(p))

        if(lower.tail) {
            m[1.0 == p] <- 0L
            p.max <- 1.0 - vismeteor::pvmgeom(0, m.max + offset, r, lower.tail = FALSE, perception.fun = perception.fun)
            idx <- p>=0.0 & p<p.max
            if (any(idx)) {
                m[idx] <- m.max + 1 + stats::qgeom(p[idx]/p.max, p.geom, lower.tail = FALSE)
            }

            idx <- p>=p.max & p<1.0
            if (any(idx)) {
                m0 <- seq(-m.max, 0, 1)
                p0 <- c(vismeteor::pvmgeom(m0, offset, r, lower.tail = TRUE, perception.fun = perception.fun), 1.0)
                m0 <- c(m0, 1L)
                p.idx <- findInterval(p[idx], p0, left.open = FALSE)
                m[idx] <- -m0[p.idx]
            }
        } else {
            m[0.0 == p] <- 0L
            p.max <- vismeteor::pvmgeom(0, m.max + offset, r, lower.tail = FALSE, perception.fun = perception.fun)
            idx <- p>p.max & p<=1.0
            if (any(idx)) {
                m[idx] <- m.max + 1L + stats::qgeom((p[idx] - p.max)/(1.0 - p.max), p.geom)
            }

            idx <- p>0.0 & p<=p.max
            if (any(idx)) {
                m0 <- seq(1, -m.max, -1)
                p0 <- vismeteor::pvmgeom(m0, offset, r, lower.tail = FALSE, perception.fun = perception.fun)
                p.idx <- findInterval(p[idx], p0, left.open = TRUE) + 1
                m[idx] <- -m0[p.idx]
                m[m<0] <- NA
            }
        }

        m
    })

    m <- lm - unsplit(m, data.f)

    if (anyNA(m)) {
        warning('NaNs produced')
    }

    as.numeric(m)
}

#' @rdname vmgeom
#' @export
rvmgeom <- function(n, lm, r, perception.fun = NULL) {
    if (anyNA(lm) | anyNA(r)) {
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

    if (is.null(perception.fun)) {
        perception.fun <- vismeteor::vmperception
    } else {
        perception.fun <- match.fun(perception.fun)
    }

    p <- stats::runif(n)
    m <- rep(NA, n)
    idx <- p < 0.5

    if (any(idx)) {
        m[idx] <- vismeteor::qvmgeom(p[idx], lm, r, lower.tail = TRUE, perception.fun = perception.fun)
    }

    if (any(!idx)) {
        m[!idx] <- vismeteor::qvmgeom(1.0 - p[!idx], lm, r, lower.tail = FALSE, perception.fun = perception.fun)
    }

    m
}

#' standardization
#'
#' @noRd
vmgeom.std <- function(m, lm) {
    m <- lm - m
    m.round <- round(m)
    offset <- rep(0.0, length(m))
    idx <- is.infinite(m)
    offset[! idx] <- m[! idx] - m.round[! idx]
    m <- m.round

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
vmgeom.norm <- function(offset, p.geom, m.max, perception.fun){
    m <- as.integer(seq(0, m.max))
    sum(stats::dgeom(m, p.geom) * perception.fun(m + offset)) +
        stats::pgeom(m.max, p.geom, lower.tail = FALSE)
}
