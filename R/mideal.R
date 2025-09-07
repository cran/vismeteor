#' @name mideal
#' @aliases dmideal
#' @title Ideal Distribution of Meteor Magnitudes
#' @description
#' Density, distribution function, quantile function and random generation
#' for the ideal distribution of meteor magnitudes.
#' @param psi numeric; the location parameter of a probability distribution.
#'     It is the only parameter of the distribution.
#' @param m numeric; meteor magnitude.
#' @param p numeric; probability.
#' @param n numeric; count of meteor magnitudes.
#' @param log logical; if `TRUE`, probabilities p are given as `log(p)`.
#' @param lower.tail logical; if `TRUE` (default) probabilities are
#'     \eqn{P[M \le m]}, otherwise, \eqn{P[M > m]}.
#' @details
#' The density of the ideal distribution of meteor magnitudes is
#' \deqn{
#'     {\displaystyle \frac{\mathrm{d}p}{\mathrm{d}m} = \frac{3}{2} \, \log(r) \sqrt{\frac{r^{3 \, \psi + 2 \, m}}{(r^\psi + r^m)^5}}}
#' }
#' where \eqn{m} is the meteor magnitude, \eqn{r = 10^{0.4} \approx 2.51189 \dots} is a constant and
#' \eqn{\psi} is the only parameter of this magnitude distribution.
#' @return
#' `dmideal` gives the density, `pmideal` gives the distribution function,
#' `qmideal` gives the quantile function and `rmideal` generates random deviates.
#'
#' The length of the result is determined by `n` for `rmideal`, and is the maximum
#' of the lengths of the numerical vector arguments for the other functions.
#'
#' `qmideal` can return `NaN` value with a warning.
#' @references Richter, J. (2018) \emph{About the mass and magnitude distributions of meteor showers}.
#'   WGN, Journal of the International Meteor Organization, vol. 46, no. 1, p. 34-38
#' @examples
#' old_par <- par(mfrow = c(2,2))
#' psi <- 5.0
#' plot(
#'     function(m) dmideal(m, psi, log = FALSE),
#'     -5, 10,
#'     main = paste0('Density of the Ideal Meteor Magnitude\nDistribution (psi = ', psi, ')'),
#'     col = "blue",
#'     xlab = 'm',
#'     ylab = 'dp/dm'
#' )
#' abline(v=psi, col="red")
#'
#' plot(
#'     function(m) dmideal(m, psi, log = TRUE),
#'     -5, 10,
#'     main = paste0('Density of the Ideal Meteor Magnitude\nDistribution (psi = ', psi, ')'),
#'     col = "blue",
#'     xlab = 'm',
#'     ylab = 'log( dp/dm )'
#' )
#' abline(v=psi, col="red")
#'
#' plot(
#'     function(m) pmideal(m, psi),
#'     -5, 10,
#'     main = paste0('Probability of the Ideal Meteor Magnitude\nDistribution (psi = ', psi, ')'),
#'     col = "blue",
#'     xlab = 'm',
#'     ylab = 'p'
#' )
#' abline(v=psi, col="red")
#'
#' plot(
#'     function(p) qmideal(p, psi),
#'     0.01, 0.99,
#'     main = paste('Quantile of the Ideal Meteor Magnitude\nDistribution (psi = ', psi, ')'),
#'     col = "blue",
#'     xlab = 'p',
#'     ylab = 'm'
#' )
#' abline(h=psi, col="red")
#'
#' # generate random meteor magnitudes
#' m <- rmideal(1000, psi)
#'
#' # log likelihood function
#' llr <- function(psi) {
#'     -sum(dmideal(m, psi, log=TRUE))
#' }
#'
#' # maximum likelihood estimation (MLE) of psi
#' est <- optim(2, llr, method='Brent', lower=0, upper=8, hessian=TRUE)
#'
#' # estimations
#' est$par # mean of psi
#' sqrt(1/est$hessian[1][1]) # standard deviation of psi
#'
#' par(old_par)

#' @rdname mideal
#' @export
dmideal <- function(m, psi = 0.0, log = FALSE) {
    a <- -base::log(10.0)/2.5
    d <- rep(NA, length(m))
    psi.exp <- 10.0
    m <- m - psi

    idx <- m > psi.exp
    if (any(idx)) {
        if (log) {
            d[idx] <- 1.5 * a * m[idx]
        } else {
            d[idx] <- exp(1.5 * a * m[idx])
        }
    }

    idx <- m < -psi.exp
    if (any(idx)) {
        if (log) {
            d[idx] <- -a * m[idx]
        } else {
            d[idx] <- exp(-a * m[idx])
        }
    }

    idx <- is.na(d)
    if (any(idx)) {
        if (log) {
            d[idx] <- (3 * a * m[idx] - 5 * base::log(1.0 + exp(a * m[idx])))/2
        } else {
            d[idx] <- base::sqrt(
                exp(3 * a * m[idx])/(1.0 + exp(a * m[idx]))^5
            )
        }
    }

    if (log) {
        base::log(-1.5 * a) + d
    } else {
        -1.5 * a * d
    }
}

#' @rdname mideal
#' @export
pmideal <- function(m, psi = 0.0, lower.tail = TRUE, log = FALSE) {
    a <- -base::log(10.0)/2.5
    psi.exp <- 10.0

    m <- m - psi
    p <- rep(NA, length(m))

    if (lower.tail) {
        spline.knods <- c(
            -8.80472534014006,
            -8.34423787637494,
            -7.88372873094501,
            -7.42318695592341,
            -6.96263200116049,
            -6.50203076601891,
            -6.04137424073419,
            -5.58062400802751,
            -5.11973760793805,
            -4.6586468529202,
            -4.19721109282873,
            -3.73524195963578,
            -3.27243690786416,
            -2.8083150258617,
            -2.34214577748568,
            -1.87281940429522,
            -1.39869703344889,
            -0.917450131705818,
            -0.425957810918723,
            0.0796223098768092,
            0.603456223455557,
            1.14927014152791,
            1.71941908808972,
            2.31420877744153,
            2.93185193543927,
            3.56905016441712,
            4.22186644876613,
            4.88646802251895,
            5.55957886852415,
            6.23863980114079,
            6.92171415166587,
            7.6074921875779,
            8.2951330347783,
            8.98380534340058,
            9.67360021970155,
            10.3632095823678,
            11.0538798267455,
            11.7465047916797,
            12.4480963334632,
            13.1258002217638,
            13.8198202144629
        )
    } else {
        spline.knods <- c(
            8.80485039014547,
            8.34431865416611,
            7.88377975423457,
            7.42322504002181,
            6.96264953286627,
            6.50204026882691,
            6.04137719922142,
            5.5806286403783,
            5.11974504501854,
            4.65864760617632,
            4.19721233788351,
            3.73524403458982,
            3.27243684618656,
            2.80831497655707,
            2.34214575402348,
            1.87281941733144,
            1.3986971549163,
            0.917450465736558,
            0.425957597606353,
            -0.0796226279169917,
            -0.603456102786167,
            -1.14927012452963,
            -1.71941902164115,
            -2.31420847179972,
            -2.93185085210491,
            -3.56905007596158,
            -4.22186620307173,
            -4.88646699897379,
            -5.55957402041309,
            -6.23861581194743,
            -6.92169294002751,
            -7.60746668941095,
            -8.295014234359,
            -8.98371489435258,
            -9.67316619598767,
            -10.3630974814401,
            -11.053325888328,
            -11.7437374092862,
            -12.4343089812847,
            -13.1248984689238,
            -13.8155074574103
        )
    }
    f.spline <- stats::splinefun(seq(-psi.exp, psi.exp, 0.5), spline.knods, method = "hyman")

    idx <- lower.tail & m < -psi.exp
    if (any(idx)) {
        p.max <- 1/(1 + exp(-f.spline(-psi.exp)))
        if (log) {
            p[idx] <- base::log(p.max) + stats::pexp(-psi.exp - m[idx], -a, lower.tail = FALSE, log = TRUE)
        } else {
            p[idx] <- p.max * stats::pexp(-psi.exp - m[idx], -a, lower.tail = FALSE, log = FALSE)
        }
    }

    idx <- !lower.tail & m > psi.exp
    if (any(idx)) {
        p.max <- 1/(1 + exp(-f.spline(psi.exp)))
        if (log) {
            p[idx] <- base::log(p.max) + stats::pexp(m[idx] - psi.exp, -1.5 * a, lower.tail = FALSE, log = TRUE)
        } else {
            p[idx] <- p.max * stats::pexp(m[idx] - psi.exp, -1.5 * a, lower.tail = FALSE, log = FALSE)
        }
    }

    idx <- is.na(p)
    if (any(idx)) {
        p[idx] <- 1/(1 + exp(-f.spline(m[idx])))
        if (log) {
            p[idx] <- base::log(p[idx])
        }
    }

    p
}

#' @rdname mideal
#' @export
qmideal <- function(p, psi = 0.0, lower.tail = TRUE) {
    a <- -base::log(10.0)/2.5
    psi.exp <- 10.0

    m <- rep(NA, length(p))
    apply.idx <- !is.na(p) & p>0.0 & p < 1.0

    if (lower.tail) {
        m[0.0 == p] <- -Inf
        m[(p + 1e-07) >= 1.0 & p <= 1.0] <- Inf
        spline.knods <- c(
            -8.80472534014006,
            -8.34423787637494,
            -7.88372873094501,
            -7.42318695592341,
            -6.96263200116049,
            -6.50203076601891,
            -6.04137424073419,
            -5.58062400802751,
            -5.11973760793805,
            -4.6586468529202,
            -4.19721109282873,
            -3.73524195963578,
            -3.27243690786416,
            -2.8083150258617,
            -2.34214577748568,
            -1.87281940429522,
            -1.39869703344889,
            -0.917450131705818,
            -0.425957810918723,
            0.0796223098768092,
            0.603456223455557,
            1.14927014152791,
            1.71941908808972,
            2.31420877744153,
            2.93185193543927,
            3.56905016441712,
            4.22186644876613,
            4.88646802251895,
            5.55957886852415,
            6.23863980114079,
            6.92171415166587,
            7.6074921875779,
            8.2951330347783,
            8.98380534340058,
            9.67360021970155,
            10.3632095823678,
            11.0538798267455,
            11.7465047916797,
            12.4480963334632,
            13.1258002217638,
            13.8198202144629
        )
        p.min <- 1/(1 + exp(-spline.knods[1]))
    } else {
        m[0.0 == p] <- Inf
        m[(p + 1e-07) >= 1.0 & p <= 1.0] <- -Inf
        spline.knods <- c(
            8.80485039014547,
            8.34431865416611,
            7.88377975423457,
            7.42322504002181,
            6.96264953286627,
            6.50204026882691,
            6.04137719922142,
            5.5806286403783,
            5.11974504501854,
            4.65864760617632,
            4.19721233788351,
            3.73524403458982,
            3.27243684618656,
            2.80831497655707,
            2.34214575402348,
            1.87281941733144,
            1.3986971549163,
            0.917450465736558,
            0.425957597606353,
            -0.0796226279169917,
            -0.603456102786167,
            -1.14927012452963,
            -1.71941902164115,
            -2.31420847179972,
            -2.93185085210491,
            -3.56905007596158,
            -4.22186620307173,
            -4.88646699897379,
            -5.55957402041309,
            -6.23861581194743,
            -6.92169294002751,
            -7.60746668941095,
            -8.295014234359,
            -8.98371489435258,
            -9.67316619598767,
            -10.3630974814401,
            -11.053325888328,
            -11.7437374092862,
            -12.4343089812847,
            -13.1248984689238,
            -13.8155074574103
        )
        p.min <- 1/(1 + exp(-spline.knods[length(spline.knods)]))
    }

    f.spline <- stats::splinefun(spline.knods, seq(-psi.exp, psi.exp, 0.5), method = "hyman")

    idx <- apply.idx & lower.tail & p <= p.min
    if (any(idx)) {
        m[idx] <- -psi.exp - stats::qexp(p[idx]/p.min, -a, lower.tail = FALSE, log = FALSE)
    }

    idx <- apply.idx & !lower.tail & p <= p.min
    if (any(idx)) {
        m[idx] <- psi.exp + stats::qexp(p[idx]/p.min, -1.5 * a, lower.tail = FALSE, log = FALSE)
    }

    idx <- apply.idx & is.na(m)
    if (any(idx)) {
        p.logit <- base::log(p[idx]/(1-p[idx]))
        m[idx] <- f.spline(p.logit)
    }

    if (anyNA(m)) {
        warning('NaNs produced')
    }

    m <- m + psi
}

#' @rdname mideal
#' @export
rmideal <- function(n, psi = 0.0) {
    p <- stats::runif(n)
    m <- rep(NA, n)

    idx <- p < 0.5
    if (any(idx)) {
        m[idx] <- vismeteor::qmideal(p[idx], psi, lower.tail = TRUE)
    }

    if (any(!idx)) {
        m[!idx] <- vismeteor::qmideal(1.0 - p[!idx], psi, lower.tail = FALSE)
    }

    m
}
