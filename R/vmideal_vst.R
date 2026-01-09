#' @name vmidealVst
#' @aliases vmidealVstFromMagn
#' @aliases vmidealVstToPsi
#' @title Variance-stabilizing Transformation for the Ideal Distribution of Visual Meteor Magnitudes
#'
#' @description
#' Applies a variance-stabilizing transformation to meteor magnitudes
#' under the assumption of the ideal magnitude distribution.
#'
#' @param m integer; the meteor magnitude.
#' @param lm numeric; limiting magnitude.
#' @param tm numeric; transformed magnitude.
#' @param deriv.degree integer; the degree of the derivative at `tm` to return
#'   instead of `psi`. Must be `0`, `1` or `2`.
#'
#' @details
#' Many linear models require the variance of visual meteor magnitudes to be
#' homoscedastic. The function `vmidealVstFromMagn` applies a transformation
#' that produces homoscedastic distributions of visual meteor magnitudes if the
#' underlying magnitudes follow the ideal magnitude distribution.
#' In this sense, the transformation acts as a normalization of meteor magnitudes
#' and yields a variance close to `1.0`.
#'
#' The ideal distribution of visual meteor magnitudes
#' depends on the [parameter psi][vismeteor::vmideal] and the limiting magnitude `lm`,
#' resulting in a two-parameter distribution. Without detection probabilities,
#' the magnitude distribution reduces to a pure [ideal magnitude distribution][vismeteor::mideal],
#' which depends only on the parameter `psi`. Since the
#' limiting magnitude `lm` is a fixed parameter and never estimated
#' statistically, the magnitudes can be transformed such that, for example,
#' the mean of the transformed magnitudes directly provides an estimate of `psi`
#' using the function `vmidealVstToPsi`.
#'
#' This transformation is valid for \eqn{-10 \le \texttt{psi} \le 9}.
#' The numerical form of the transformation is version-specific and may change
#' substantially in future releases. Do not rely on equality of transformed
#' values across package versions.
#'
#' @return
#' * `vmidealVstFromMagn`: a numeric value, the transformed meteor magnitude.
#' * `vmidealVstToPsi`: a numeric value of the parameter `psi`, derived from
#'   the mean of `tm`.
#' The argument `deriv.degree` can be used to apply the delta method.
#'
#' @note
#' The internal approximations used here are derived from the perception
#' probabilities produced by [vismeteor::vmperception].
#' For details on the derivation, see the script `inst/derivation/vmideal_vst.R` in the
#' package's source code.
#'
#' @seealso [vismeteor::vmideal] [vismeteor::mideal] [vismeteor::vmperception]
#' @examples
#' N <- 100
#' psi <- 5.0
#' limmag <- 6.3
#'
#' # Simulate magnitudes
#' m <- rvmideal(N, limmag, psi)
#'
#' # Variance-stabilizing transformation
#' tm <- vmidealVstFromMagn(m, limmag)
#' tm.mean <- mean(tm)
#' tm.var  <- var(tm)
#'
#' # Estimator for psi from the transformed mean
#' psi.hat  <- vmidealVstToPsi(tm.mean, limmag)
#'
#' # Derivative d(psi)/d(tm) at tm.mean (needed for the delta method)
#' dpsi_dtm <- vmidealVstToPsi(tm.mean, limmag, deriv.degree = 1L)
#'
#' # Variance of the sample mean of tm
#' var_tm.mean <- tm.var / N
#'
#' # Delta method: variance and standard error of psi.hat
#' var_psi.hat <- (dpsi_dtm^2) * var_tm.mean
#' se_psi.hat  <- sqrt(var_psi.hat)
#'
#' # Results
#' print(psi.hat)
#' print(se_psi.hat)

#' @rdname vmidealVst
#' @export
vmidealVstFromMagn <- function(m, lm) {
    offset <- lm - round(lm)
    if (1L == length(lm)) {
        limmag <- rep(lm, length(m))
        offset <- rep(offset, length(m))
    }

    arg.data <- data.frame(
        x = lm - m,
        offset = offset
    )

    data.f <- as.factor(arg.data$offset)
    data.s <- split(arg.data, data.f)
    y <- lapply(data.s, function(data) {
        x <- data$x
        offset <- data$offset[1]
        params <- .vmidealVstFromMagn.params
        sy <- c(
            params$p1.fun(offset),
            params$p2.fun(offset),
            params$p3.fun(offset),
            params$p4.fun(offset),
            params$p5.fun(offset),
            params$p6.fun(offset)
        )
        sx <- params$sx
        f.spline <- stats::splinefun(sx, sy)
        y <- rep(NA, length(x))

        idx <- x < sx[1]
        if (any(idx)) {
            y[idx] <- f.spline(sx[1], deriv = 1) * (x[idx] - sx[1]) + sy[1]
        }

        idx <- x > sx[length(sx)]
        if (any(idx)) {
            y[idx] <- f.spline(sx[length(sx)], deriv = 1) * (x[idx] - sx[length(sx)]) + sy[length(sy)]
        }

        idx <- x >= sx[1] & x <= sx[length(sx)]
        if (any(idx)) {
            y[idx] <- f.spline(x[idx])
        }

        y - params$p0.fun(offset)
    })

    unsplit(y, data.f)
}

#' @rdname vmidealVst
#' @export
vmidealVstToPsi <- function(tm, lm, deriv.degree = 0L) {
    poly.coef0 <- c(-2.97086442804517, -2.72858043751615, -0.683184284791628, -0.1971973188227,
        -0.0707494993259986, -0.0267813318470729, -0.00482163891332698,
        -1.03296948649241e-05, 6.01391927719714e-05
    )
    names(poly.coef0) <- seq_along(poly.coef0) - 1 # exponents

    # x min 0.016 (psi approx 9 at limiting maginitde of 5.5)
    # x max 8.22(psi approx -10 at limiting maginitde of 6.5)
    tm[tm < 0.016 | tm > 8.22] <- NA
    y <- log(tm)

    if(deriv.degree > 0L) {
        poly.coef1 <- f.polynomial.coef(poly.coef0, deriv.degree = 1L)
    }
    if(deriv.degree > 1L) {
        poly.coef2 <- f.polynomial.coef(poly.coef1, deriv.degree = 1L)
    }
    if(deriv.degree > 2L) {
       stop(paste('deriv.degree', deriv.degree, 'not implemented!'))
    }

    if(2L == deriv.degree) {
        (f.polynomial(y, poly.coef2) - f.polynomial(y, poly.coef1))/tm^2
    } else if(1L == deriv.degree) {
        f.polynomial(y, poly.coef1)/tm
    } else {
        lm + f.polynomial(y, poly.coef0)
    }
}

#' @keywords internal
.vmidealVstFromMagn.params <- (function() {
    param.df <- data.frame(
        'p1' = c(0.379578706193683, 0.380865978506213,
        0.383646581863156, 0.386594955357005, 0.384516871380943, 0.378119364136071,
        0.39783016804395, 0.373393379287508, 0.366393502192895, 0.345593857555874,
        0.367086720870886, 0.359195376100079, 0.347082537633701, 0.312420722892609,
        0.312232925285191), 'p2' = c(0.758955184180252, 0.757434042449928,
        0.756545241530695, 0.754597324053272, 0.746096530130275, 0.736273688899513,
        0.754472950700777, 0.729839979801878, 0.724090251762525, 0.706272163875503,
        0.732424006749156, 0.730856947686901, 0.722541149425633, 0.690367122987635,
        0.691926213129086), 'p3' = c(1.89943177726206, 1.89844069948433,
        1.89825430762927, 1.89725244849407, 1.88997593310439, 1.88072563899649,
        1.89899661183851, 1.87396958675163, 1.86749053248022, 1.84882849631795,
        1.87410789172307, 1.87172629242088, 1.86305301012132, 1.83068592572223,
        1.83212684285442), 'p4' = c(3.27014476551314, 3.26892787903128,
        3.26844009590108, 3.26702893292816, 3.25921320867758, 3.24969559483436,
        3.26789784387066, 3.24296948683872, 3.23673983810267, 3.21843295052407,
        3.24413734533724, 3.24221462073629, 3.23377186059092, 3.20154167652271,
        3.20307252532895), 'p5' = c(4.50424034207103, 4.50305721713432,
        4.50261471059499, 4.5012654394704, 4.49353145072154, 4.48405231361262,
        4.50225758259143, 4.47729738300568, 4.47099710887246, 4.45259464315401,
        4.47819774376931, 4.4761817895577, 4.46769920746914, 4.43544755919314,
        4.43696568503357), 'p6' = c(5.65036712145262, 5.64927790140628,
        5.64896064069461, 5.6477824229466, 5.64028041416121, 5.63092781192841,
        5.64918249953455, 5.62420382559607, 5.61781703615347, 5.5992741104315,
        5.62469988859743, 5.62248449312786, 5.61389842883615, 5.58158409117071,
        5.58306055049816), 'offset' = c(-0.5, -0.48, -0.45, -0.4, -0.3,
        -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.48, 0.5), 'intercept' = c(1.4811934820941,
        1.47927325731333, 1.47783119789222, 1.47503283992891, 1.46516880013637,
        1.45440969975367, 1.47200033556686, 1.44700869631271, 1.44106727057982,
        1.42347451350523, 1.45029075335968, 1.44986985073928, 1.44229503571524,
        1.41062845239295, 1.41258415137388)
    )

    list(
        sx = c(1.0, 2.0, 4.0, 6.0, 8.0, 10.0),
        p0.fun = stats::approxfun(param.df$offset, param.df[['intercept']]),
        p1.fun = stats::approxfun(param.df$offset, param.df[['p1']]),
        p2.fun = stats::approxfun(param.df$offset, param.df[['p2']]),
        p3.fun = stats::approxfun(param.df$offset, param.df[['p3']]),
        p4.fun = stats::approxfun(param.df$offset, param.df[['p4']]),
        p5.fun = stats::approxfun(param.df$offset, param.df[['p5']]),
        p6.fun = stats::approxfun(param.df$offset, param.df[['p6']])
    )
})()
