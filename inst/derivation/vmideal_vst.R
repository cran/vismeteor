##
## Script for deriving/validating the variance-stabilizing transformation
## for the ideal magnitude model (vmideal).
##
## Purpose
## - Provide a smooth approximation t(m; limmag) for vmideal using a
##   spline in x = limmag − m. The spline shape depends on the fractional
##   component (offset) of limmag.
## - Choose parameters so that Var[t] under vmideal is approximately constant
##   across a wide range of psi, and align an intercept to a convenient
##   reference (via regression/optimization).
## - Fit a regression mapping from log(E[t]) to psi − limmag for quick
##   inversion (myVmidealVstToPsi).
##
## Contents
## - myVmidealVstFromMagn(): variance‑stabilizing transform t for vmideal,
##   built from a cubic spline with anchors at sx = {1, 2, 4, 6, 8, 10}.
##   Coefficients depend on limmag offset via linear interpolation.
## - param.df: precomputed coefficients for offsets derived from limmag in
##   [5.5, 6.5]; columns '1'..'6' are spline values at sx; 'intercept' is the
##   alignment term; 'offset' is limmag − round(limmag).
## - data.t.fun(): generates E[t] and Var[t] over grids of psi and limmag.
## - Optimization blocks (guarded by if(FALSE)) to refine coefficients.
## - myVmidealVstToPsi(): polynomial mapping from log(x) with x = E[t] to
##   psi, including optional derivative for Jacobians.
##
## Usage / workflow
## 1. Execute the script sequentially. Activate the `if (FALSE)` sections only
##    when you need to recompute spline anchors or intercepts; they can be slow.
## 2. Copy the updated `param.df` (columns `1`..`6`, `intercept`, `offset`) and
##    the regression coefficients `lm.coeff` into `R/vmideal_vst.R` (see
##    `vmidealVstFromMagn` and `vmidealVstToPsi`). The package keeps these
##    values inline for speed and reproducibility.
## 3. Validate the new values with `devtools::load_all()` followed by the
##    relevant unit tests (for example `test_vmideal_vst.R`) before committing
##    the changes.
##
## Output locations in the package
## - `param.df` → embedded table inside `R/vmideal_vst.R`.
## - `lm.coeff` → polynomial used by `vmidealVstToPsi` in
##   `R/vmideal_vst.R`.
## - Intermediate diagnostics (`data.t`, `lm.res`) are for exploratory use and
##   are not saved automatically.
##
## Dependencies
## - Requires 'vismeteor' for dvmideal and helper polynomial functions,
##   and base R.
##
## Notes on execution
## - if (FALSE) sections contain time‑consuming optimization and are not run by
##   default. Coefficients below are precomputed and used as is.
## - myVmidealVstToPsi expects polynomial coefficients as produced by the
##   regression at the end of this script; keep the valid x‑domain checks in
##   mind when reusing the mapping.
##
# Approximation of the Variance-stabilizing transformation
library(vismeteor)

m <- seq(6, -200, -1)

#' Variance‑stabilizing transform approximation for vmideal
#'
#' Builds t(m; limmag) from a cubic spline in x = limmag − m with anchor
#' points sx = {1, 2, 4, 6, 8, 10}. The spline values at sx depend on the
#' fractional component (offset) of limmag via linear interpolation of
#' param.df columns '1'..'6'. An offset‑specific intercept is subtracted so
#' that the transformation is centered consistently across offsets.
#'
#' Arguments
#' - m: numeric vector of magnitudes.
#' - limmag: scalar or vector of limiting magnitudes corresponding to m.
#' - param.df: data.frame with columns 'offset', '1'..'6', and 'intercept'.
#'
#' Returns
#' - Numeric vector t with the same length as m.
myVmidealVstFromMagn <- function(m, limmag, param.df) {
    offset <- limmag - round(limmag)
    if (1L == length(limmag)) {
        limmag <- rep(limmag, length(m))
        offset <- rep(offset, length(m))
    }

    sx <- c(1.0, 2.0, 4.0, 6.0, 8.0, 10.0)
    p0.fun <- approxfun(param.df$offset, param.df[['intercept']])
    p1.fun <- approxfun(param.df$offset, param.df[['1']])
    p2.fun <- approxfun(param.df$offset, param.df[['2']])
    p3.fun <- approxfun(param.df$offset, param.df[['3']])
    p4.fun <- approxfun(param.df$offset, param.df[['4']])
    p5.fun <- approxfun(param.df$offset, param.df[['5']])
    p6.fun <- approxfun(param.df$offset, param.df[['6']])

    # x is the excess relative to the limiting magnitude
    arg.data <- data.frame(
        x = limmag - m,
        offset = offset
    )

    data.f <- as.factor(arg.data$offset)
    data.s <- split(arg.data, data.f)
    y <- lapply(data.s, function(data) {
        x <- data$x
        offset <- data$offset[1]
        sy <- c(
            p1.fun(offset),
            p2.fun(offset),
            p3.fun(offset),
            p4.fun(offset),
            p5.fun(offset),
            p6.fun(offset)
        )
        f.spline <- splinefun(sx, sy)
        y <- rep(NA, length(x))

        # Linear extrapolation below the first anchor using the local slope
        idx <- x < sx[1]
        if (any(idx)) {
            y[idx] <- f.spline(sx[1], deriv = 1) * (x[idx] - sx[1]) + sy[1]
        }

        # Linear extrapolation above the last anchor
        idx <- x > sx[length(sx)]
        if (any(idx)) {
            y[idx] <- f.spline(sx[length(sx)], deriv = 1) * (x[idx] - sx[length(sx)]) + sy[length(sy)]
        }

        # Spline interpolation within the anchor interval
        idx <- x >= sx[1] & x <= sx[length(sx)]
        if (any(idx)) {
            y[idx] <- f.spline(x[idx])
        }

        y - p0.fun(offset)
    })

    unsplit(y, data.f)
}

# NOTE: Mirror any edits of this table into `R/vmideal_vst.R::vmidealVstFromMagn`.
# The runtime code uses a literal copy of these coefficients.
param.df <- structure(list(limmag = c(5.5, 5.52, 5.55, 5.6, 5.7, 5.8, 5.9,
    6, 6.1, 6.2, 6.3, 6.4, 6.45, 6.48, 6.5), `1` = c(0.379578706193683,
    0.380865978506213, 0.383646581863156, 0.386594955357005, 0.384516871380943,
    0.378119364136071, 0.39783016804395, 0.373393379287508, 0.366393502192895,
    0.345593857555874, 0.367086720870886, 0.359195376100079, 0.347082537633701,
    0.312420722892609, 0.312232925285191), `2` = c(0.758955184180252,
    0.757434042449928, 0.756545241530695, 0.754597324053272, 0.746096530130275,
    0.736273688899513, 0.754472950700777, 0.729839979801878, 0.724090251762525,
    0.706272163875503, 0.732424006749156, 0.730856947686901, 0.722541149425633,
    0.690367122987635, 0.691926213129086), `3` = c(1.89943177726206,
    1.89844069948433, 1.89825430762927, 1.89725244849407, 1.88997593310439,
    1.88072563899649, 1.89899661183851, 1.87396958675163, 1.86749053248022,
    1.84882849631795, 1.87410789172307, 1.87172629242088, 1.86305301012132,
    1.83068592572223, 1.83212684285442), `4` = c(3.27014476551314,
    3.26892787903128, 3.26844009590108, 3.26702893292816, 3.25921320867758,
    3.24969559483436, 3.26789784387066, 3.24296948683872, 3.23673983810267,
    3.21843295052407, 3.24413734533724, 3.24221462073629, 3.23377186059092,
    3.20154167652271, 3.20307252532895), `5` = c(4.50424034207103,
    4.50305721713432, 4.50261471059499, 4.5012654394704, 4.49353145072154,
    4.48405231361262, 4.50225758259143, 4.47729738300568, 4.47099710887246,
    4.45259464315401, 4.47819774376931, 4.4761817895577, 4.46769920746914,
    4.43544755919314, 4.43696568503357), `6` = c(5.65036712145262,
    5.64927790140628, 5.64896064069461, 5.6477824229466, 5.64028041416121,
    5.63092781192841, 5.64918249953455, 5.62420382559607, 5.61781703615347,
    5.5992741104315, 5.62469988859743, 5.62248449312786, 5.61389842883615,
    5.58158409117071, 5.58306055049816), offset = c(-0.5, -0.48,
    -0.45, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.45, 0.48,
    0.5), intercept = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0)), row.names = c("5.5", "5.52", "5.55", "5.6", "5.7", "5.8",
    "5.9", "6", "6.1", "6.2", "6.3", "6.4", "6.45", "6.48", "6.5"
    ), class = "data.frame")

if (FALSE) {
    # Estimate the anchor values (columns '1'..'6') per limmag by minimizing
    # deviation of t-variance from a target (here 1.0) across psi.
    param.df <- with(new.env(), {
        limmag <- c(5.5, 5.52, 5.55, seq(5.6, 6.4, 0.1), 6.45, 6.48, 6.5)
        psi <- c(-100, seq(-4, 9, 0.25), 100)

        f.model <- function(param.df, psi, limmag) {
            model <- expand.grid(psi = psi, limmag = limmag)
            model <- do.call(
                rbind.data.frame,
                mapply(function(psi, limmag) {
                    p <- dvmideal(m, limmag, psi)
                    t <- myVmidealVstFromMagn(m, limmag, param.df)
                    t.mean <- sum(p * t)
                    t.var <- sum(p * (t - t.mean)^2)

                    list(
                        t.mean = t.mean,
                        t.var = t.var
                    )
                }, model$psi, model$limmag , SIMPLIFY = FALSE)
            )
        }

        param.list <- list()
        for (l in limmag) {
            print(paste('limmag:', l))
            params <- as.vector(unlist(param.df[as.character(l), as.character(seq(1, 6))]))
            opt.res <- optim(params, fn = function(params, t.var.mean) {
                param.row <- as.data.frame(t(params))
                param.df <- rbind.data.frame(param.row, param.row)
                colnames(param.df) <- seq(1, length(params))
                param.df$intercept <- 0.0
                param.df$offset <- l - round(l)
                param.df$offset[1] <- param.df$offset[1] - 0.01
                param.df$offset[2] <- param.df$offset[2] + 0.01

                df <- f.model(param.df, psi, l)
                sum((df$t.var - t.var.mean)^2)/nrow(df)
            }, t.var.mean = 1.0)
            params  <- c(l, opt.res[["par"]])
            param.list[[as.character(l)]] <- params
        }
        param.df <- do.call(
            rbind.data.frame,
            as.list(param.list)
        )
        colnames(param.df) <- c('limmag', as.character(seq(1, 6)))
        param.df <- param.df[order(param.df$limmag),]
        row.names(param.df) <- as.character(param.df$limmag)
        param.df
    })
    param.df$offset <- round(param.df$limmag - round(param.df$limmag), 2)
    param.df$intercept <- 0.0
    dput(param.df)
}

#' Generate E[t] and Var[t] over grids of psi and limmag for vmideal
#'
#' Arguments
#' - param.df: parameter frame with columns 'offset', '1'..'6', 'intercept'.
#' - limmag: vector of limiting magnitudes.
#' - psi: vector of vmideal shape parameters.
#'
#' Returns
#' - data.frame with columns psi, limmag, t.mean, t.var.
data.t.fun <- function(param.df, limmag, psi) {
    model <- expand.grid(psi = psi, limmag = limmag)
    model <- do.call(
        rbind.data.frame,
        mapply(function(psi, limmag) {
            p <- dvmideal(m, limmag, psi)
            t <- myVmidealVstFromMagn(m, limmag, param.df)
            t.mean <- sum(p * t)
            t.var <- sum(p * (t - t.mean)^2)

            list(
                psi = psi,
                limmag = limmag,
                t.mean = t.mean,
                t.var = t.var
            )
        }, model$psi, model$limmag , SIMPLIFY = FALSE)
    )
}

# Intercept calibration (example MSE ~ 8.696629e-05 during tuning)
# Intercepts are offset-specific shifts; keep them synchronised with the
# constant column in `vmidealVstFromMagn`.
param.df$intercept <- c(0.0318278656993037, 0.0299076409185304, 0.0284655814974199,
    0.025667223534113, 0.0158031837415663, 0.00504408335887503, 0.0226347191720595,
    -0.00235692008208782, -0.0082983458149825, -0.0258911028895704,
    0.000925136964883697, 0.000504234344480941, -0.00707058067955687,
    -0.0387371640018533, -0.0367814650209238
)
limmag <- c(5.5, 5.52, 5.55, seq(5.6, 6.4, 0.1), 6.45, 6.48, 6.5)
psi <- c(-100, seq(-10, 9, 0.25))

if (FALSE) {
    # Fine‑tune the intercept per offset by minimizing squared error in the
    # regression psi - limmag ~ poly(log(E[t]), 8).
    opt.res <- optim(param.df$intercept, fn = function(intercept) {
        param.df$intercept <- intercept
        data.t <- data.t.fun(param.df, limmag, psi)
        data.t$x <- log(data.t$t.mean - min(data.t.fun(param.df, limmag, 100)$t.mean) + 1e-05)
        data.t$y <- data.t$psi - data.t$limmag
        lm.res <- lm(y ~ poly(x, 8, raw=TRUE), data = data.t)
        data.t$y.predicted <- predict(lm.res)
        print(sum((data.t$y - data.t$y.predicted)^2)/nrow(data.t))
    }, control = list(maxit = 250))
    dput(opt.res[["par"]])
    param.df$intercept <- opt.res[["par"]]
}
# Align intercepts by adding the minimal E[t] at high psi so that log(E[t])
# is well‑behaved for the subsequent regression/inversion mapping.
param.df$intercept <- param.df$intercept + min(data.t.fun(param.df, limmag, 1000)$t.mean)

#' Map x = E[t] to psi via polynomial in log(x)
#'
#' Arguments
#' - x: positive numeric vector (E[t]). Values outside [0.016, 8.22] are set NA.
#' - limmag: corresponding limiting magnitudes.
#' - poly.coef: numeric coefficients of a polynomial in log(x). The k-th entry
#'   corresponds to exponent (k-1).
#' - deriv: if TRUE, return derivative d psi / d x (scaled by 1/x) using the
#'   derivative polynomial.
#'
#' Returns
#' - psi (same shape as x) or its derivative if deriv = TRUE.
myVmidealVstToPsi <- function(x, limmag, poly.coef, deriv = FALSE) {
    names(poly.coef) <- seq_along(poly.coef) - 1 # exponents
    # x min 0.016 (psi approx 9 at limiting maginitde of 5.5)
    # x max 8.22(psi approx -10 at limiting maginitde of 6.5)
    x[x < 0.016 | x > 8.22] <- NA
    x <- log(x)

    if(deriv) {
        poly.coef1 <- vismeteor:::f.polynomial.coef(poly.coef, deriv.degree = 1L)
        vismeteor:::f.polynomial(x, poly.coef1) / x
    } else {
        limmag + vismeteor:::f.polynomial(x, poly.coef)
    }
}

# evaluate vmideal expectations with the current spline/intercept set.
data.t <- data.t.fun(param.df, limmag, psi)
data.t$psi.delta <- data.t$psi - data.t$limmag
data.t$x <- log(data.t$t.mean)
data.t$y <- data.t$psi - data.t$limmag
lm.res <- lm(y ~ poly(x, 8, raw=TRUE), data = data.t)
lm.coeff <- coefficients(lm.res)
names(lm.coeff) <- NULL
dput(lm.coeff)  # Copy this numeric vector into `R/vmideal_vst.R::vmidealVstToPsi`.

data.t$y.predicted <- predict(lm.res)
data.t$y.diff <- data.t$y.predicted - data.t$y
data.t$psi.est <- myVmidealVstToPsi(data.t$t.mean, data.t$limmag, lm.coeff)
#data.t$psi.est <- vmidealVstToPsi(data.t$t.mean, data.t$limmag)
data.t$psi.est.diff <- round(data.t$psi.est - data.t$psi, 3)
#View(data.t)
