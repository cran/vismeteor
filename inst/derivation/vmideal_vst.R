# Approximation of the Variance-stabilizing transformation
library(vismeteor)

m <- seq(6, -200, -1)

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

        y - p0.fun(offset)
    })

    unsplit(y, data.f)
}

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

# 8.696629e-05
param.df$intercept <- c(0.0318278656993037, 0.0299076409185304, 0.0284655814974199,
    0.025667223534113, 0.0158031837415663, 0.00504408335887503, 0.0226347191720595,
    -0.00235692008208782, -0.0082983458149825, -0.0258911028895704,
    0.000925136964883697, 0.000504234344480941, -0.00707058067955687,
    -0.0387371640018533, -0.0367814650209238
)
limmag <- c(5.5, 5.52, 5.55, seq(5.6, 6.4, 0.1), 6.45, 6.48, 6.5)
psi <- c(-100, seq(-10, 9, 0.25))

if (FALSE) {
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
param.df$intercept <- param.df$intercept + min(data.t.fun(param.df, limmag, 1000)$t.mean)

myVmidealVstToPsi <- function(x, limmag, poly.coef, deriv = FALSE) {
    names(poly.coef) <- seq_along(poly.coef) - 1 # exponents
    # x min 0.016 (psi approx 9 at limiting maginitde of 5.5)
    # x max 8.22(psi approx -10 at limiting maginitde of 6.5)
    x[x < 0.016 | x > 8.22] <- NA
    x <- log(x)

    if(deriv) {
        poly.coef1 <- f.polynomial.coef(poly.coef, deriv.degree = 1L)
        f.polynomial(x, poly.coef1) / x
    } else {
        limmag + f.polynomial(x, poly.coef)
    }
}

data.t <- data.t.fun(param.df, limmag, psi)
data.t$psi.delta <- data.t$psi - data.t$limmag
data.t$x <- log(data.t$t.mean)
data.t$y <- data.t$psi - data.t$limmag
lm.res <- lm(y ~ poly(x, 8, raw=TRUE), data = data.t)
lm.coeff <- coefficients(lm.res)
names(lm.coeff) <- NULL
dput(lm.coeff)

data.t$y.predicted <- predict(lm.res)
data.t$y.diff <- data.t$y.predicted - data.t$y
data.t$psi.est <- myVmidealVstToPsi(data.t$t.mean, data.t$limmag, lm.coeff)
#data.t$psi.est <- vmidealVstToPsi(data.t$t.mean, data.t$limmag)
data.t$psi.est.diff <- round(data.t$psi.est - data.t$psi, 3)
#View(data.t)
