# Approximation of the Variance-stabilizing transformation
library(vismeteor)
library(ggplot2)

r <- seq(1.4, 3.5, 0.1)
m <- seq(-200, 6, 1)

# dput(param.df)
param.df <- structure(list(limmag = c(5.5, 5.52, 5.55, 5.6, 5.7, 5.8, 5.9,
    6, 6.1, 6.2, 6.3, 6.4, 6.45, 6.48, 6.5), b = c(2.176, 2.077,
    2.004, 1.945, 1.897, 1.885, 1.889, 1.906, 1.933, 1.97, 2.02,
    2.087, 2.128, 2.156, 2.176), c = c(-0.317, -0.278, -0.249, -0.225,
    -0.206, -0.2, -0.201, -0.207, -0.217, -0.231, -0.252, -0.279,
    -0.296, -0.308, -0.317), d = c(0.669, 0.748, 0.813, 0.872, 0.922,
    0.934, 0.928, 0.91, 0.881, 0.843, 0.796, 0.739, 0.706, 0.684,
    0.669)), row.names = c("5.5", "5.52", "5.55", "5.6", "5.7", "5.8",
    "5.9", "6", "6.1", "6.2", "6.3", "6.4", "6.45", "6.48", "6.5"
    ), class = "data.frame")

if (FALSE) {
    param.df <- with(new.env(), {
        limmag <- c(5.52, 5.55, seq(5.6, 6.4, 0.1), 6.45, 6.48, 6.5)
        myVmgeomVstFromMagn <- function(m, limmag, params) {
            a <- 0.0
            b <- params[1]
            c <- params[2]
            d <- params[3]

            x <- limmag - m
            a - exp(b + c * (x + 0.5)^d)
        }

        f.model <- function(params, limmag) {
            model <- expand.grid(r = r, limmag = limmag)
            model <- do.call(
                rbind.data.frame,
                mapply(function(r, limmag) {
                    p <- dvmgeom(m, limmag, r)
                    t <- myVmgeomVstFromMagn(m, limmag, params)
                    t.mean <- sum(p * t)
                    t.var <- sum(p * (t - t.mean)^2)

                    list(
                        r = r,
                        q = log(r),
                        limmag = limmag,
                        t.mean = t.mean,
                        t.var = t.var
                    )
                }, model$r, model$limmag , SIMPLIFY = FALSE)
            )
        }

        param.list <- list()
        for (l in limmag) {
            print(paste('limmag:', l))
            params <- as.vector(unlist(param.df[as.character(l), c('b', 'c', 'd')]))
            opt.res <- optim(params, fn = function(params, t.var.mean) {
                df <- f.model(params, l)
                sum((df$t.var - t.var.mean)^2)/nrow(df)
            }, t.var.mean = 1.0)
            params  <- c(l, opt.res[["par"]])
            param.list[[as.character(l)]] <- params
        }
        param.df <- do.call(
            rbind.data.frame,
            as.list(param.list)
        )
        colnames(param.df) <- c('limmag', 'b', 'c', 'd')

        d55 <- param.df[param.df == 6.5,]
        d55$limmag <- 5.5
        param.df <- rbind.data.frame(d55, param.df)
        param.df$b <- round(param.df$b, 3)
        param.df$c <- round(param.df$c, 3)
        param.df$d <- round(param.df$d, 3)
        param.df <- param.df[order(param.df$limmag),]
        row.names(param.df) <- as.character(param.df$limmag)
        param.df
    })
    dput(param.df)
}

param.df$offset <- round(param.df$limmag - round(param.df$limmag), 2)

myVmgeomVstFromMagn <- function(m, limmag, param.df) {
    offset <- limmag - round(limmag)
    if (1L == length(limmag)) {
        limmag <- rep(limmag, length(m))
        offset <- rep(offset, length(m))
    }

    pa.fun <- approxfun(param.df$offset, param.df$a)
    pb.fun <- approxfun(param.df$offset, param.df$b)
    pc.fun <- approxfun(param.df$offset, param.df$c)
    pd.fun <- approxfun(param.df$offset, param.df$d)

    a <- pa.fun(offset)
    b <- pb.fun(offset)
    c <- pc.fun(offset)
    d <- pd.fun(offset)

    x <- limmag - m
    a - exp(b + c * (x + 0.5)^d)
}

limmag <- c(5.5, 5.52, 5.55, seq(5.6, 6.4, 0.1), 6.45, 6.48, 6.5)
data.t.fun <- function(param.df) {
    model <- expand.grid(r = r, limmag = limmag)
    model <- do.call(
        rbind.data.frame,
        mapply(function(r, limmag) {
            p <- dvmgeom(m, limmag, r)

            m.mean <- sum(p * (limmag - m))
            m.var <- sum(p * (limmag - m - m.mean)^2)

            t <- myVmgeomVstFromMagn(m, limmag, param.df)
            t.mean <- sum(p * t)
            t.var <- sum(p * (t - t.mean)^2)

            list(
                r = r,
                q = log(r),
                limmag = limmag,
                m.mean = m.mean,
                m.var = m.var,
                t.mean = t.mean,
                t.var = t.var
            )
        }, model$r, model$limmag , SIMPLIFY = FALSE)
    )
}
data.t <- data.t.fun(param.df)

param.df$a <- c(8.83660725413019, 8.52231677831387, 8.31587431846188, 8.16602303316079,
    8.05750794169318, 8.0501384519307, 8.07216870699935, 8.1189428855403,
    8.19069184757109, 8.29154614666007, 8.40974020292365, 8.58969657552978,
    8.70398711417616, 8.78298253370545, 8.83710213031958)

if (FALSE) {
    opt.res <- optim(param.df$a, fn = function(a) {
        param.df$a <- a
        data.t <- data.t.fun(param.df)
        param.df$a <- round(param.df$a - myVmgeomVstFromMagn(6, 5.5, param.df), 6) + 1e-06
        data.t <- data.t.fun(param.df)
        lm.res <- lm(q ~ poly(t.mean, 4, raw = TRUE), data = data.t)
        data.t$q.predicted <- predict(lm.res)
        print(sum((data.t$q - data.t$q.predicted)^2)/nrow(data.t))
    })
    dput(opt.res[["par"]])
    param.df$a <- opt.res[["par"]]
}

param.df$a <- round(param.df$a - myVmgeomVstFromMagn(6, 5.5, param.df), 6) + 1e-06
print(paste(c('should be zero:', myVmgeomVstFromMagn(6, 5.5, param.df))))
dput(param.df)
data.t <- data.t.fun(param.df)

lm.res.t <- lm(q ~ poly(t.mean, 4, raw = TRUE), data = data.t)
lm.res.t.coef <- coef(lm.res.t)
print(lm.res.t.coef)
data.t$pred <- exp(predict(lm.res.t))
data.t$pred.diff <- data.t$r - data.t$pred

if (TRUE) {
    with(new.env(), {
        m <- seq(7.0, -4.5, -0.1)
        limmag <- seq(5.6, 6.5, 0.1)

        p <- ggplot(data.t) + theme_bw()

        for (l in limmag) {
            plot.data <- expand.grid(limmag=limmag, m=m)
            plot.data$dm <- plot.data$limmag - plot.data$m
            plot.data <- plot.data[order(plot.data$dm),]
            plot.data$x <-  myVmgeomVstFromMagn(m = plot.data$m, limmag = plot.data$limmag, param.df = param.df)
            plot.data <- subset(plot.data, plot.data$x > 0)
            p <- p + geom_line(aes(x = dm, y = x), color='blue', data = subset(plot.data, l == plot.data$limmag))
        }

        p <- p +
            scale_x_continuous(
                name = "m",
                breaks = seq(0, 11, 1.0)
            ) +
            scale_y_continuous(
                name = "t",
                breaks = seq(0, 7, 1)
            )
        print(p)
    })
}

if (TRUE) {
    with(new.env(), {
        data.t <- data.t[order(data.t$q, data.t$limmag),]
        p <- ggplot(data.t) + theme_bw()

        for (l in limmag) {
            plot.data <- subset(data.t, data.t$limmag == l)
            p <- p + geom_line(aes(x = r, y = t.mean), color='blue', data = plot.data)
        }

        p <- p +
            scale_x_continuous(
                name = "r",
                breaks = seq(1.4, 3.3, 0.1)
            ) +
            scale_y_continuous(
                name = "t.mean",
                breaks = seq(4.2, 5.9, 0.1)
            )
        print(p)
    })
}
