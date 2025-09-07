# Approximation of perception probabilities
library(vismeteor)
library(ggplot2)

# exact like vmperception(), but with polynomial coefficients as argument
perception.fun <- function(poly.coef, m, deriv.degree = 0L) {
    names(poly.coef) <- seq(along = poly.coef) # exponents

    m <- m + 0.5
    p <- rep(0.0, length(m))
    idx <- m > .Machine$double.eps
    if (any(idx)) {
        f0 <- f.polynomial(m[idx], poly.coef)
        if (0L == deriv.degree) {
            p[idx] <- 1.0 - exp(-f0)
        } else if (1L == deriv.degree) {
            poly.coef1 <- f.polynomial.coef(poly.coef, deriv.degree = 1L)
            f1 <- f.polynomial(m[idx], poly.coef1)
            p[idx] <- exp(-f0) * f1
        } else {
            stop('Not implemented')
        }
    }

    p
}

data <- data.frame(
    m = seq(-0.4, 7.6, 0.2),
    p = c(
        0.00046, 0.0011, 0.0023, 0.0046, 0.0081, 0.0122, 0.0182,
        0.0257, 0.03444, 0.04365, 0.05495, 0.06918, 0.08511, 0.10351,
        0.13, 0.16, 0.2, 0.24, 0.29, 0.34674, 0.4, 0.46, 0.52, 0.57544,
        0.63096, 0.67764, 0.71, 0.74, 0.77, 0.79, 0.81, 0.83, 0.85,
        0.87, 0.89, 0.91, 0.93, 0.94, 0.96, 0.98, 1.0
    )
)
data <- rbind(data.frame(m=-0.5, p=0.0), data)
p.fun <- approxfun(data$m, data$p, yleft = 0.0, yright = 1.0)
limmag <- seq(5.6, 6.4, 0.2)

# coef.model <- c(0.00373590946848783, 0.00189710356299022, 0.00271083620131325, 0.000899903791081514)
coef.model <- c(0.0037, 0.0019, 0.00271, 0.0009)

# estimate parameters
if (FALSE) {
    coef.model <- with(new.env(), {
        r <- seq(1.4, 3.5, 0.1)
        psi <- seq(-4, 9, 0.5)
        m <- seq(-200, 6, 1)

        model.geom.fun <- function(poly.coef, r) {
            vmperception.local <- function(m, deriv.degree = 0L) {
                perception.fun(poly.coef, m, deriv.degree)
            }

            model.geom <- expand.grid(r=r, limmag=limmag)
            model.geom <- do.call(
                rbind.data.frame,
                mapply(function(r, limmag) {
                    p.org <- dvmgeom(m, limmag, r, perception.fun = p.fun)
                    p.est <- dvmgeom(m, limmag, r, perception.fun = vmperception.local)
                    m.mean.org <- sum(p.org * m)
                    m.mean.est <- sum(p.est * m)

                    list(
                        m.mean.org = m.mean.org,
                        m.mean.est = m.mean.est
                    )
                }, model.geom$r, model.geom$limmag, SIMPLIFY = FALSE)
            )
        }

        model.ideal.fun <- function(poly.coef, psi) {
            vmperception.local <- function(m, deriv.degree = 0L) {
                perception.fun(poly.coef, m, deriv.degree)
            }

            model.ideal <- expand.grid(psi=psi, limmag=limmag)
            model.ideal <- do.call(
                rbind.data.frame,
                mapply(function(psi, limmag) {
                    p.org <- dvmideal(m, limmag, psi, perception.fun = p.fun)
                    p.est <- dvmideal(m, limmag, psi, perception.fun = vmperception.local)
                    m.mean.org <- sum(p.org * m)
                    m.mean.est <- sum(p.est * m)

                    list(
                        m.mean.org = m.mean.org,
                        m.mean.est = m.mean.est
                    )
                }, model.ideal$psi, model.ideal$limmag, SIMPLIFY = FALSE)
            )
        }

        f <- function(params) {
            model.geom <- model.geom.fun(exp(params), r = r)
            model.ideal <- model.ideal.fun(exp(params), psi = psi)
            print(
                sum((model.geom$m.mean.org - model.geom$m.mean.est)^2) +
                    sum((model.ideal$m.mean.org - model.ideal$m.mean.est)^2)
            )
        }

        optim(log(coef.model), f, control = list('reltol' = 1e-06))$par
    })
    coef.model <- exp(coef.model)
}
names(coef.model) <- seq(along = coef.model) # exponents
print(paste(c('coef model:', paste(coef.model, collapse = ', '))))

model.geom.fun <- function(poly.coef, r) {
    m <- seq(-200, 6, 1)
    vmperception.local <- function(m, deriv.degree = 0L) {
        perception.fun(poly.coef, m, deriv.degree)
    }

    model.geom <- expand.grid(r=r, limmag=limmag)
    model.geom <- model.geom[order(model.geom$r, model.geom$limmag),]
    do.call(
        rbind.data.frame,
        mapply(function(r, limmag) {
            p.org <- dvmgeom(m, limmag, r, perception.fun = p.fun)
            p.est <- dvmgeom(m, limmag, r, perception.fun = vmperception.local)

            # maximum likelihood estimation (MLE) of r
            llr <- function(r) {
                -sum(p.org * dvmgeom(m, limmag, r, log=TRUE, perception.fun = vmperception.local))
            }
            opt.result <- optim(r, llr, method='Brent', lower=1.1, upper=5, hessian=TRUE)
            r.est <- opt.result$par
            r.var <- opt.result$hessian[1][1]
            m.mean.org <- sum(p.org * (limmag - m))
            m.mean.est <- sum(p.est * (limmag - m))

            list(
                limmag = limmag,
                r = r,
                r.est = r.est,
                r.var = r.var,
                m.mean.org = m.mean.org,
                m.mean.est = m.mean.est
            )
        }, model.geom$r, model.geom$limmag, SIMPLIFY = FALSE)
    )
}

model.ideal.fun <- function(poly.coef, psi) {
    m <- seq(-200, 6, 1)
    vmperception.local <- function(m, deriv.degree = 0L) {
        perception.fun(poly.coef, m, deriv.degree)
    }

    model.ideal <- expand.grid(psi=psi, limmag=limmag)
    model.ideal <- model.ideal[order(model.ideal$psi, model.ideal$limmag),]
    do.call(
        rbind.data.frame,
        mapply(function(psi, limmag) {
            p.org <- dvmideal(m, limmag, psi, perception.fun = p.fun)
            p.est <- dvmideal(m, limmag, psi, perception.fun = vmperception.local)

            # maximum likelihood estimation (MLE) of r
            llr <- function(psi) {
                -sum(p.org * dvmideal(m, limmag, psi, log=TRUE, perception.fun = vmperception.local))
            }
            opt.result <- optim(psi, llr, method='Brent', lower=-15, upper=15, hessian=TRUE)
            psi.est <- opt.result$par
            psi.var <- opt.result$hessian[1][1]
            m.mean.org <- sum(p.org * (limmag - m))
            m.mean.est <- sum(p.est * (limmag - m))

            list(
                limmag = limmag,
                psi = psi,
                psi.est = psi.est,
                psi.var = psi.var,
                m.mean.org = m.mean.org,
                m.mean.est = m.mean.est
            )
        }, model.ideal$psi, model.ideal$limmag, SIMPLIFY = FALSE)
    )
}

if (TRUE) {
    with(new.env(), {
        result.geom <- model.geom.fun(coef.model, r = seq(1.4, 3.5, 0.1))
        result.geom <- result.geom[order(result.geom$r, result.geom$limmag),]
        print(result.geom)
        result.ideal <- model.ideal.fun(coef.model, psi = seq(-4, 9, 0.5))
        result.ideal <- result.ideal[order(result.ideal$psi, result.ideal$limmag),]
        print(result.ideal)
    })
}

# print g
if (TRUE) {
    with(new.env(), {
        m.idx <- data$m > -0.5
        g.model <- perception.fun(coef.model, data$m[m.idx])
        plot.data <- data.frame(
            m = data$m[m.idx],
            p = data$p[m.idx],
            g.model = g.model
        )
        print(plot.data)

        p <- ggplot(plot.data) +
            theme_bw() +
            geom_point(aes(x = m, y = p), color='black') +
            geom_line(aes(x = m, y = g.model), color='red') +
            scale_x_continuous(
                name = "m",
                breaks = seq(-0.5, 9.0, 0.5)
            ) +
            scale_y_continuous(
                name = "g",
                trans = "log",
                #labels = scales::comma,
                breaks = c(0.0002, 0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0)
            )
        print(p)
    })
}

# print g'
if (TRUE) {
    with(new.env(), {
        m <- c(-0.499, seq(-0.45, 8.5, 0.05))
        g.model <- perception.fun(coef.model, m, deriv.degree = 1L)
        plot.data <- data.frame(
            m = m,
            g.model = g.model
        )

        p <- ggplot(plot.data) +
            theme_bw() +
            geom_line(aes(x = m, y = g.model), color='red') +
            scale_x_continuous(
                name = "m",
                breaks = seq(-0.5, 9.0, 0.5)
            ) +
            scale_y_continuous(
                name = "g'",
                trans = "log",
                #labels = scales::comma,
                breaks = c(0.001, 0.002,0.005, 0.01, 0.02, 0.05, 0.1, 0.2)
            )
        print(p)
    })
}

# print q
if (TRUE) {
    with(new.env(), {
        m <- seq(-0.4, 8.4, 0.05)
        q <- perception.fun(coef.model, m, deriv.degree = 1L)/perception.fun(coef.model, m)
        plot.data <- data.frame(
            m = m,
            q = q
        )

        p <- ggplot(plot.data) +
            theme_bw() +
            geom_line(aes(x = m, y = q), color='blue') +
            scale_x_continuous(
                name = "m",
                limits = c(-0.5, 8.5),
                breaks = seq(-0.5, 8.5, 0.5)
            ) +
            scale_y_continuous(
                name = "q",
                trans = "log",
                #labels = scales::comma,
                limits = c(0.002, 12),
                breaks = c(0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10)
            )
        print(p)
    })
}
