test_that("vmperception", {
    # test approximation of `KOSREN90` perceptions published in
    # J. Rendtel, Handbook for meteor observers 2022 edition, p. 149, (c) International Meteor Organization
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

    # ignore lower magnitudes
    data0 <- subset(data, data$m >= 0.29)

    # check perceptions
    p <- vmperception(data0$m)
    expect_true(all(p < 1.0))
    expect_true(all(p > 0.0))
    rdiff <- (p - data0$p)/data0$p
    expect_lt(mean(abs(rdiff)), 0.06)
    expect_true(all(abs(rdiff) < 0.121))

    # test first derivative
    f <- function(m) {
        vmperception(m, deriv.degree = 1L)
    }
    p <- vmperception(4.0) - vmperception(1.0)
    expect_true(abs(p - stats::integrate(f, 1.0, 4.0)$value) < 1e-10)

    # test second derivative
    f <- function(m) {
        vmperception(m, deriv.degree = 2L)
    }
    p <- vmperception(4.0, deriv.degree = 1L) - vmperception(1.0, deriv.degree = 1L)
    expect_true(abs(p - stats::integrate(f, 1.0, 4.0)$value) < 1e-10)

    # Checks whether the approximation formula changes the r-value
    data.q <- with(new.env(), {
        r <- seq(1.3, 3.6, 0.1)
        limmag <- seq(5.6, 6.5, 0.2)
        m <- seq(-100, 6, 1)

        do.call(
            rbind.data.frame,
            sapply(r, function(r) {
                obs <- do.call(
                    rbind.data.frame,
                    sapply(limmag, function(limmag) {
                        p <- dgeom(6 - m, 1.0 - 1.0/r) * p.fun(limmag - m)
                        p <- p/sum(p)

                        list(
                            limmag = rep(limmag, length(m)),
                            m = m,
                            p = p,
                            q = p * vmperception(limmag - m, TRUE)/vmperception(limmag - m)
                        )
                    }, simplify = FALSE)
                )

                r.approx <- with(obs,{
                    llr <- function(r) {
                        -sum(p * dvmgeom(m, limmag, r, log=TRUE))
                    }

                    optim(2, llr, method='Brent', lower=1.1, upper=5, hessian=FALSE)
                })

                list(
                    q = log(r),
                    qr = log(r.approx$par),
                    qq = sum(obs$q)/length(limmag)
                )
            }, simplify = FALSE)
        )
    })

    qdiff <- (data.q$q - data.q$qr)
    expect_lt(mean(abs(qdiff)), 0.0007)
    expect_true(all(abs(qdiff) < 0.0025))

    qdiff <- (data.q$q - data.q$qq)
    expect_lt(mean(abs(qdiff)), 0.001)
    expect_true(all(abs(qdiff) < 0.0016))
})

# Approximation of perception probabilities
if(FALSE) {
    #library(vismeteor)
    #library(ggplot2)

    # exact like vmperception(), but with polynomial coefficients as argument
    perception.fun <- function(poly.coef, m, deriv.degree = 0L) {
        names(poly.coef) <- seq(along = poly.coef) # exponents

        m <- m + 0.5
        p <- rep(0.0, length(m))
        idx <- m > .Machine$double.eps
        if (any(idx)) {
            f0 <- vismeteor:::f.polynomial(m[idx], poly.coef)
            if (0L == deriv.degree) {
                p[idx] <- 1.0 - exp(-f0)
            } else if (1L == deriv.degree) {
                poly.coef1 <- vismeteor:::f.polynomial.coef(poly.coef, deriv.degree = 1L)
                f1 <- vismeteor:::f.polynomial(m[idx], poly.coef1)
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

    if (FALSE) {
        coef.model <- with(new.env(), {

            optim.fun <- function(params) {
                q <- seq(0.2, 1.4, 0.02)
                r <- exp(q)

                limmag <- seq(5.6, 6.5, 0.2)
                m <- seq(-100, 6, 1)
                poly.coef <- exp(params)

                vmperception.local <- function(m, deriv.degree = 0L) {
                    perception.fun(poly.coef, m, deriv.degree)
                }

                result <- do.call(
                    rbind.data.frame,
                    sapply(r, function(r) {
                        obs <- do.call(
                            rbind.data.frame,
                            sapply(limmag, function(limmag) {
                                p <- dvmgeom(m, limmag, r, perception.fun = p.fun)

                                list(
                                    limmag = rep(limmag, length(m)),
                                    m = m,
                                    p = p,
                                    q = p * vmperception.local(limmag - m, TRUE)/vmperception.local(limmag - m)
                                )
                            }, simplify = FALSE)
                        )

                        r.approx <- with(obs,{
                            # log likelihood function
                            llr <- function(r) {
                                -sum(p * dvmgeom(m, limmag, r, log=TRUE, perception.fun = vmperception.local))
                            }

                            # maximum likelihood estimation (MLE) of r
                            optim(r, llr, method='Brent', lower=1.1, upper=5, hessian=FALSE)
                        })

                        list(
                            q = log(r),
                            qr = log(r.approx$par),
                            qq = sum(obs$q)/length(limmag)
                        )
                    }, simplify = FALSE)
                )

                print(with(result, {
                    sum((q - qr)^2 + (q - qq)^2)
                }))
            }

            optim(log(c(0.00028, 0.0081, 1e-07, 0.0012)), optim.fun, control = list('reltol' = 1e-06))$par
        })

        coef.model <- exp(coef.model)
    } else {
        # round
        # 0.00303959663899688, 0.00558016073459467, 4.98709880344163e-09, 0.00140590886426205
        coef.model <- c(0.003, 0.0056, 0, 0.0014)
    }
    names(coef.model) <- seq(along = coef.model) # exponents
    print(paste(c('coef model:', paste(coef.model, collapse = ', '))))

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
                    limits = c(0.0003, 12),
                    breaks = c(0.0005, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10)
                )
            print(p)
        })
    }

    if (TRUE) {
        with(new.env(), {
            r <- seq(1.2, 4.0, 0.1)
            limmag <- seq(5.6, 6.5, 0.2)
            m <- seq(-100, 6, 1)
            vmperception.local <- function(m, deriv.degree = 0L) {
                perception.fun(coef.model, m, deriv.degree)
            }

            result <- do.call(
                rbind.data.frame,
                sapply(r, function(r) {
                    obs <- do.call(
                        rbind.data.frame,
                        sapply(limmag, function(limmag) {
                            p <- dvmgeom(m, limmag, r, perception.fun = p.fun)
                            #p <- dvmgeom(m, limmag, r, perception.fun = vmperception.local)

                            list(
                                limmag = rep(limmag, length(m)),
                                m = m,
                                p = p,
                                q = p * vmperception.local(limmag - m, TRUE)/vmperception.local(limmag - m)
                            )
                        }, simplify = FALSE)
                    )

                    r.approx <- with(obs,{
                        # log likelihood function
                        llr <- function(r) {
                            -sum(p * dvmgeom(m, limmag, r, log=TRUE, perception.fun = vmperception.local))
                        }

                        # maximum likelihood estimation (MLE) of r
                        optim(r, llr, method='Brent', lower=1.1, upper=5, hessian=FALSE)
                    })

                    list(
                        r = r,
                        r.est = r.approx$par,
                        q.exp = exp(sum(obs$q)/length(limmag))
                    )
                }, simplify = FALSE)
            )

            print(result)
        })
    }
}
