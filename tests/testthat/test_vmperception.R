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
    pdiff <- (p - data0$p)/data0$p
    expect_lt(mean(abs(pdiff)), 0.057)
    expect_true(all(abs(pdiff) < 0.14))

    # Checks whether the approximation formula changes the r-value
    model <- with(new.env(), {
        r <- seq(1.6, 3.0, 0.1)
        limmag <- seq(5.6, 6.4, 0.2)
        m <- seq(-100, 6, 1)

        df <- expand.grid(r=r, limmag=limmag)
        do.call(
            rbind.data.frame,
            mapply(function(r, limmag) {
                p <- dvmgeom(m, limmag, r, perception.fun = p.fun)

                # maximum likelihood estimation (MLE) of r
                llr <- function(r) {
                    -sum(p * dvmgeom(m, limmag, r, log=TRUE))
                }
                r.est <- optim(r, llr, method='Brent', lower=1.1, upper=5, hessian=FALSE)$par

                list(
                    r = r,
                    limmag = limmag,
                    r.est = r.est
                )
            }, df$r, df$limmag, SIMPLIFY = FALSE)
        )
    })

    rdiff <- (model$r - model$r.est)
    expect_lt(mean(abs(rdiff)), 0.0041)
    expect_true(all(abs(rdiff) < 0.022))

    # lower r values
    model.rlow <- subset(model, model$r < 2.7)

    rdiff <- (model.rlow$r - model.rlow$r.est)
    expect_lt(mean(abs(rdiff)), 0.0021)
    expect_true(all(abs(rdiff) < 0.009))
})
