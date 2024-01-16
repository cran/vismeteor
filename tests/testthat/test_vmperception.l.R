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

    data.laplace <- with(new.env(), {
        limmag <- seq(5.6, 6.5, 0.2)
        m <- seq(-200, 6, 1)
        q <- log(seq(1.2, 4.0, 0.1))

        data.combined <- expand.grid(limmag = limmag, m = m)
        do.call(
            rbind.data.frame,
            sapply(q, function(q) {
                f0 <- function(m) {
                    q * exp(-q*m) * p.fun(m)
                }
                qL <- stats::integrate(f0, -0.5, Inf)$value

                m.mean <- with(data.combined, {
                    sum((limmag - m) * dvmgeom(m, limmag, exp(q)))
                })/length(limmag)

                m.var <- with(data.combined, {
                    x <- (limmag - m - m.mean)^2
                    sum(x * dvmgeom(m, limmag, exp(q)))
                })/length(limmag)

                list(
                    q = q,
                    qL = qL,
                    m.mean = m.mean,
                    m.var = m.var
                )
            }, simplify = FALSE)
        )
    })

    # test LaplaceTrans(perception)
    L <- log(vmperception.l(data.laplace$q))
    expect_true(all(abs(L - log(data.laplace$qL/data.laplace$q)) < 0.021))

    # test q * LaplaceTrans(perception)
    qL <- data.laplace$q * vmperception.l(data.laplace$q)
    expect_true(all(abs(qL - data.laplace$qL) < 0.00189))

    # test mean of (limmag - m)
    m.mean.log <- log(-vmperception.l(data.laplace$q, deriv.degree = 1L)/vmperception.l(data.laplace$q))
    expect_true(all(
        abs(m.mean.log - log(data.laplace$m.mean)
    ) < 0.013))

    # test variance of (limmag - m)
    m.var.log <- log(
        vmperception.l(data.laplace$q, deriv.degree = 2L)/vmperception.l(data.laplace$q) -
            (-vmperception.l(data.laplace$q, deriv.degree = 1L)/vmperception.l(data.laplace$q))^2
    )
    expect_true(all(
        abs(m.var.log - log(data.laplace$m.var)
    ) < 0.03))

    # test first derivative
    f <- function(s) {
        vmperception.l(s, deriv.degree = 1L)
    }
    res <- vmperception.l(4.0) - vmperception.l(0.1)
    expect_true(abs(res - stats::integrate(f, 0.1, 4.0)$value) < 1e-08)

    # test second derivative
    f <- function(s) {
        vmperception.l(s, deriv.degree = 2L)
    }
    res <- vmperception.l(4.0, deriv.degree = 1L) - vmperception.l(0.1, deriv.degree = 1L)
    expect_true(abs(res - stats::integrate(f, 0.1, 4.0)$value) < 1e-08)
    expect_true(abs(res - stats::integrate(f, 0.1, 4.0)$value) < 1e-08)
})
