test_that("dvmideal", {
    lm <- 6.5
    psi <- 4.0

    # from documentation
    dp.fun <- function(m, psi) {
        r <- 10^0.4
        1.5 * log(r) * sqrt(r^(3 * psi + 2 * m)/((r^psi + r^m)^5))
    }

    for (current.psi in c(-10.0, psi, 5)) {
        m <- seq(6, current.psi - 25, -1)
        expected_p <- sapply(m, function(m) {
            stats::integrate(function(m) dp.fun(m, current.psi), m - 0.5, m + 0.5)$value *
                vismeteor::vmperception(lm - m)
        })
        expected_p <- expected_p/sum(expected_p)

        # density of meteor magnitudes
        p <- vismeteor::dvmideal(m, lm, current.psi)
        expect_type(p, 'double')
        expect_length(p, length(m))
        expect_false(any(abs(log(p) - log(expected_p)) > 0.03), label = paste('psi =', current.psi))

        # log density of meteor magnitudes
        p <- vismeteor::dvmideal(m, lm, current.psi, log = TRUE)
        expect_type(p, 'double')
        expect_length(p, length(m))
        expect_false(any(abs(p - log(expected_p)) > 0.03), label = paste('psi =', current.psi))

        # Is the density sum of meteor magnitudes equal to 1.0 ?
        p <- sum(vismeteor::dvmideal(m, lm, current.psi))
        expect_equal(round(p, 6), 1.0, label = paste('psi =', current.psi))
    }

    p <- vismeteor::dvmideal(7, lm, psi)
    expect_type(p, 'double')
    expect_length(p, 1)
    expect_equal(p, 0.0)

    # log density of meteor magnitudes
    p <- vismeteor::dvmideal(7, lm, psi, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, 1)
    expect_equal(p, -Inf)

    p <- vismeteor::dvmideal(c(6, 6, 6), lm = c(5.4, 5.5, 5.6), psi)
    expect_type(p, 'double')
    expect_length(p, 3)
    expect_equal(p[1], 0.0)
    expect_equal(p[2], 0.0)
    expect_gt(p[3], 0.0)

    # density of infinite meteor magnitudes
    p <- vismeteor::dvmideal(c(-Inf, Inf), lm, psi)
    expect_type(p, 'double')
    expect_length(p, 2)
    expect_equal(p[1], 0.0)
    expect_equal(p[2], 0.0)

    # log density of infinite meteor magnitudes
    p <- vismeteor::dvmideal(c(-Inf, Inf), lm, psi, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, 2)
    expect_equal(p[1], -Inf)
    expect_equal(p[2], -Inf)

    # test order of densities of meteor magnitudes
    lm <- seq(2.5, 6.5, 0.1)
    p <- vismeteor::dvmideal(rep(3, length(lm)), lm, psi)
    expect_equal(p[order(p)], p)

    # test maximum likelihood estimation (MLE) of psi
    lm <- 5.8
    m <- as.integer(seq(-20, 6))
    p <- vismeteor::dvmideal(m, lm, psi)
    ll <- function(psi) { # log likelihood function
        -sum(p * vismeteor::dvmideal(m, lm, psi, log = TRUE))
    }
    est <- optim(1, ll, method='Brent', lower=0, upper=7, hessian=TRUE)
    expect_equal(round(est$par, 6), psi)

    # density of meteor magnitudes equals geometric distribution
    lm <- 6.3
    m <- as.integer(seq(-20, 6))

    p <- vismeteor::dvmideal(m, lm, 16.25, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, length(m))
    diff.p <- abs(p - dvmgeom(m, lm, 10^0.4, log = TRUE))
    expect_true(any(diff.p > 0.0))
    expect_false(any(diff.p > 1e-08))

    p <- vismeteor::dvmideal(m, lm, 16.35, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_equal(p, dvmgeom(m, lm, 10^0.4, log = TRUE))
})
