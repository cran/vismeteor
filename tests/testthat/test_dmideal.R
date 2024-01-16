test_that("dmideal", {
    psi <- 4.0

    # from documentation
    dp.fun <- function(m, psi) {
        r <- 10^0.4
        1.5 * log(r) * sqrt(r^(3 * psi + 2 * m)/((r^psi + r^m)^5))
    }

    expected_p <- dp.fun(3, psi)

    # density of meteor magnitudes
    p <- vismeteor::dmideal(c(-Inf, Inf, 3), psi)
    expect_type(p, 'double')
    expect_length(p, 3)
    expect_equal(p[1], 0.0)
    expect_equal(p[2], 0.0)
    expect_equal(p[3], expected_p)

    # density of meteor magnitudes (shifting property)
    p <- vismeteor::dmideal(c(2, 3, 4), c(psi-1, psi, psi+1))
    expect_type(p, 'double')
    expect_length(p, 3)
    expect_equal(p[1], expected_p)
    expect_equal(p[2], expected_p)
    expect_equal(p[3], expected_p)

    # log density of meteor magnitudes
    p <- vismeteor::dmideal(c(-Inf, Inf, 3), psi, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, 3)
    expect_equal(p[1], -Inf)
    expect_equal(p[2], -Inf)
    expect_equal(p[3], log(expected_p))

    # test maximum likelihood estimation (MLE) of psi
    m <- seq(-20, 20, 0.001)
    p <- vismeteor::dmideal(m, psi)
    ll <- function(psi) { # log likelihood function
        -sum(p * vismeteor::dmideal(m, psi, log = TRUE))
    }
    est <- optim(2, ll, method='Brent', lower=0, upper=7, hessian=TRUE)
    expect_equal(round(est$par, 6), psi)
})
