test_that("dvmgeom", {
    r <- 1.8
    lm <- 6.5

    # from documentation
    dp.fun <- function(m) {
        stats::dgeom(m, 1 - 1/r) * vismeteor::vmperception(m + 0.5)
    }
    m <- seq(0, 40)
    expected_p <- dp.fun(m)
    expected_p <- expected_p/sum(expected_p)
    m <- 6.0 - m

    # density of meteor magnitudes
    p <- vismeteor::dvmgeom(m, lm, r)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_true(all(abs(p - expected_p) < 1e-10))
    expect_true(all(abs(log(p) - log(expected_p)) < 1e-09))

    p <- vismeteor::dvmgeom(c(6, 6, 6), lm = c(5.4, 5.5, 5.6), r)
    expect_type(p, 'double')
    expect_length(p, 3)
    expect_equal(p[1], 0.0)
    expect_equal(p[2], 0.0)
    expect_gt(p[3], 0.0)

    # log density of meteor magnitudes
    p <- vismeteor::dvmgeom(m, lm, r, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_true(all(abs(p - log(expected_p)) < 1e-09))

    # density of infinite meteor magnitudes
    p <- vismeteor::dvmgeom(c(-Inf, Inf), lm, r)
    expect_type(p, 'double')
    expect_length(p, 2)
    expect_equal(p[1], 0.0)
    expect_equal(p[2], 0.0)

    # log density of infinite meteor magnitudes
    p <- vismeteor::dvmgeom(c(-Inf, Inf), lm, r, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, 2)
    expect_equal(p[1], -Inf)
    expect_equal(p[2], -Inf)

    # test order of densities of meteor magnitudes
    lm <- seq(2.5, 6.5, 0.1)
    p <- vismeteor::dvmgeom(rep(3L, length(lm)), lm, r)
    expect_equal(p[order(p)], p)

    # test densities of meteor magnitudes with different parameters
    p1 <- vismeteor::dvmgeom(3, 6.5, 1.8)
    p2 <- vismeteor::dvmgeom(4, 6.3, 2.0)
    p <- vismeteor::dvmgeom(c(3, 4), c(6.5, 6.3), c(1.8, 2.0))
    expect_equal(p, c(p1, p2))

    # Is the density sum of meteor magnitudes equal to 1.0 ?
    for (lm in seq(5.4, 6.4, 0.1)) {
        p <- sum(vismeteor::dvmgeom(as.integer(seq(6, -100, -1)), lm, r))
        expect_equal(round(p, 6), 1.0, label = paste('lm =', lm))
    }

    # test mean of meteor magnitudes
    q <- log(r)
    f <- function(m) {
        exp(-q * m) * vmperception(m)
    }
    norm <- stats::integrate(f, -0.5, Inf)$value
    f <- function(m) {
        m * exp(-q * m) * vmperception(m)/norm
    }
    expected_m.mean <- stats::integrate(f, -0.5, Inf)$value
    expect_equal(round(expected_m.mean, 2), 4.42)

    f <- function(m) {
        m * vismeteor::dvmgeom(m, lm, r)
    }
    m.mean <- sum(f(as.integer(6 - seq(0, 100, 1))))
    expect_equal(round(lm - m.mean, 2), round(expected_m.mean, 2))

    # test maximum likelihood estimation (MLE) of q of meteor magnitudes
    lm <- 5.8
    m <- as.integer(seq(6, -60, -1))
    p <- vismeteor::dvmgeom(m, lm, r)
    ll <- function(q) { # log likelihood function
        -sum(p * vismeteor::dvmgeom(m, lm, exp(q), log = TRUE))
    }
    est <- optim(1, ll, method='Brent', lower=0.01, upper=4, hessian=TRUE)
    expect_equal(round(est$par, 6), round(base::log(r), 6))

    # test q value of meteor magnitudes
    lm <- 5.8
    m <- as.integer(seq(8, -100, -1))
    p <- vismeteor::dvmgeom(m, lm, r)
    q.fun <- function(m) {
        q <- rep(0.0, length(m))
        idx <- m > -0.5
        q[idx] <- vismeteor::vmperception(m[idx], deriv = TRUE)/vismeteor::vmperception(m[idx])
        q
    }
    expect_true(abs(r - exp(sum(p * q.fun(lm - m)))) < 0.01)

    # density of meteor magnitudes equals geometric distribution
    perception.const <- function(m, log = FALSE) {
        rep(ifelse(log, 0.0, 1.0), length(m))
    }

    m <- as.integer(seq(6, -30, -1))
    p <- vismeteor::dvmgeom(m, lm, r, perception.fun = perception.const)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_equal(stats::dgeom(6 - m, 1 - 1/r), p)
})
