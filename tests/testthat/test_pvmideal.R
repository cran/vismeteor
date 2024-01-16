test_that("pvmideal", {
    lm <- 6.5
    psi <- 4.0

    # from documentation
    dp.fun <- function(m, psi) {
        r <- 10^0.4
        1.5 * log(r) * sqrt(r^(3 * psi + 2 * m)/((r^psi + r^m)^5))
    }

    for (current.psi in c(-10.0, psi, 5)) {
        # test probabilities lower tail

        m <- seq(current.psi - 25, 6)
        expected_p <- sapply(m, function(m) {
            stats::integrate(function(m) dp.fun(m, current.psi), m - 0.5, m + 0.5)$value *
                vismeteor::vmperception(lm - m)
        })
        expected_p[1] <- expected_p[1] + vismeteor::pmideal(m[1] - 0.5, current.psi, lower.tail = TRUE)
        expected_p <- cumsum(expected_p)/sum(expected_p) # lower tail

        p <- vismeteor::pvmideal(m + 1, lm, current.psi, lower.tail = TRUE)
        expect_type(p, 'double')
        expect_length(p, length(m))
        expect_false(any(abs(log(p) - log(expected_p)) > 0.02), label = paste('psi =', current.psi))

        p <- vismeteor::pvmideal(m + 1, lm, current.psi, lower.tail = TRUE, log = TRUE)
        expect_type(p, 'double')
        expect_length(p, length(m))
        expect_false(any(abs(p - log(expected_p)) > 0.02), label = paste('psi =', current.psi))

        # test probabilities upper tail

        m <- seq(6, current.psi - 25, -1)
        expected_p <- sapply(m, function(m) {
            stats::integrate(function(m) dp.fun(m, current.psi), m - 0.5, m + 0.5)$value *
                vismeteor::vmperception(lm - m)
        })
        expected_p <- cumsum(expected_p)/sum(expected_p) # upper tail

        p <- vismeteor::pvmideal(m, lm, current.psi, lower.tail = FALSE)
        expect_type(p, 'double')
        expect_length(p, length(m))
        expect_false(any(abs(log(p) - log(expected_p)) > 0.02), label = paste('psi =', current.psi))

        p <- vismeteor::pvmideal(m, lm, current.psi, lower.tail = FALSE, log = TRUE)
        expect_type(p, 'double')
        expect_length(p, length(m))
        expect_false(any(abs(p - log(expected_p)) > 0.02), label = paste('psi =', current.psi))

        # test sum is 1.0
        m <- seq(current.psi - 25, 6)
        p <- vismeteor::pvmideal(m, lm, current.psi, lower.tail = TRUE) +
            vismeteor::pvmideal(m, lm, current.psi, lower.tail = FALSE)
        expect_type(p, 'double')
        expect_length(p, length(m))
        expect_false(any(abs(1.0 - p) > 1e-06), label = paste('psi =', current.psi))
    }

    # probability of m >= 3
    expected_p <- round(sum(vismeteor::dvmideal(as.integer(seq(3, 6)), 6.3, psi)), 6)
    expect_equal(expected_p, 0.348833)

    # probability
    p <- vismeteor::pvmideal(c(3, 7, -100, -Inf), 6.3, psi, lower.tail = TRUE)
    expect_type(p, 'double')
    expect_length(p, 4)
    expect_equal(round(p[1], 6), 1 - expected_p)
    expect_equal(p[2], 1.0)
    expect_equal(p[3], 0.0)
    expect_equal(p[4], 0.0)

    p <- vismeteor::pvmideal(c(3, 7, -100, -Inf), 6.3, psi, lower.tail = FALSE)
    expect_type(p, 'double')
    expect_length(p, 4)
    expect_equal(round(p[1], 6), expected_p)
    expect_equal(p[2], 0.0)
    expect_equal(p[3], 1.0)
    expect_equal(p[4], 1.0)

    # log probability
    p <- vismeteor::pvmideal(c(3, 7, -Inf), 6.3, psi, lower.tail = TRUE, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, 3)
    expect_equal(round(p[1], 4), round(log(1 - expected_p), 4))
    expect_equal(p[2], 0.0)
    expect_equal(p[3], -Inf)

    p <- vismeteor::pvmideal(c(3, 7, -Inf), 6.3, psi, lower.tail = FALSE, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, 3)
    expect_equal(round(p[1], 4), round(log(expected_p), 4))
    expect_equal(p[2], -Inf)
    expect_equal(p[3], 0.0)

    # test order of probabilities
    lm <- seq(2.5, 6.5, 0.1)
    p <- vismeteor::pvmideal(rep(3L, length(lm)), lm, psi, lower.tail = TRUE)
    expect_equal(p, p[order(p, decreasing = TRUE)])

    lm <- seq(2.5, 6.5, 0.1)
    p <- vismeteor::pvmideal(rep(3L, length(lm)), lm, psi, lower.tail = FALSE)
    expect_equal(p, p[order(p, decreasing = FALSE)])

    # test lower probabilities
    lm <- seq(5.4, 6.4, 0.1)
    p1 <- vismeteor::dvmideal(rep(6L, length(lm)), lm, psi)
    p2 <- vismeteor::pvmideal(rep(6L, length(lm)), lm, psi, lower.tail = FALSE)
    expect_equal(p1, p2)

    # test sum lower probabilities
    lm <- seq(5.4, 6.4, 0.1)
    p1 <- vismeteor::dvmideal(rep(6L, length(lm)), lm, psi) +
        vismeteor::dvmideal(rep(5L, length(lm)), lm, psi)
    p2 <- vismeteor::pvmideal(rep(5L, length(lm)), lm, psi, lower.tail = FALSE)
    expect_equal(p1, p2)

    # density of meteor magnitudes equals geometric distribution
    lm <- 6.3
    m <- as.integer(seq(-20, 6))

    p <- vismeteor::pvmideal(m, lm, 16.25, lower.tail = TRUE, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_equal(p, pvmgeom(m, lm, 10^0.4, lower.tail = TRUE, log = TRUE))

    p <- vismeteor::pvmideal(m, lm, 16.25, lower.tail = FALSE, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_equal(p, pvmgeom(m, lm, 10^0.4, lower.tail = FALSE, log = TRUE))

    p <- vismeteor::pvmideal(m, lm, 16.35, lower.tail = TRUE, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_equal(p, pvmgeom(m, lm, 10^0.4, lower.tail = TRUE, log = TRUE))

    p <- vismeteor::pvmideal(m, lm, 16.35, lower.tail = FALSE, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_equal(p, pvmgeom(m, lm, 10^0.4, lower.tail = FALSE, log = TRUE))
})
