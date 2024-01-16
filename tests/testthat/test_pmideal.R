test_that("pmideal", {
    psi <- 4.0

    # from documentation
    dp.fun <- function(m, psi) {
        r <- 10^0.4
        1.5 * log(r) * sqrt(r^(3 * psi + 2 * m)/((r^psi + r^m)^5))
    }

    # probability of m >= 3
    expected_p <- stats::integrate(function(m) dp.fun(m, psi), -Inf, 3)$value
    expect_equal(vismeteor::pmideal(3, psi, lower.tail = TRUE), expected_p)

    # probability
    p <- vismeteor::pmideal(c(3, -Inf, -100, 100, Inf), psi, lower.tail = TRUE)
    expect_type(p, 'double')
    expect_length(p, 5)
    expect_equal(p[1], expected_p)
    expect_equal(p[2], 0.0)
    expect_equal(p[3], 0.0)
    expect_equal(p[4], 1.0)
    expect_equal(p[5], 1.0)

    p <- vismeteor::pmideal(c(3, -Inf, -100, 100, Inf), psi, lower.tail = FALSE)
    expect_type(p, 'double')
    expect_length(p, 5)
    expect_equal(round(p[1], 6), round(1.0 - expected_p, 6))
    expect_equal(p[2], 1.0)
    expect_equal(p[3], 1.0)
    expect_equal(p[4], 0.0)
    expect_equal(p[5], 0.0)

    # log probability
    p <- vismeteor::pmideal(c(3, -Inf, -100, 100, Inf), psi, lower.tail = TRUE, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, 5)
    expect_equal(round(p[1], 6), round(log(expected_p), 6))
    expect_equal(p[2], -Inf)
    expect_equal(round(p[3], 6), -95.382075)
    expect_equal(p[4], 0)
    expect_equal(p[5], 0)

    p <- vismeteor::pmideal(c(3, -Inf, -100, 100, Inf), psi, lower.tail = FALSE, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, 5)
    expect_equal(round(p[1], 6), round(log(1.0 - expected_p), 6))
    expect_equal(p[2], 0.0)
    expect_equal(p[3], 0.0)
    expect_equal(round(p[4], 6), -132.628899)
    expect_equal(p[5], -Inf)

    # density of meteor magnitudes equals exponential distribution
    m <- seq(-20, -10, 1)
    p <- vismeteor::pmideal(m, psi, lower.tail = TRUE)
    p.max <- stats::integrate(function(m) dp.fun(m, psi), -Inf, -10)$value
    expect_equal(
        round(log(p), 3),
        round(log(p.max * stats::pexp(-m - 10, 0.4 * log(10), lower.tail = FALSE)), 3)
    )

    p <- vismeteor::pmideal(m, psi, lower.tail = TRUE, log = TRUE)
    expect_equal(
        round(p, 3),
        round(log(p.max) + stats::pexp(-m - 10, 0.4 * log(10), lower.tail = FALSE, log = TRUE), 3)
    )

    m <- seq(15, 20, 1)
    p <- vismeteor::pmideal(m, psi, lower.tail = FALSE)
    p.max <- stats::integrate(function(m) dp.fun(m, psi), 15, Inf)$value
    expect_equal(
        round(log(p), 3),
        round(log(p.max * stats::pexp(m - 15, 1.5 * 0.4 * log(10), lower.tail = FALSE)), 3)
    )

    p <- vismeteor::pmideal(m, psi, lower.tail = FALSE, log = TRUE)
    expect_equal(
        round(p, 3),
        round(log(p.max) + stats::pexp(m - 15, 1.5 * 0.4 * log(10), lower.tail = FALSE, log = TRUE), 3)
    )
})
