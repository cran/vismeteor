test_that("cvmideal", {
    psi <- 4.0
    lm <- 6.5

    # from documentation
    dp.fun <- function(m, psi) {
        r <- 10^0.4
        1.5 * log(r) * sqrt(r^(3 * psi + 2 * m)/((r^psi + r^m)^5))
    }

    # test convolution
    m <- seq(6, -15, -1)
    expected_p <- sum(sapply(m, function(m) {
        stats::integrate(function(m) dp.fun(m, psi), m - 0.5, m + 0.5)$value *
            vismeteor::vmperception(lm - m)
    }))

    p <- vismeteor::cvmideal(lm, psi)
    expect_type(p, 'double')
    expect_length(p, 1)
    expect_equal(round(p, 5), round(expected_p, 5))
    expect_equal(round(log(p), 3), round(log(expected_p), 3))

    p <- vismeteor::cvmideal(lm, psi, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, 1)
    expect_equal(round(p, 3), round(log(expected_p), 3))

    p <- vismeteor::cvmideal(c(-15, 20), psi)
    expect_type(p, 'double')
    expect_length(p, 2)
    expect_equal(round(p, 5), c(0.0, 1.0))

    p <- vismeteor::cvmideal(lm, c(-15, 20))
    expect_type(p, 'double')
    expect_length(p, 2)
    expect_equal(round(p, 5), c(1.0, 0.0))

    p <- suppressWarnings(
        vismeteor::cvmideal(c(-Inf, -Inf, Inf, Inf), c(-Inf, Inf, -Inf, Inf))
    )
    expect_type(p, 'double')
    expect_length(p, 4)
    expect_true(is.na(p[1]))
    expect_equal(p[2], 0.0)
    expect_equal(p[3], 1.0)
    expect_true(is.na(p[4]))

    p <- suppressWarnings(
        vismeteor::cvmideal(c(-Inf, -Inf, Inf, Inf), c(-Inf, Inf, -Inf, Inf), log = TRUE)
    )
    expect_type(p, 'double')
    expect_length(p, 4)
    expect_true(is.na(p[1]))
    expect_equal(p[2], -Inf)
    expect_equal(p[3], 0.0)
    expect_true(is.na(p[4]))
})
