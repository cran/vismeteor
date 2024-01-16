test_that("qmideal", {
    psi <- 4.0

    p <- vismeteor::pmideal(6, psi, lower.tail = FALSE)
    m <- suppressWarnings(
        vismeteor::qmideal(c(-0.1, 0.0, p, 1.0, 1.2), psi, lower.tail = FALSE)
    )
    expect_length(m, 5)
    expect_equal(m, c(NA, Inf, 6, -Inf, NA))

    p <- vismeteor::pmideal(6, psi, lower.tail = TRUE)
    m <- suppressWarnings(
        vismeteor::qmideal(c(-0.1, 0.0, p, 1.0, 1.2), psi, lower.tail = TRUE)
    )
    expect_length(m, 5)
    expect_equal(m, c(NA, -Inf, 6, Inf, NA))

    # quantile of meteor magnitudes equals exponential distribution

    # from documentation
    dp.fun <- function(m, psi) {
        r <- 10^0.4
        1.5 * log(r) * sqrt(r^(3 * psi + 2 * m)/((r^psi + r^m)^5))
    }

    m_expected <- seq(-20, -10, 1)
    p <- vismeteor::pmideal(m_expected, psi, lower.tail = TRUE)
    m <- vismeteor::qmideal(p, psi, lower.tail = TRUE)
    expect_equal(m, m_expected)
    p.max <- stats::integrate(function(m) dp.fun(m, psi), -Inf, -9)$value
    m <- -9 - stats::qexp(p/p.max, 0.4 * log(10), lower.tail = FALSE)
    expect_equal(round(m, 3), m_expected)

    m_expected <- seq(15, 20, 1)
    p <- vismeteor::pmideal(m_expected, psi, lower.tail = FALSE)
    m <- vismeteor::qmideal(p, psi, lower.tail = FALSE)
    expect_equal(m, m_expected)
    p.max <- stats::integrate(function(m) dp.fun(m, psi), 14, Inf)$value
    m <- 14 + stats::qexp(p/p.max, 1.5 * 0.4 * log(10), lower.tail = FALSE)
    expect_equal(round(m, 3), m_expected)
})
