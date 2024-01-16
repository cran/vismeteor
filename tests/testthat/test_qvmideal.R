test_that("qvmideal", {
    lm <- 6.0
    psi <- 4.0

    for (current.psi in c(-10.0, psi)) {
        if (psi == current.psi) {
            m <- seq(-10, 6)
        } else {
            m <- current.psi + seq(-10, 10)
        }

        # lower tail
        p <- vismeteor::pvmideal(m, lm, current.psi, lower.tail = FALSE)
        p <- c(-1.0, 0.0, p - 1e-07, p + 1e-07, 1.0)
        q <- suppressWarnings(
            vismeteor::qvmideal(p, lm, current.psi, lower.tail = FALSE)
        )
        expect_equal(q, c(NA, 6, m, m - 1, -Inf), label = paste('psi =', current.psi))

        # upper tail
        p <- vismeteor::pvmideal(m, lm, current.psi, lower.tail = TRUE)
        q <- suppressWarnings(
            vismeteor::qvmideal(c(-1.0, 1.0, p + 1e-07, p - 1e-07, 0.0), lm, current.psi, lower.tail = TRUE)
        )
        expect_equal(q, c(NA, 6, m, m - 1, -Inf), label = paste('psi =', current.psi))
    }

    lms <- seq(5.6, 6.5, 0.1)
    for (current.lm in lms) {
        q <- suppressWarnings(
            vismeteor::qvmideal(c(-1.0, 0.0, 1.0), current.lm, psi, lower.tail = TRUE)
        )
        expect_equal(q, c(NA, -Inf, 6), label = paste('lm =', current.lm))

        q <- suppressWarnings(
            vismeteor::qvmideal(c(-1.0, 0.0, 1.0), current.lm, psi, lower.tail = FALSE)
        )
        expect_equal(q, c(NA, 6, -Inf), label = paste('lm =', current.lm))
    }

    for (current.lm in lms) {
        p <- vismeteor::pvmideal(6, current.lm, psi, lower.tail = TRUE)
        q <- vismeteor::qvmideal(c(p + 1e-07, p -1e-07), current.lm, psi, lower.tail = TRUE)
        expect_equal(q, c(6, 5), label = paste('lm =', current.lm))

        p <- vismeteor::pvmideal(6, current.lm, psi, lower.tail = FALSE)
        q <- vismeteor::qvmideal(c(p - 1e-07, p + 1e-07), current.lm, psi, lower.tail = FALSE)
        expect_equal(q, c(6, 5), label = paste('lm =', current.lm))
    }

    for (current.lm in lms) {
        p <- vismeteor::pvmideal(5, current.lm, psi, lower.tail = TRUE)
        q <- vismeteor::qvmideal(c(p + 1e-07, p - 1e-07), current.lm, psi, lower.tail = TRUE)
        expect_equal(q, c(5, 4), label = paste('lm =', current.lm))

        p <- vismeteor::pvmideal(5, current.lm, psi, lower.tail = FALSE)
        q <- vismeteor::qvmideal(c(p - 1e-07, p + 1e-07), current.lm, psi, lower.tail = FALSE)
        expect_equal(q, c(5, 4), label = paste('lm =', current.lm))
    }

    # quantile of meteor magnitudes equals geometric distribution
    psi <- 30
    m <- as.integer(seq(-10, 6, 1))
    p <- round(vismeteor::pvmideal(m, lm, psi, lower.tail = TRUE), 6)
    m <- vismeteor::qvmideal(p, lm, psi, lower.tail = TRUE)
    expect_type(m, 'double')
    expect_length(m, length(m))
    expect_equal(qvmgeom(p, lm, 10^0.4, lower.tail = TRUE), m)

    m <- as.integer(seq(-10, 6, 1))
    p <- round(vismeteor::pvmideal(m, lm, psi, lower.tail = FALSE), 6)
    m <- vismeteor::qvmideal(p, lm, psi, lower.tail = FALSE)
    expect_type(m, 'double')
    expect_length(m, length(m))
    expect_equal(vismeteor::qvmgeom(p, lm, 10^0.4, lower.tail = FALSE), m)
})
