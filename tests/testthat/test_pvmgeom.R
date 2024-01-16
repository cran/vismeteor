test_that("pvmgeom", {
    r <- 1.8
    lm <- 6.5

    # from documentation
    dp.fun <- function(m) {
        stats::dgeom(m, 1 - 1/r) * vismeteor::vmperception(m + 0.5)
    }

    # probability lower tail

    m <- seq(40, 0)
    expected_p <- dp.fun(m)
    expected_p[1] <- expected_p[1] + stats::pgeom(m[1], 1 - 1/r, lower.tail = FALSE)
    expected_p <- expected_p/sum(expected_p)
    expected_p <- cumsum(expected_p) # lower tail
    m <- 6.0 - m

    p <- vismeteor::pvmgeom(m + 1, lm, r, lower.tail = TRUE)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_equal(p, expected_p)
    expect_equal(log(p), log(expected_p))

    p <- vismeteor::pvmgeom(m + 1, lm, r, lower.tail = TRUE, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_equal(p, log(expected_p))

    # probability upper tail

    m <- seq(0, 40)
    expected_p <- dp.fun(m)
    expected_p <- expected_p/sum(expected_p)
    expected_p <- cumsum(expected_p) # upper tail
    m <- 6.0 - m

    # probability
    p <- vismeteor::pvmgeom(m, lm, r, lower.tail = FALSE)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_equal(p, expected_p)
    expect_equal(log(p), log(expected_p))

    p <- vismeteor::pvmgeom(m, lm, r, lower.tail = FALSE, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_equal(p, log(expected_p))

    # sum lower.tail and upper.tail
    p <- vismeteor::pvmgeom(m, lm, r, lower.tail = FALSE) + vismeteor::pvmgeom(m, lm, r, lower.tail = TRUE)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_equal(p, rep(1.0, length(m)))
    expect_equal(log(p), rep(0.0, length(m)))

    p <- vismeteor::pvmgeom(6, c(5.4, 5.5, 5.6, 100), r, lower.tail = FALSE)
    expect_type(p, 'double')
    expect_length(p, 4)
    expect_equal(p[1], 0.0)
    expect_equal(p[2], 0.0)
    expect_equal(p[3], vismeteor::dvmgeom(6, 5.6, r))
    expect_equal(p[4], 1.0)

    p <- vismeteor::pvmgeom(6, c(5.4, 5.5, 5.6, 100), r, lower.tail = FALSE, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, 4)
    expect_equal(p[1], -Inf)
    expect_equal(p[2], -Inf)
    expect_equal(p[3], vismeteor::dvmgeom(6, 5.6, r, log = TRUE))
    expect_equal(p[4], 0.0)

    p <- vismeteor::pvmgeom(6, c(5.4, 5.5, 5.6, 100), r, lower.tail = TRUE)
    expect_type(p, 'double')
    expect_length(p, 4)
    expect_equal(p[1], 1.0)
    expect_equal(p[2], 1.0)
    expect_equal(p[3], 1 - vismeteor::dvmgeom(6, 5.6, r))
    expect_equal(p[4], 0.0)

    p <- vismeteor::pvmgeom(6, c(5.4, 5.5, 5.6, 100), r, lower.tail = TRUE, log = TRUE)
    expect_type(p, 'double')
    expect_length(p, 4)
    expect_equal(p[1], 0.0)
    expect_equal(p[2], 0.0)
    expect_equal(p[3], log(1 - vismeteor::dvmgeom(6, 5.6, r)))
    expect_lt(p[4], -50)

    # test probabilities of meteor magnitudes with different parameters
    p1 <- vismeteor::pvmgeom(3, 6.5, 1.8)
    p2 <- vismeteor::pvmgeom(4, 6.3, 2.0)
    p <- vismeteor::pvmgeom(c(3, 4), c(6.5, 6.3), c(1.8, 2.0))
    expect_equal(p, c(p1, p2))

    # test order of probabilities
    lm <- seq(2.5, 6.5, 0.1)
    p <- vismeteor::pvmgeom(rep(3L, length(lm)), lm, r, lower.tail = FALSE)
    expect_equal(p, p[order(p, decreasing = FALSE)])

    lm <- seq(2.5, 6.5, 0.1)
    p <- vismeteor::pvmgeom(rep(3L, length(lm)), lm, r, lower.tail = TRUE)
    expect_equal(p, p[order(p, decreasing = TRUE)])

    # test lower probabilities
    lm <- seq(5.4, 6.4, 0.1)
    p1 <- vismeteor::dvmgeom(rep(6L, length(lm)), lm, r)
    p2 <- vismeteor::pvmgeom(rep(6L, length(lm)), lm, r, lower.tail = FALSE)
    expect_equal(p1, p2)

    # test sum lower probabilities
    lm <- seq(5.4, 6.4, 0.1)
    p1 <- vismeteor::dvmgeom(rep(6L, length(lm)), lm, r) +
        vismeteor::dvmgeom(rep(5L, length(lm)), lm, r)
    p2 <- vismeteor::pvmgeom(rep(5L, length(lm)), lm, r, lower.tail = FALSE)
    expect_equal(p1, p2)

    # probability of meteor magnitudes equals geometric distribution

    perception.const <- function(m, log = FALSE) {
        rep(ifelse(log, 0.0, 1.0), length(m))
    }

    m <- seq(0, 30, 1)
    p <- vismeteor::pvmgeom(6 - m, 6.5, r, lower.tail = FALSE, perception.fun = perception.const)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_equal(p, stats::pgeom(m, 1 - 1/r, lower.tail = TRUE))

    m <- seq(0, 30, 1)
    p <- vismeteor::pvmgeom(6 - m, 6.5, r, lower.tail = TRUE, perception.fun = perception.const)
    expect_type(p, 'double')
    expect_length(p, length(m))
    expect_equal(p, stats::pgeom(m, 1 - 1/r, lower.tail = FALSE))
})
