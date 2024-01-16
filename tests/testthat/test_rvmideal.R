test_that("rvmideal", {
    psi <- 4.0

    with_seed <- function(seed, code) {
        code <- substitute(code)
        orig.seed <- .Random.seed
        on.exit(.Random.seed <<- orig.seed)
        set.seed(seed)
        eval.parent(code)
    }

    # meteor magnitudes <= 6.0
    m <- vismeteor::rvmideal(1000, 6.5, psi)
    expect_type(m, 'double')
    expect_length(m, 1000)
    expect_false(anyNA(m))
    expect_false(any(is.infinite(m)))
    expect_true(all(m <= 6.0))
    expect_equal(as.integer(m), m)

    # meteor magnitudes equals geometric distribution
    lm <- 6.3
    m <- with_seed(8, vismeteor::rvmideal(100000, lm, 30))
    expect_type(m, 'double')
    expect_length(m, 100000)
    expect_false(anyNA(m))
    expect_false(any(is.infinite(m)))
    expect_true(all(m <= 6))
    expect_equal(as.integer(m), m)
    m.geom <- with_seed(8, rvmgeom(100000, lm, 10^(1/2.5)))
    expect_equal(round(mean(m.geom), 2), round(mean(m), 2))
})