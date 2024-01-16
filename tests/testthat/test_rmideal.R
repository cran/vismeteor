test_that("rmideal", {
    psi <- 4.0

    with_seed <- function(seed, code) {
        code <- substitute(code)
        orig.seed <- .Random.seed
        on.exit(.Random.seed <<- orig.seed)
        set.seed(seed)
        eval.parent(code)
    }

    m <- vismeteor::rmideal(1000, psi)
    expect_type(m, 'double')
    expect_length(m, 1000)
    expect_false(anyNA(m))
    expect_false(any(is.infinite(m)))

    # MLE of psi
    m <- with_seed(8, vismeteor::rmideal(100000, psi))
    expect_type(m, 'double')
    expect_length(m, 100000)
    expect_false(anyNA(m))
    expect_false(any(is.infinite(m)))
    llr <- function(psi) {
        -sum(vismeteor::dmideal(m, psi, log=TRUE))
    }
    est <- optim(2, llr, method='Brent', lower=0, upper=8, hessian=TRUE)
    expect_equal(round(est$par, 2), psi)
})