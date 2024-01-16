test_that("freq.quantile", {

    # simple example I
    freq <- c(1,2,3,4,5,6,7,8,9)
    f <- vismeteor::freq.quantile(freq, 10)
    expect_type(f, 'integer')
    expect_true(is.factor(f))
    expect_equal(as.integer(f), c(1, 1, 1, 1, 2, 2, 3, 3, 3))
    a <- sapply(split(freq, f), sum)
    expect_true(all(a >= 10))

    # simple example II
    freq <- c(9,8,7,6,5,4,3,2,1)
    f <- vismeteor::freq.quantile(freq, 10)
    expect_type(f, 'integer')
    expect_true(is.factor(f))
    expect_equal(as.integer(f), c(1, 1, 2, 2, 3, 3, 3, 3, 3))
    a <- sapply(split(freq, f), sum)
    expect_true(all(a >= 10))

    # simple example III
    freq <- rep(10, 6)
    f <- vismeteor::freq.quantile(freq, 10)
    expect_type(f, 'integer')
    expect_true(is.factor(f))
    expect_equal(as.integer(f), c(1, 2, 3, 4, 5, 6))
    a <- sapply(split(freq, f), sum)
    expect_true(all(a >= 10))

    # simple example IV
    freq <- rep(5, 12)
    f <- vismeteor::freq.quantile(freq, 10)
    expect_type(f, 'integer')
    expect_true(is.factor(f))
    expect_equal(as.integer(f), rep(c(1, 2, 3, 4, 5, 6), each=2))
    a <- sapply(split(freq, f), sum)
    expect_true(all(a >= 10))

    # test overflow
    freq <- c(1, 10)
    f <- vismeteor::freq.quantile(freq, 10)
    expect_type(f, 'integer')
    expect_true(is.factor(f))
    expect_equal(as.integer(f), c(1, 1))
    a <- sapply(split(freq, f), sum)
    expect_true(all(a >= 10))

    # test single large frequency
    freq <- 12
    f <- vismeteor::freq.quantile(freq, 10)
    expect_type(f, 'integer')
    expect_true(is.factor(f))
    expect_equal(as.integer(f), 1)

    # test single small frequency
    freq <- 2
    f <- vismeteor::freq.quantile(freq, 10)
    expect_type(f, 'integer')
    expect_true(is.factor(f))
    expect_equal(as.integer(f), 1)
})
