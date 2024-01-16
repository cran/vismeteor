test_that("vmtable", {

    test.mt <- function(mt) {
        mt.int <- vismeteor::vmtable(mt)
        expect_type(mt.int, 'integer')
        expect_true(isa(mt.int, 'table'))
        expect_equal(
            as.vector(margin.table(mt, 1L)),
            as.vector(margin.table(mt.int, 1L))
        )

        diff <- as.vector(margin.table(mt, 2L)) -
            as.vector(margin.table(mt.int, 2L))
        expect_true(any(abs(diff) <= 0.5))
        expect_true(0.0 == sum(mt - mt.int))
        expect_false(any(abs(mt - mt.int) > 1.0))

        return(mt.int)
    }

    # null observation
    mt <- as.table(matrix(
        c(
            0.0, 0.0, 0.0, 0.0
        ), nrow = 1, ncol = 4, byrow = TRUE
    ))
    mt.int <- test.mt(mt)
    expected <- as.table(matrix(
        c(
            0L, 0L, 0L, 0L
        ), nrow = 1, ncol = 4, byrow = TRUE
    ))
    expect_equal(expected, mt.int)

    # single observations
    mt <- as.table(matrix(
        c(
            1.0, 0.5, 0.5, 0.0
        ), nrow = 1, ncol = 4, byrow = TRUE
    ))
    mt.int <- test.mt(mt)
    expected <- as.table(matrix(
        c(
            1L, 0L, 1L, 0L
        ), nrow = 1, ncol = 4, byrow = TRUE
    ))
    expect_equal(expected, mt.int)

    mt <- as.table(matrix(
        c(
            1.5, 0.5, 0.5, 0.5
        ), nrow = 1, ncol = 4, byrow = TRUE
    ))
    mt.int <- test.mt(mt)
    expected <- as.table(matrix(
        c(
            2L, 0L, 1L, 0L
        ), nrow = 1, ncol = 4, byrow = TRUE
    ))
    expect_equal(expected, mt.int)

    # multi observation I
    mt <- as.table(matrix(
        c(
            0.0, 0.0, 0.0, 1.0,
            0.0, 0.0, 0.5, 0.5,
            0.0, 0.5, 0.5, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.5, 0.5, 0.0, 0.0,
            1.0, 0.0, 0.0, 0.0
        ), nrow = 6, ncol = 4, byrow = TRUE
    ))
    mt.int <- test.mt(mt)
    expected <- as.table(matrix(
        c(
            0L, 0L, 0L, 1L,
            0L, 0L, 1L, 0L,
            0L, 1L, 0L, 0L,
            0L, 0L, 0L, 0L,
            1L, 0L, 0L, 0L,
            1L, 0L, 0L, 0L
        ), nrow = 6, ncol = 4, byrow = TRUE
    ))
    expect_equal(expected, mt.int)

    # multi observation II
    mt <- as.table(matrix(
        c(
            0.0, 0.0, 0.0, 1.0,
            0.5, 0.5, 0.5, 0.5,
            0.0, 0.5, 1.0, 0.5,
            0.0, 0.0, 0.0, 0.0,
            0.5, 1.0, 0.5, 0.0,
            1.0, 2.0, 3.0, 4.0
        ), nrow = 6, ncol = 4, byrow = TRUE
    ))
    mt.int <- test.mt(mt)
    expected <- as.table(matrix(
        c(
            0L, 0L, 0L, 1L,
            1L, 0L, 1L, 0L,
            0L, 1L, 0L, 1L,
            0L, 0L, 0L, 0L,
            0L, 1L, 1L, 0L,
            1L, 2L, 3L, 4L
        ), nrow = 6, ncol = 4, byrow = TRUE
    ))
    expect_equal(expected, mt.int)

    # Perseids margin check
    mt <- PER_2015_magn$magnitudes
    mt.margin <- as.data.frame(margin.table(mt, 1))
    colnames(mt.margin)[2] <- 'freq0'
    mt.int <- test.mt(mt)

    parts <- 20L
    # split by row
    f <- floor(parts * (seq(1L, nrow(mt)) - 1L)/nrow(mt))
    for (i in (seq_len(parts) - 1L)) {
        mti <- mt[i == f,]
        test.mt(mti)
        mti.int <- mt.int[i == f,]
        diff <- as.vector(margin.table(mti, 2L)) - as.vector(margin.table(mti.int, 2L))
        expect_true(any(abs(diff) <= 0.5))
        diff <- sum(mti - mti.int)
        expect_true(0.0 == diff)
    }

    if (FALSE) {
        # create meteors from gamma distribution
        obs.n <- 200L
        m.n <- 1000L
        set.seed(13)
        m.rnd <- rgamma(m.n, 3.0)
        m.gamma.mean <- mean(m.rnd)
        m.gamma.var <- var(m.rnd)
        expect_true(abs(m.gamma.mean - 3.0) < 0.03)
        expect_true(abs(m.gamma.var - 3.0) < 0.1)
        m.rnd <- round(2.0 * m.rnd) / 2.0 # round 0.5
        m.obs <- seq(1, obs.n)
        names(m.rnd) <- rep(m.obs, each=m.n %/% obs.n)
        m.df <- do.call(
            rbind.data.frame,
            tapply(m.rnd, names(m.rnd), function(m){
                obs <- as.integer(names(m)[1])
                m <- c(floor(m), floor(m + 0.5))
                t <- table(m) / 2.0
                df <- as.data.frame(t)
                df$m <- as.integer(levels(df$m))
                df$obs <- rep(obs, nrow(df))
                df
            })
        )
        m.df$m <- factor(as.character(m.df$m), levels = as.character(sort(unique(m.df$m))), ordered = TRUE)
        mt <- xtabs(Freq ~ obs + m, data=m.df)
        mt.r <- margin.table(mt, 1)
        mt.c <- margin.table(mt, 2)
        mt.mean <- sum(as.integer(names(mt.c)) * mt.c)/m.n
        mt.var <- sum((as.integer(names(mt.c)) - mt.mean)^2 * mt.c)/(m.n - 1)
        expect_true(abs(mt.mean - m.gamma.mean) < 0.01)
        expect_true(abs(mt.var - m.gamma.var) < 0.2)

        mt.int <- test.mt(mt)
        mt.int.c <- margin.table(mt.int, 2)
        mt.int.mean <- sum(as.integer(names(mt.int.c)) * mt.int.c)/m.n
        mt.int.var <- sum((as.integer(names(mt.int.c)) - mt.mean)^2 * mt.int.c)/(m.n - 1)
        expect_true(abs(mt.int.mean - m.gamma.mean) < 0.01)
        expect_true(abs(mt.int.var - m.gamma.var) < 0.2)
        expect_true(abs(mt.int.mean - mt.mean) < 0.01)
        expect_true(abs(mt.int.var - mt.var) < 0.03)
    }
})