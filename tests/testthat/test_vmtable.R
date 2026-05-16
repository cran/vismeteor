test_that("vmtable", {
    test_mt <- function(mt) {
        mt_int <- vismeteor::vmtable(mt)
        expect_type(mt_int, "integer")
        expect_true(isa(mt_int, "table"))
        expect_equal(
            as.vector(margin.table(mt, 1L)),
            as.vector(margin.table(mt_int, 1L))
        )

        diff <- as.vector(margin.table(mt, 2L)) -
            as.vector(margin.table(mt_int, 2L))
        expect_true(any(abs(diff) <= 0.5))
        expect_true(0.0 == sum(mt - mt_int))
        expect_false(any(abs(mt - mt_int) > 1.0))

        mt_int
    }

    # null observation
    mt <- as.table(matrix(
        c(
            0.0, 0.0, 0.0, 0.0
        ),
        nrow = 1, ncol = 4, byrow = TRUE
    ))
    mt_int <- test_mt(mt)
    expected <- as.table(matrix(
        c(
            0L, 0L, 0L, 0L
        ),
        nrow = 1, ncol = 4, byrow = TRUE
    ))
    expect_equal(expected, mt_int)

    # single observations
    mt <- as.table(matrix(
        c(
            1.0, 0.5, 0.5, 0.0
        ),
        nrow = 1, ncol = 4, byrow = TRUE
    ))
    mt_int <- test_mt(mt)
    expected <- as.table(matrix(
        c(
            1L, 0L, 1L, 0L
        ),
        nrow = 1, ncol = 4, byrow = TRUE
    ))
    expect_equal(expected, mt_int)

    mt <- as.table(matrix(
        c(
            1.5, 0.5, 0.5, 0.5
        ),
        nrow = 1, ncol = 4, byrow = TRUE
    ))
    mt_int <- test_mt(mt)
    expected <- as.table(matrix(
        c(
            2L, 0L, 1L, 0L
        ),
        nrow = 1, ncol = 4, byrow = TRUE
    ))
    expect_equal(expected, mt_int)

    # multi observation I
    mt <- as.table(matrix(
        c(
            0.0, 0.0, 0.0, 1.0,
            0.0, 0.0, 0.5, 0.5,
            0.0, 0.5, 0.5, 0.0,
            0.0, 0.0, 0.0, 0.0,
            0.5, 0.5, 0.0, 0.0,
            1.0, 0.0, 0.0, 0.0
        ),
        nrow = 6, ncol = 4, byrow = TRUE
    ))
    mt_int <- test_mt(mt)
    expected <- as.table(matrix(
        c(
            0L, 0L, 0L, 1L,
            0L, 0L, 1L, 0L,
            0L, 1L, 0L, 0L,
            0L, 0L, 0L, 0L,
            1L, 0L, 0L, 0L,
            1L, 0L, 0L, 0L
        ),
        nrow = 6, ncol = 4, byrow = TRUE
    ))
    expect_equal(expected, mt_int)

    # multi observation II
    mt <- as.table(matrix(
        c(
            0.0, 0.0, 0.0, 1.0,
            0.5, 0.5, 0.5, 0.5,
            0.0, 0.5, 1.0, 0.5,
            0.0, 0.0, 0.0, 0.0,
            0.5, 1.0, 0.5, 0.0,
            1.0, 2.0, 3.0, 4.0
        ),
        nrow = 6, ncol = 4, byrow = TRUE
    ))
    mt_int <- test_mt(mt)
    expected <- as.table(matrix(
        c(
            0L, 0L, 0L, 1L,
            1L, 0L, 1L, 0L,
            0L, 1L, 0L, 1L,
            0L, 0L, 0L, 0L,
            0L, 1L, 1L, 0L,
            1L, 2L, 3L, 4L
        ),
        nrow = 6, ncol = 4, byrow = TRUE
    ))
    expect_equal(expected, mt_int)

    # Perseids margin check
    mt <- PER_2015_magn$magnitudes
    mt_margin <- as.data.frame(margin.table(mt, 1))
    colnames(mt_margin)[2] <- "freq0"
    mt_int <- test_mt(mt)

    parts <- 20L
    # split by row
    f <- floor(parts * (seq(1L, nrow(mt)) - 1L) / nrow(mt))
    for (i in (seq_len(parts) - 1L)) {
        mti <- mt[i == f, ]
        test_mt(mti)
        mti_int <- mt_int[i == f, ]
        diff <- as.vector(margin.table(mti, 2L)) - as.vector(margin.table(mti_int, 2L))
        expect_true(any(abs(diff) <= 0.5))
        diff <- sum(mti - mti_int)
        expect_true(0.0 == diff)
    }

    if (FALSE) {
        # create meteors from gamma distribution
        obs_n <- 200L
        m_n <- 1000L
        set.seed(13)
        m_rnd <- rgamma(m_n, 3.0)
        m_gamma_mean <- mean(m_rnd)
        m_gamma_var <- var(m_rnd)
        expect_true(abs(m_gamma_mean - 3.0) < 0.03)
        expect_true(abs(m_gamma_var - 3.0) < 0.1)
        m_rnd <- round(2.0 * m_rnd) / 2.0 # round 0.5
        m_obs <- seq(1, obs_n)
        names(m_rnd) <- rep(m_obs, each = m_n %/% obs_n)
        m_df <- do.call(
            rbind.data.frame,
            tapply(m_rnd, names(m_rnd), function(m) {
                obs <- as.integer(names(m)[1])
                m <- c(floor(m), floor(m + 0.5))
                t <- table(m) / 2.0
                df <- as.data.frame(t)
                df$m <- as.integer(levels(df$m))
                df$obs <- rep(obs, nrow(df))
                df
            })
        )
        m_df$m <- factor(as.character(m_df$m), levels = as.character(sort(unique(m_df$m))), ordered = TRUE)
        mt <- xtabs(Freq ~ obs + m, data = m_df)
        mt_r <- margin.table(mt, 1)
        mt_c <- margin.table(mt, 2)
        mt_mean <- sum(as.integer(names(mt_c)) * mt_c) / m_n
        mt_var <- sum((as.integer(names(mt_c)) - mt_mean)^2 * mt_c) / (m_n - 1)
        expect_true(abs(mt_mean - m_gamma_mean) < 0.01)
        expect_true(abs(mt_var - m_gamma_var) < 0.2)

        mt_int <- test_mt(mt)
        mt_int_c <- margin.table(mt_int, 2)
        mt_int_mean <- sum(as.integer(names(mt_int_c)) * mt_int_c) / m_n
        mt_int_var <- sum((as.integer(names(mt_int_c)) - mt_mean)^2 * mt_int_c) / (m_n - 1)
        expect_true(abs(mt_int_mean - m_gamma_mean) < 0.01)
        expect_true(abs(mt_int_var - m_gamma_var) < 0.2)
        expect_true(abs(mt_int_mean - mt_mean) < 0.01)
        expect_true(abs(mt_int_var - mt_var) < 0.03)
    }
})
