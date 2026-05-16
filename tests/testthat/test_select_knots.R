test_that("select_knots", {
    # synthetic, deterministic signal
    set.seed(1)
    n <- 200
    x <- seq(0, 10, length.out = n)
    y <- sin(x) + stats::rnorm(n, sd = 0.2)
    dat <- data.frame(x = x, y = y)

    fit <- \(d, knots) {
        f <- if (length(knots) == 0L) {
            y ~ splines::ns(x)
        } else {
            y ~ splines::ns(x, knots = knots)
        }
        stats::lm(f, data = d)
    }
    score_bic <- \(d, knots) stats::BIC(fit(d, knots))
    cand <- seq(1, 9, by = 0.5)

    # forward selection: result has the expected shape
    res <- select_knots(dat, cand, score_bic)
    expect_named(res, c(
        "backward", "knots", "fixed_knots",
        "score", "n_steps", "history"
    ))
    expect_false(res$backward)
    expect_true(is.numeric(res$knots))
    expect_true(is.numeric(res$fixed_knots))
    expect_equal(res$knots, sort(res$knots))
    expect_equal(res$fixed_knots, sort(res$fixed_knots))
    expect_true(all(res$knots %in% cand))
    expect_length(intersect(res$knots, res$fixed_knots), 0L)
    expect_null(res$n_steps)
    expect_equal(res$score, tail(res$history$score, 1L))
    expect_true(is.data.frame(res$history))
    expect_named(
        res$history,
        c("step", "n_knots", "changed_knot", "score", "extra")
    )
    # NULL-mode never records a worsening step
    expect_false(any(res$history$extra))

    # n_steps = 2L: exactly 2 total iterations (NOT 2 past optimum). In
    # forward mode without bulk that is 2 added knots and a 3-row history.
    res2 <- select_knots(dat, cand, score_bic, n_steps = 2L)
    expect_identical(res2$n_steps, 2L)
    expect_length(res2$knots, 2L)
    expect_equal(nrow(res2$history), 3L)

    # A large n_steps lets the search run through the optimum and onward;
    # the best point along that trajectory is at least as good as the
    # NULL-mode optimum (both follow the same path up to it).
    res_long <- select_knots(dat, cand, score_bic,
        n_steps = length(cand) + 5L
    )
    expect_lte(min(res_long$history$score), res$score)

    # Single-step exploration past the optimum: take one controlled step
    # from the previous optimum and check that exactly one new knot appears.
    anchor <- c(res$knots, res$fixed_knots)
    res_step1 <- select_knots(dat, cand, score_bic,
        fixed_knots = anchor, n_steps = 1L
    )
    expect_length(res_step1$knots, 1L)
    expect_false(res_step1$knots %in% anchor)
    expect_identical(sort(res_step1$fixed_knots), sort(anchor))

    # backward elimination produces a valid run. Greedy stepwise selection is
    # path-dependent (see the "When to use this" section of ?select_knots),
    # so backward and forward need not reach the same knot set; only the
    # structural invariants must hold.
    res_b <- select_knots(dat, cand, score_bic, backward = TRUE)
    expect_true(res_b$backward)
    expect_true(all(res_b$knots %in% cand))
    expect_equal(res_b$knots, sort(res_b$knots))
    expect_lt(res_b$score, score_bic(dat, numeric(0)))

    # bulk_gap = 0L disables bulk; it is a different trajectory than the
    # default bulk run, so optima can differ, but both must be valid.
    res_b0 <- select_knots(dat, cand, score_bic,
        backward = TRUE, bulk_gap = 0L
    )
    expect_true(all(res_b0$knots %in% cand))
    expect_lt(res_b0$score, score_bic(dat, numeric(0)))

    # score_fun that errors out is treated as Inf and the search still runs
    score_err <- \(d, knots) {
        if (length(knots) > 2L) stop("boom")
        score_bic(d, knots)
    }
    res_err <- select_knots(dat, cand, score_err)
    expect_true(is.finite(res_err$score))
    expect_lte(length(res_err$knots), 2L)

    # fixed_knots, forward: echoed back in $fixed_knots, never in $knots.
    res_ff <- select_knots(dat, cand, score_bic, fixed_knots = 5.25)
    expect_true(5.25 %in% res_ff$fixed_knots)
    expect_false(5.25 %in% res_ff$knots)
    expect_length(intersect(res_ff$knots, res_ff$fixed_knots), 0L)
    expect_equal(res_ff$knots, sort(res_ff$knots))

    # fixed_knots, backward: never removed -> still in $fixed_knots, not in $knots
    res_fb <- select_knots(dat, cand, score_bic,
        backward = TRUE, fixed_knots = 5.25
    )
    expect_true(5.25 %in% res_fb$fixed_knots)
    expect_false(5.25 %in% res_fb$knots)
    expect_equal(res_fb$knots, sort(res_fb$knots))

    # fixed_knots = cand in backward mode: nothing is removable, the search
    # terminates immediately. $knots is empty; $fixed_knots holds the whole pool.
    res_fall <- select_knots(dat, cand, score_bic,
        backward = TRUE, fixed_knots = cand
    )
    expect_length(res_fall$knots, 0L)
    expect_equal(sort(c(res_fall$knots, res_fall$fixed_knots)), sort(cand))

    # Edge case: forward starting from a fully-saturated fixed_knots set
    # (overfit start, no improving move possible). NULL-mode must terminate
    # at the starting state without recording any worsening step.
    res_sat <- select_knots(dat, cand, score_bic, fixed_knots = cand)
    expect_length(res_sat$knots, 0L)
    expect_equal(nrow(res_sat$history), 1L)
    expect_false(any(res_sat$history$extra))

    # Same start, but explicit n_steps = 2L: two worsening steps are taken.
    # (No candidates remain to add in forward mode -- the loop exits via the
    # `length(remaining) == 0` guard at the top. Test that branch instead by
    # using a set that overfits but still has candidates to spare.)
    overfit_anchor <- cand[seq(1L, length(cand) - 2L)]
    res_over_null <- select_knots(dat, cand, score_bic,
        fixed_knots = overfit_anchor
    )
    res_over_2 <- select_knots(dat, cand, score_bic,
        fixed_knots = overfit_anchor, n_steps = 2L
    )
    # With n_steps = 2L the loop runs both iterations regardless of improvement.
    expect_equal(nrow(res_over_2$history), 3L)
    # If the NULL run stopped at the start (no improvement possible), then the
    # explicit n_steps run necessarily explores worsening territory.
    if (nrow(res_over_null$history) == 1L) {
        expect_true(any(res_over_2$history$extra))
    }

    # score_fun receives sorted knots (no defensive sort needed in the user fit)
    seen_unsorted <- FALSE
    score_probe <- \(d, knots) {
        if (length(knots) >= 2L && is.unsorted(knots)) seen_unsorted <<- TRUE
        score_bic(d, knots)
    }
    select_knots(dat, cand, score_probe)
    expect_false(seen_unsorted)
    select_knots(dat, cand, score_probe, backward = TRUE)
    expect_false(seen_unsorted)

    # Names on knot_candidates and fixed_knots are preserved through the
    # returned $knots and $fixed_knots.
    cand_named <- setNames(cand, paste0("k", seq_along(cand)))
    res_named <- select_knots(dat, cand_named, score_bic)
    expect_false(is.null(names(res_named$knots)))
    expect_true(all(names(res_named$knots) %in% names(cand_named)))

    fk_named <- c(anchor = 5.25)
    res_fnamed <- select_knots(dat, cand_named, score_bic,
        fixed_knots = fk_named
    )
    expect_true("anchor" %in% names(res_fnamed$fixed_knots))
    expect_equal(
        unname(res_fnamed$fixed_knots[
            names(res_fnamed$fixed_knots) == "anchor"
        ]),
        5.25
    )

    # input validation
    expect_error(select_knots(dat, cand, score_bic, backward = NA))
    expect_error(select_knots(dat, cand, score_bic, n_steps = -1L))
    expect_error(select_knots(dat, cand, score_bic, n_steps = 0L))
    expect_error(select_knots(dat, cand, score_bic, n_steps = 1.5))
    expect_error(select_knots(dat, cand, score_bic, n_steps = NA_integer_))
    expect_error(select_knots(dat, cand, score_bic, n_cores = 0L))
    expect_error(select_knots(dat, cand, score_fun = "not a function"))
    expect_error(select_knots(dat, cand, score_bic, fixed_knots = c(1, NA)))
})
