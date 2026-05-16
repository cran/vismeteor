#' @title Quantiles with a minimum frequency
#' @description
#' This function generates quantiles with a minimum frequency.
#' These quantiles are formed from a vector `freq` of frequencies.
#' Each quantile then has the minimum total frequency `min`.
#' @param freq integer; A vector of frequencies.
#' @param min integer; Minimum total frequency per quantile.
#' @details
#' The frequencies `freq` are grouped in the order in which they
#' are passed as a vector.
#' The minimum `min` must be greater than `0`.
#' @return
#' A factor of indices is returned.
#' The index references the corresponding passed frequency `freq`.
#' @examples
#' freq <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
#' cumsum(freq)
#' (f <- freq_quantile(freq, 10))
#' tapply(freq, f, sum)
#' @export
freq_quantile <- function(freq, min) {
    if (0 >= min) {
        stop(paste0('min must be greater than 0 instead of "', min, '"!'))
    }

    n_sum <- 0
    id_last <- 1L
    id <- integer(length(freq))
    for (i in seq_along(freq)) {
        n_i <- freq[i]
        if ((n_i + n_sum) < min) {
            id[i] <- id_last
            n_sum <- n_i + n_sum
        } else {
            id[i] <- id_last
            n_sum <- 0
            id_last <- id_last + 1L
        }
    }

    if (0 == n_sum) {
        id_last <- id_last - 1L
    }

    id_max <- id[length(id)]
    if (id_max > 1 && n_sum > 0 && sum(freq[id == id_last]) < min) {
        id[id == id_max] <- id_max - 1L
    }

    factor(id, ordered = TRUE)
}

#' @title Rounds a contingency table of meteor magnitude frequencies
#' @description
#' The meteor magnitude contingency table of VMDB contains half meteor counts (e.g. `3.5`).
#' This function converts these frequencies to integer values.
#' @param mt table; A two-dimensional contingency table of meteor magnitude frequencies.
#' @details
#' The contingency table of meteor magnitudes `mt` must be two-dimensional.
#' The row names refer to the magnitude observations.
#' Column names must be integer meteor magnitude values.
#' Also, the columns must be sorted in ascending or descending order of meteor magnitude.
#'
#' A sum-preserving algorithm is used for rounding.
#' It ensures that the total frequency of meteors per observation is preserved.
#' The marginal frequencies of the magnitudes are also preserved with
#' the restriction that the deviation is at most \eqn{\pm 0.5}.
#' If the total sum of a meteor magnitude is integer,
#' then the deviation is \eqn{\pm 0}.
#'
#' The algorithm is unbiased: for a fixed observation order it preserves the original
#' totals without introducing systematic drift, even though each run follows the
#' deterministic sequence dictated by the observed counts and their ordering.
#'
#' @note Internally the counts are doubled to half-meteor units, leftover halves are
#' alternated between rows so column margins stay within \eqn{\pm 0.5}, and when the
#' grand total is odd the matrix is temporarily mirrored so the unavoidable surplus
#' meteor originates from the opposite end of the magnitude scale rather than always
#' favouring the faintest bin. The mirroring is only the initial condition; the loop
#' then processes the table cell by cell so the rounding direction alternates between
#' bright and faint magnitudes depending on the current row and column state.
#'
#' @return
#' A rounded contingency table of meteor magnitudes is returned.
#' @examples
#' # For example, create a contingency table of meteor magnitudes
#' mt <- as.table(matrix(
#'     c(
#'         0.0, 0.0, 2.5, 0.5, 0.0, 1.0,
#'         0.0, 1.5, 2.0, 0.5, 0.0, 0.0,
#'         1.0, 0.0, 0.0, 3.0, 2.5, 0.5
#'     ),
#'     nrow = 3, ncol = 6, byrow = TRUE
#' ))
#' colnames(mt) <- seq(6)
#' rownames(mt) <- c("A", "B", "C")
#' mt
#' margin.table(mt, 1)
#' margin.table(mt, 2)
#'
#' # contingency table with integer values
#' (mt_int <- vmtable(mt))
#' margin.table(mt_int, 1)
#' margin.table(mt_int, 2)
#' @export
vmtable <- function(mt) {
    if (!methods::is(mt, "table")) {
        stop(paste0("Magnitude table is not a table!"))
    }

    if (2L != length(dim(mt))) {
        stop(paste0("Magnitude table is not two-dimensional!"))
    }

    mt_m <- as.matrix(mt) # work with raw numeric matrix to control rounding precisely
    ncol_mt <- ncol(mt_m)
    # "Fair" rounding. When the total count is odd we mirror the columns so that the surplus
    # meteor is taken from the opposite end of the magnitude scale instead of consistently
    # favouring the faintest magnitudes during the rounding pass
    from_right <- 0L != sum(mt_m) %% 2L
    if (from_right) {
        mt_m <- t(apply(mt_m, 1L, rev))
    }
    mt2c <- round(2L * mt_m) # work in half-meteor units so integer rounding keeps the original totals

    # phase 1: round column margin and add dummy row
    margin_v <- as.integer(as.vector(margin.table(mt2c, 2L)))
    # dummy row absorbs the 0.5 remainder per column so every column margin stays intact after rounding
    dummy_row <- diff(c(0L, sapply(cumsum(margin_v), \(freq) {
        2L * (freq %/% 2L) # round down
    }))) - margin_v
    mt2c <- rbind(mt2c, dummy_row) # add dummy row
    nrow_mt2c <- nrow(mt2c)

    # phase 2: cumsum column-wise and round
    mt2c_v <- as.integer(as.vector(mt2c)) # traverse cells column-by-column to match VMDB layout
    mt2c_v_cs <- cumsum(mt2c_v) # cumulative sum acts as running total for the sum-preserving rounding logic
    # retain row membership so remainders alternate per observation
    mt2c_row_ids <- rep(seq_len(nrow_mt2c), times = ncol_mt)

    rows_reminder <- rep(0L, nrow_mt2c) # tracks whether a row already kept the previous 0.5 remainder
    mt2c_f <- integer(length(mt2c_v_cs))
    # walk through the cumulative frequencies and distribute every odd remainder to alternating rows
    for (i in seq_along(mt2c_v_cs)) {
        row_id <- mt2c_row_ids[i]
        freq <- mt2c_v_cs[i] # running total for this cell including all previous entries of the column
        freq_quotient <- freq %/% 2L # largest whole-meteor multiple we can commit without breaking the total
        freq_reminder <- freq %% 2L # 1 when the cumulative sum would otherwise leave a dangling half meteor

        freq <- 2L * freq_quotient
        if (0L != freq_reminder) {
            row_reminder_org <- rows_reminder[row_id]
            row_reminder_new <- 1L
            if (1L == row_reminder_org) {
                freq <- freq + 2L # give the extra meteor to the opposite row on alternating encounters
                row_reminder_new <- 0L
            }

            if (row_reminder_org != row_reminder_new) {
                rows_reminder[row_id] <- row_reminder_new # remember which rows have already claimed the half meteor
            }
        }
        mt2c_f[i] <- freq # store the adjusted cumulative total for later differencing
    }

    mt2c_v <- diff(append(mt2c_f, 0L, after = 0L)) # undo cumulative form to recover per-cell counts
    mt2c_m <- matrix(mt2c_v, ncol = ncol_mt) # restore table shape after column-wise walk
    mt2c_m <- utils::head(mt2c_m, -1L) # remove dummy_row
    if (from_right) {
        mt2c_m <- t(apply(mt2c_m, 1L, rev)) # reapply original magnitude ordering if we flipped earlier
    }

    result <- as.table(mt2c_m %/% 2L) # collapse half-meteor counts back to whole meteors
    dimnames(result) <- dimnames(mt) # restore dimnames

    result
}
