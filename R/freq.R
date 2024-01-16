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
#' freq <- c(1,2,3,4,5,6,7,8,9)
#' cumsum(freq)
#' (f <- freq.quantile(freq, 10))
#' tapply(freq, f, sum)
#' @export
freq.quantile <- function(freq, min) {
    if (0 >= min) {
        stop(paste0('min must be greater than 0 instead of "', min, '"!'))
    }

    n.sum <- 0
    id.last <- 1L
    id <- integer(length(freq))
    for (i in seq_along(freq)) {
        n.i <- freq[i]
        if ((n.i + n.sum) < min) {
            id[i] <- id.last
            n.sum <- n.i + n.sum
        } else {
            id[i] <- id.last
            n.sum <- 0
            id.last <- id.last + 1L
        }
    }

    if (0 == n.sum) {
        id.last <- id.last - 1L
    }

    id.max <- id[length(id)]
    if (id.max > 1 & n.sum > 0 & sum(freq[id == id.last]) < min) {
        id[id == id.max] <- id.max - 1L
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
#' The algorithm is asymptotic. This means that the more meteors the table contains,
#' the more unbiased is the result of the rounding.
#' @return
#' A rounded contingency table of meteor magnitudes is returned.
#' @examples
#' # For example, create a contingency table of meteor magnitudes
#' mt <- as.table(matrix(
#'     c(
#'         0.0, 0.0, 2.5, 0.5, 0.0, 1.0,
#'         0.0, 1.5, 2.0, 0.5, 0.0, 0.0,
#'         1.0, 0.0, 0.0, 3.0, 2.5, 0.5
#'     ), nrow = 3, ncol = 6, byrow = TRUE
#' ))
#' colnames(mt) <- seq(6)
#' rownames(mt) <- c('A', 'B', 'C')
#' mt
#' margin.table(mt, 1)
#' margin.table(mt, 2)
#'
#' # contingency table with integer values
#' (mt.int <- vmtable(mt))
#' margin.table(mt.int, 1)
#' margin.table(mt.int, 2)
#' @export
vmtable <- function(mt) {
    if (! methods::is(mt, 'table')) {
        stop(paste0('Magnitude table is not a table!'))
    }

    if (2L != length(dim(mt))) {
        stop(paste0('Magnitude table is not two-dimensional!'))
    }

    mt.m <- as.matrix(mt)
    ncol.mt <- ncol(mt.m)
    # "fair" rounding. Determinate the direction
    from.right <- 0L != sum(mt.m) %% 2L
    if (from.right) {
        mt.m <- t(apply(mt.m, 1L, rev))
    }
    mt2c <- round(2L * mt.m) # "half" meteors


    # phase 1: round column margin and add dummy row
    margin.v <- as.integer(as.vector(margin.table(mt2c, 2L)))
    dummy.row <- diff(c(0L, sapply(cumsum(margin.v), function(freq) {
        2L * (freq %/% 2L) # round down
    }))) - margin.v
    mt2c <- rbind(mt2c, dummy.row) # add dummy row
    nrow.mt2c <- nrow(mt2c)

    # phase 2: cumsum column-wise and round
    mt2c.v <- as.integer(as.vector(mt2c)) # column-wise
    mt2c.cs <- cumsum(mt2c.v) # due to sum preserving rounding
    mt2c.cs <- mapply(function(freq, row.id) { # add row index
        c(freq = freq, row.id = row.id)
    }, mt2c.cs, rep(seq_len(nrow.mt2c), times=ncol.mt), SIMPLIFY = FALSE)
    rows.reminder <- rep(0L, nrow.mt2c) # row reminder
    # apply each frequency column-wise, starting from upper left
    mt2c.f <- sapply(mt2c.cs, function(cs) {
        row.id <- cs['row.id']
        freq <- cs['freq']
        freq.quotient <- freq %/% 2L
        freq.reminder <- freq %% 2L

        freq <- 2L * freq.quotient
        if (0L == freq.reminder) {
            return(freq)
        }

        row.reminder.org <- rows.reminder[row.id]
        row.reminder.new <- 1L
        if (1L == row.reminder.org) {
            freq <- freq + 2L
            row.reminder.new <- 0L
        }

        if (row.reminder.org != row.reminder.new) {
            rows.reminder[row.id] <<- row.reminder.new
        }

        return(freq)
    })

    mt2c.v <- diff(append(mt2c.f, 0L, after = 0L)) # inverse of cumsum
    mt2c.m <- matrix(mt2c.v, ncol = ncol.mt)
    mt2c.m <- utils::head(mt2c.m, -1L) # remove dummy.row
    if (from.right) {
        mt2c.m <- t(apply(mt2c.m, 1L, rev))
    }

    result <- as.table(mt2c.m %/% 2L)
    dimnames(result) <- dimnames(mt) # restore dimnames

    return(result)
}
