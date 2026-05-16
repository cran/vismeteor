#' @title Greedy stepwise knot selection for a regression spline
#' @description
#' Selects a parsimonious subset of \code{knot_candidates} as interior knots
#' for a regression spline by greedy forward addition or backward elimination,
#' scored by a user-supplied criterion (BIC, AIC, or any cross-validation
#' function). The function is model-agnostic: it only chooses which
#' knot positions to include and passes that vector to your \code{score_fun};
#' your \code{score_fun} decides which basis (\code{splines::ns},
#' \code{splines::bs} of any degree, …) and which model family to fit.
#' Typical use: fitting a smooth, parsimonious trend through noisy data — for
#' example an activity profile of a meteor shower across solar longitude —
#' where the number and position of knots should remain small and inspectable.
#'
#' @details
#' The typical use case --- meteor shower rates over the solar longitude or
#' time, for example --- is a smooth, slowly varying signal whose shape is
#' not known in advance and whose curvature changes locally: too rigid for
#' a low-order polynomial, too noisy for a histogram. A regression spline
#' (typically cubic, e.g. \code{splines::ns} or
#' \code{splines::bs(degree = 3)}) with a modest number of well-placed
#' interior knots gives a smooth, locally flexible fit with interpretable
#' degrees of freedom. Picking that number is a bias/variance trade-off:
#' too many knots overfit (ringing, unstable derivatives, slow fits), too
#' few introduce bias. \code{select_knots()} automates the trade-off by
#' greedy forward addition or backward elimination, scored by a
#' user-defined criterion; you control which knots are even allowed
#' (\code{knot_candidates}) and what counts as "better" (\code{score_fun}).
#'
#' Your \code{score_fun} must take \code{(data, knots)} and return a single
#' \code{numeric} value; \emph{lower is better}. Typically it fits a model
#' with these interior knots (\code{length(knots) == 0L} means "no interior
#' knots") and reports an information criterion (\code{stats::BIC},
#' \code{stats::AIC}) or a held-out score. \code{knots} arrives sorted, so
#' the fit code does not need to sort it again. A divergent or failed fit
#' should ideally return \code{Inf}; the function additionally wraps each
#' call in \code{tryCatch()} and treats errors as \code{Inf} so the search
#' continues robustly. See \code{Examples} for a runnable template.
#'
#' Backward elimination can drop multiple knots per round in "bulk" mode.
#' When \code{backward = TRUE} and \code{bulk_gap >= 1L}, each round removes
#' several well-separated knots at once (minimum index gap \code{bulk_gap}
#' in the current sorted knot list). For a B-spline basis of degree
#' \code{d}, each basis function has support over \code{d + 1} consecutive
#' knot intervals, so removing knots that are at least \code{d + 1}
#' positions apart has nearly additive effect on the score --- giving
#' roughly a \code{bulk_gap}-fold speed-up at the cost of a small
#' approximation (mitigated by a verify-fit after each bulk round). The
#' default \code{4L} matches \code{d + 1} for cubic bases such as
#' \code{splines::ns} or \code{splines::bs(degree = 3)}; for other bases
#' pick \code{bulk_gap = degree + 1}, or set \code{bulk_gap = 0L} for
#' strict one-knot-per-round behaviour.
#'
#' The result splits the final knot set into two disjoint vectors:
#' \code{knots} (positions the algorithm itself selected) and
#' \code{fixed_knots} (the user-supplied anchors, echoed back, deduplicated
#' and sorted). The full vector to fit a model on is
#' \code{c(knots, fixed_knots)} (or its sorted form). With the default
#' \code{n_steps = NULL} the loop stops as soon as no further move improves
#' the score, so the end state \emph{is} the score-optimum; with
#' \code{n_steps > 0L} the loop runs that many iterations regardless of
#' improvement, and the score-best point along the trajectory is recoverable
#' from \code{history} via the row at \code{which.min(history$score)}.
#'
#' If the starting state already is (locally) optimal --- typical when
#' \code{fixed_knots} alone overfits in forward mode, or when
#' \code{knot_candidates} is so small that the full pool already minimises
#' the score in backward mode --- the \code{n_steps = NULL} run terminates
#' immediately with \code{history} containing only the initial row and
#' \code{score} equal to the starting score; an explicit \code{n_steps = N}
#' run will instead take \code{N} worsening steps so the post-optimum
#' landscape is still observable.
#'
#' Knots sit next to extrema, not on them. \code{select_knots()} chooses
#' knots for a \emph{good fit}, not for detecting extrema of the response.
#' Knots typically land next to peaks or troughs --- where the curvature is
#' highest --- and not on them. Read local extrema off the \emph{shape} of
#' the fitted curve, not from \code{knots}. A knot supplied
#' through \code{fixed_knots} is the exception: it sits wherever the user
#' puts it, since it is a constraint imposed on the search rather than a
#' finding produced by it.
#'
#' @section When to use this:
#' \code{select_knots()} is for situations where knot \emph{positions} are
#' meaningful and you want a small, inspectable set of interior knots scored
#' under a criterion you choose. If you instead want a continuous roughness
#' penalty over a dense knot grid, the \pkg{mgcv} package's penalised
#' B-splines (\code{bs = "ps"}) or adaptive smoothers (\code{bs = "ad"}) are
#' usually a better fit; for greedy stepwise selection on a hinge-function
#' basis, see the \pkg{earth} package.
#'
#' Like every greedy stepwise procedure, \code{select_knots()} performs a
#' \emph{local} search: in each round it commits to the locally best
#' add/drop. It therefore returns a local optimum of the score, which is
#' usually but not provably the global optimum --- a knot \emph{combination}
#' that beats every individually-best move can be unreachable once an
#' earlier, locally-attractive knot has been picked. The only way to
#' guarantee globality would be exhaustive enumeration over all
#' \eqn{2^{|knot\_candidates|}} subsets, which is exponential and
#' impractical for any realistic candidate pool. In practice the greedy
#' optimum is close; repeating the search in the opposite direction
#' (\code{backward = TRUE}) is a cheap robustness check.
#'
#' @param data Data passed unchanged to \code{score_fun}; typically a
#'   \code{data.frame}. The function itself does not inspect or mutate it.
#' @param knot_candidates Numeric vector of candidate interior-knot positions.
#'   Duplicates and unsorted input are tolerated (sorted/uniq'd internally).
#' @param score_fun \code{\\(data, knots) -> numeric scalar}. Lower is
#'   better. The caller is responsible for fitting whatever model they want
#'   and returning a single criterion value (BIC, AIC, cross-validation
#'   error, ...).
#' @param backward Logical. \code{FALSE} (default) = forward selection,
#'   \code{TRUE} = backward elimination.
#' @param n_steps \code{NULL} (default) or a positive integer.
#'   \code{NULL} runs the greedy search until the next move would not
#'   improve the score and stops there --- so the end state of the search
#'   IS the score optimum (under the current direction and constraints).
#'   A positive integer \code{N} runs exactly \code{N} iterations regardless
#'   of whether each one improves the score; this is a deliberate
#'   exploration mode, intended for inspecting the score landscape just
#'   past a previously-found optimum. The canonical recipe is to first
#'   call with \code{n_steps = NULL} to find the optimum, then re-call with
#'   \code{fixed_knots = prev$knots} and
#'   \code{n_steps = 1L} (or 2L, 3L) to take one or more controlled steps
#'   onward and see how rapidly the score deteriorates. The value
#'   \code{0L} is not allowed (use \code{NULL} for "no steps past the
#'   optimum"); non-integer or negative values error out.
#' @param bulk_gap Integer (>= 0). Minimum index gap between knots removed in
#'   the same round of backward elimination. Ignored when
#'   \code{backward = FALSE}. \code{0L} disables bulk removal. The default
#'   \code{4L} matches the support of a cubic spline basis
#'   (\code{degree + 1} for cubic B-splines; the same value also fits
#'   \code{splines::ns}, which is cubic by construction); for higher-degree
#'   B-splines pick \code{degree + 1} accordingly. Note that
#'   \code{bulk_gap = 1L} imposes no real separation constraint and therefore
#'   removes \emph{every} improving knot in one round -- a very aggressive
#'   setting that defeats the per-round verify-fit's role as a safety net.
#' @param fixed_knots Numeric vector of knot positions that must be present
#'   in every fitted model during the search. In forward mode they are
#'   set from the start, and further knots may be added on top; in
#'   backward mode they are never proposed for removal, neither singly nor
#'   in bulk. Need not be a subset of \code{knot_candidates} -- included
#'   regardless. Duplicates and unsorted input are tolerated. Default
#'   \code{numeric(0)} (no fixed knots).
#' @param verbose Logical. If \code{TRUE}, prints per-round progress
#'   (\code{cat()} to stdout). Default \code{FALSE}.
#' @param n_cores Integer (>= 1). Performance-only knob: the per-round
#'   candidate scoring runs in parallel across \code{n_cores} workers. Default
#'   \code{1L} (serial, no extra dependency). When \code{n_cores > 1L} the
#'   base-R package \pkg{parallel} is loaded and \code{parallel::mclapply()}
#'   is used (fork-based on macOS/Linux; falls back to serial on Windows).
#'   Quick recipes:
#'   \itemize{
#'     \item \code{n_cores = 1L} -- safe default, no extra package loaded.
#'     \item \code{n_cores = max(1L, parallel::detectCores() - 1L)} -- use
#'       all cores except one; good default for interactive multi-tasking.
#'     \item \code{n_cores = parallel::detectCores()} -- use every core;
#'       fastest but resource-hungry (the whole machine is busy).
#'   }
#'   The function errors out with a clear message if \code{n_cores > 1L} but
#'   \pkg{parallel} is not installed.
#'
#'   Reproducibility deserves a note when \code{n_cores > 1L}:
#'   \code{mclapply()} is fork-based and inherits the parent RNG state. If
#'   your \code{score_fun} uses randomness (e.g. cross-validation splits),
#'   set \code{RNGkind("L'Ecuyer-CMRG")} and seed via \code{set.seed()}
#'   \emph{before} calling \code{select_knots()} to get reproducible
#'   per-worker streams; otherwise results can differ run to run and across
#'   \code{n_cores}.
#'
#' @return A list with elements (in this order):
#' \describe{
#'   \item{\code{backward}}{Logical — the direction used.}
#'   \item{\code{knots}}{Sorted numeric vector — the interior knots the
#'     algorithm itself selected, with \code{fixed_knots} \emph{excluded}.
#'     Disjoint from \code{fixed_knots} by construction. The full knot
#'     vector to fit a model on is \code{c(knots, fixed_knots)} (or its
#'     sorted form). With the default \code{n_steps = NULL} this is the
#'     score-optimal selection; with \code{n_steps > 0L} it can be a state
#'     past the optimum.}
#'   \item{\code{fixed_knots}}{Sorted numeric vector — the user-supplied
#'     fixed knots echoed back (deduplicated and sorted). Empty vector if
#'     the caller did not supply any.}
#'   \item{\code{score}}{Numeric — the score at the end state, i.e. at
#'     \code{c(knots, fixed_knots)}.}
#'   \item{\code{n_steps}}{The \code{n_steps} value used (\code{NULL} or a
#'     positive integer).}
#'   \item{\code{history}}{\code{data.frame} of per-round records:
#'     \code{step}, \code{n_knots} (total interior knot count including
#'     \code{fixed_knots}), \code{changed_knot} (the knot added or removed
#'     in that round; \code{NA} for the initial state and for bulk-removal
#'     rounds), \code{score}, \code{extra} (boolean: \code{TRUE} when that
#'     step worsened the score). For \code{n_steps = NULL} runs no
#'     worsening step is taken, so \code{extra} is always \code{FALSE}.}
#' }
#'
#' @examples
#' \dontrun{
#' # Greedy knot selection on a simple synthetic signal.
#' set.seed(1)
#' n <- 200
#' x <- seq(0, 10, length.out = n)
#' y <- stats::rpois(n, lambda = 50 + 30 * sin(x))
#' dat <- data.frame(x = x, y = y)
#'
#' # score_fun: fit a Poisson GLM with a natural cubic spline, return BIC.
#' # Lower is better.
#' fit <- \(d, knots) {
#'     f <- if (length(knots) == 0L) {
#'         y ~ splines::ns(x)
#'     } else {
#'         y ~ splines::ns(x, knots = knots)
#'     }
#'     stats::glm(f, family = stats::poisson(), data = d)
#' }
#' score_bic <- \(d, knots) stats::BIC(fit(d, knots))
#'
#' cand <- seq(1, 9, by = 0.5)
#' res <- select_knots(dat, cand, score_bic, verbose = TRUE)
#' # Full knot vector to fit the final model on:
#' final_knots <- sort(c(res$knots, res$fixed_knots))
#' }
#'
#' @seealso \code{\link[splines]{ns}}, \code{\link[splines]{bs}},
#'   \code{\link[stats]{BIC}}, \code{\link[stats]{AIC}}
#' @export
select_knots <- \(data,
    knot_candidates,
    score_fun,
    backward = FALSE,
    n_steps = NULL,
    bulk_gap = 4L,
    fixed_knots = numeric(0),
    verbose = FALSE,
    n_cores = 1L) {
    stopifnot(
        is.function(score_fun),
        is.logical(backward), length(backward) == 1L, !is.na(backward),
        is.logical(verbose), length(verbose) == 1L, !is.na(verbose),
        is.numeric(bulk_gap), length(bulk_gap) == 1L, bulk_gap >= 0L,
        is.null(n_steps) || (is.numeric(n_steps) &&
            length(n_steps) == 1L &&
            !is.na(n_steps) &&
            n_steps >= 1L &&
            n_steps == as.integer(n_steps)),
        is.numeric(n_cores), length(n_cores) == 1L, n_cores >= 1L,
        is.numeric(fixed_knots), !anyNA(fixed_knots),
        all(is.finite(fixed_knots))
    )

    forward <- !backward
    # !duplicated() + boolean indexing instead of unique(): for atomic
    # vectors unique() strips names, which would silently drop the labels
    # of a named knot_candidates / fixed_knots input. The order() variant
    # of sort() likewise preserves names.
    cand_pool <- knot_candidates[!duplicated(knot_candidates)]
    cand_pool <- cand_pool[order(cand_pool)]
    fixed_knots <- fixed_knots[!duplicated(fixed_knots)]
    fixed_knots <- fixed_knots[order(fixed_knots)]
    dir_lbl <- if (backward) "backward" else "forward"
    verb <- if (forward) "add" else "remove"
    bulk_active <- backward && (bulk_gap >= 1L) # bulk only for backward; may flip off mid-run

    # --- private helpers ---
    # Score for one knot set (Inf if the user's score_fun errors out).
    # Knots are sorted before they reach score_fun so the user's fit code
    # never has to defensively sort.
    score_of <- \(knots) {
        tryCatch(
            suppressWarnings(score_fun(data, sort(knots))),
            error = \(e) Inf
        )
    }
    # Scores for many candidate knot sets. Serial (no extra dep) by default;
    # parallel only when explicitly requested via n_cores > 1L.
    par_score <- \(items, make_knots) {
        if (n_cores > 1L) {
            if (!requireNamespace("parallel", quietly = TRUE)) {
                stop("select_knots: n_cores > 1 requires the 'parallel' package; ",
                    "install it or use n_cores = 1L.",
                    call. = FALSE
                )
            }
            res <- parallel::mclapply(items, \(it) score_of(make_knots(it)),
                mc.cores = n_cores
            )
        } else {
            res <- lapply(items, \(it) score_of(make_knots(it)))
        }
        # a crashed worker returns a try-error object -> treat it as Inf
        vapply(res, \(x) if (is.numeric(x) && length(x) == 1L) x else Inf, numeric(1))
    }

    if (forward) {
        current_knots <- fixed_knots
        # cand_pool[!cand_pool %in% fixed_knots] instead of setdiff(): the
        # latter strips names.
        remaining <- cand_pool[!cand_pool %in% fixed_knots]
    } else {
        combined <- c(cand_pool, fixed_knots)
        combined <- combined[!duplicated(combined)]
        current_knots <- combined[order(combined)]
        remaining <- NULL
    }
    current_score <- score_of(current_knots)

    # best_score is tracked internally to control the bulk-removal gate (bulk
    # is only worthwhile while we are still improving). It is not returned;
    # the returned score is the score at the loop's end state.
    best_score <- current_score

    history <- data.frame(
        step         = 0,
        n_knots      = length(current_knots),
        changed_knot = NA_real_,
        score        = current_score,
        extra        = FALSE
    )

    if (verbose) {
        fixed_lbl <- if (length(fixed_knots) > 0L) {
            sprintf(" (%d fixed)", length(fixed_knots))
        } else {
            ""
        }
        cat(sprintf(
            "Start (%s, %d core%s): %d knots%s,  score = %.2f\n",
            dir_lbl, n_cores, if (n_cores == 1) "" else "s",
            length(current_knots), fixed_lbl, current_score
        ))
    }

    step <- 1

    repeat {
        # -- Score of every one-step neighbour of the current knot set (parallel) --
        # In backward mode, fixed knots are excluded from the removal candidates;
        # `removable_idx` maps from score_candidates positions back into
        # current_knots positions so the bulk-gap distance check (which is in
        # current_knots index space, where the d+1 basis-support argument
        # lives) and the single-knot removal both stay correct.
        if (forward) {
            if (length(remaining) == 0) break
            score_candidates <- par_score(remaining, \(k) c(current_knots, k))
            cand_ids <- remaining
            removable_idx <- NULL
        } else {
            removable_idx <- which(!current_knots %in% fixed_knots)
            if (length(removable_idx) == 0) break
            score_candidates <- par_score(
                removable_idx,
                \(i) current_knots[-i]
            )
            cand_ids <- current_knots[removable_idx]
        }

        # -- Bulk removal: drop multiple well-separated improving knots at once.
        #    Active only when going backward, bulk_gap >= 1, and we are still
        #    below the score minimum. Falls through to the single-knot path if
        #    fewer than 2 non-adjacent improving candidates exist, or if the
        #    verify-fit shows no improvement (in which case bulk mode is
        #    permanently disabled for the rest of the run).
        if (bulk_active && min(score_candidates) < best_score) {
            ord <- order(score_candidates)
            improving <- ord[score_candidates[ord] < best_score]
            # picked stores positions in current_knots (mapped via
            # removable_idx) so the gap check matches the d+1 basis-support
            # argument from the docs.
            picked <- integer(0)
            for (idx in improving) {
                pos <- removable_idx[idx]
                if (length(picked) == 0L || all(abs(pos - picked) >= bulk_gap)) {
                    picked <- c(picked, pos)
                }
            }
            if (length(picked) >= 2L) {
                proposed_knots <- current_knots[-picked]
                verify_score <- score_of(proposed_knots)
                if (verify_score < best_score) {
                    if (verbose) {
                        cat(sprintf(
                            "Step %2d: remove %d knots [bulk]  ->  score %.2f -> %.2f  (delta = %.4f)\n",
                            step, length(picked), current_score, verify_score,
                            verify_score - current_score
                        ))
                    }
                    current_knots <- proposed_knots
                    current_score <- verify_score
                    best_score <- current_score
                    history <- rbind(history, data.frame(
                        step         = step,
                        n_knots      = length(current_knots),
                        changed_knot = NA_real_,
                        score        = current_score,
                        extra        = FALSE
                    ))
                    step <- step + 1
                    if (!is.null(n_steps) && step > n_steps) break
                    next
                } else {
                    if (verbose) {
                        cat(sprintf(
                            paste0(
                                "Step %2d: bulk removal of %d knots did not improve ",
                                "(delta = %.4f) -> switching to single-knot mode\n"
                            ),
                            step, length(picked), verify_score - current_score
                        ))
                    }
                    bulk_active <- FALSE
                    # fall through to single-knot path (uses the same score_candidates)
                }
            }
            # length(picked) < 2 -> fall through to single-knot path; keep bulk_active for next round
        }

        best_idx <- which.min(score_candidates)
        cand_score <- score_candidates[best_idx]
        changed_knot <- cand_ids[best_idx]
        improved <- cand_score < best_score

        # NULL-mode: exit before applying a non-improving step, so the loop
        # end state IS the score optimum. The failed candidate is not recorded
        # in history -- use n_steps = 1L for an explicit one-step look-ahead.
        if (is.null(n_steps) && !improved) {
            if (verbose) {
                cat(sprintf(
                    "Step %2d: best %s would worsen score by %.4f -> stopping at optimum\n",
                    step, verb, cand_score - current_score
                ))
            }
            break
        }

        if (verbose) {
            step_lbl <- if (!is.null(n_steps)) {
                sprintf("  [step %d/%d]", step, n_steps)
            } else {
                ""
            }
            cat(sprintf(
                "Step %2d: %s knot %.4f  ->  score %.2f -> %.2f  (delta = %.4f)%s\n",
                step, verb, changed_knot, current_score, cand_score, cand_score - current_score,
                step_lbl
            ))
        }

        if (forward) {
            current_knots <- c(current_knots, changed_knot)
            remaining <- remaining[-best_idx]
        } else {
            current_knots <- current_knots[-removable_idx[best_idx]]
        }
        current_score <- cand_score

        if (improved) best_score <- cand_score

        history <- rbind(history, data.frame(
            step         = step,
            n_knots      = length(current_knots),
            changed_knot = changed_knot,
            score        = current_score,
            extra        = !improved
        ))

        step <- step + 1

        if (!is.null(n_steps) && step > n_steps) break
    }

    free_knots <- current_knots[!current_knots %in% fixed_knots]

    list(
        backward    = backward,
        knots       = free_knots[order(free_knots)],
        fixed_knots = fixed_knots,
        score       = current_score,
        n_steps     = n_steps,
        history     = history
    )
}
