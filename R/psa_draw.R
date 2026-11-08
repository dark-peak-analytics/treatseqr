# Drawing PSA parameters and building per-iteration transition probability
# lists.
#
# draw_psa() performs the RNG for every iteration in one go. Survival parameters
# are drawn on the estimated scale and then mapped to natural. Any injected
# transition probabilities supplied by the modeller are validated and merged
# with the \code{flexsurv}-powered rows. tp_source_psa() then builds the
# tp_source for each iteration, which goes straight into generate_m_list() to
# estimate each probabilistic trace efficiently.

#' Draw PSA parameters for a treatseqr model
#'
#' Draws the survival-model parameters described in a
#' \code{\link{specify_psa}} specification for every PSA iteration, validates
#' any pre-drawn transition probabilities supplied via \code{tp_draws}, and
#' returns everything \code{\link{tp_source_psa}} needs to build the
#' transition probability list for a given iteration.
#'
#' @param spec A \code{treatseqr_psa_spec} from \code{\link{specify_psa}}.
#' @param n_draws Integer. Number of PSA iterations.
#' @param seed Integer seed for \code{\link{draw_uniforms}}. Provide one of
#' \code{seed} or \code{rng} when the specification includes survival models or
#' neither when it does not.
#' @param rng Optional pre-drawn matrix of uniform streams (see
#' \code{\link{draw_uniforms}}) with \code{n_draws} rows and at least the
#' columns named in \code{spec$stream_names}. Supplying it lets you keep one set
#' of streams across scenarios (common random numbers).
#' @param tp_draws Optional list of pre-drawn transition probabilities to
#' inject, in the same list structure as the rest of the package. Each
#' destination is an \code{n_cycles[state] x n_draws} numeric matrix of
#' probabilities (one column per iteration). Use \code{NA} for cycles that
#' should come from the survival models instead. See Details for how rows are
#' merged.
#' @param tp_point Required when \code{tp_draws} is given. Deterministic values
#' of the injected probabilities, in the same list structure with one numeric
#' vector of length \code{n_cycles[state]} per destination. \code{NA} must
#' appear in exactly the same cycles as in \code{tp_draws}.
#'
#' @return An object of class \code{treatseqr_psa_draws}: a list with
#' \code{spec}, \code{n_draws}, \code{rng}, \code{params}, \code{tp_draws} and
#' \code{tp_point} and \code{injected_rows}.
#'
#' @details
#' Within one row, a state's exit probabilities form a single probability
#' vector. Hazards and injected probabilities cannot be mixed within it. For
#' each state and cycle, if \emph{any} destination in \code{tp_draws} has a
#' non-\code{NA} value, the whole row is taken from \code{tp_draws} and any
#' destination left as \code{NA} in that row becomes 0. Otherwise the row comes
#' from the survival models. Every cycle of every transient state must be
#' covered by one or the other.
#'
#' For holding states, a fixed-duration tunnel whose final cycle is
#' survival-powered re-applies that cycle's hazard increment indefinitely to
#' anyone still in the tunnel (see \code{vignette("fixed-duration-tunnels",
#' package = "treatseqr")}). If that increment comes from a non-exponential
#' distribution a warning is issued: the frozen increment is unlikely to be what
#' you intend. Either use \code{dist = "exp"} for those transitions, or inject a
#' decision node for the final cycle via \code{tp_draws}: a Dirichlet row from
#' \code{\link{dirichlet_draws}} with \code{NA} in every earlier cycle, as
#' shown in \code{vignette("probabilistic", package = "treatseqr")}.
#'
#' @seealso \code{\link{specify_psa}}, \code{\link{tp_source_psa}},
#'   \code{\link{draw_uniforms}}
#' @export
#'
#' @examplesIf requireNamespace("flexsurv", quietly = TRUE)
#' th <- 20L
#' m_spec <- specify_m(tunnel_lengths = c(3L, th), pre_tunnel_states = 1L)
#' state_names <- c("on_trt", "surgery", "responder", "die")
#' nv <- function(se, nm) matrix(se^2, dimnames = list(nm, nm))
#' spec <- specify_psa(
#'   surv = list(
#'     on_trt = list(
#'       surgery = list(dist = "exp", coefs = c(rate = log(0.05)), vcov = nv(0.1, "rate")),
#'       die = list(dist = "exp", coefs = c(rate = log(0.01)), vcov = nv(0.2, "rate"))
#'     ),
#'     surgery = list(
#'       die = list(dist = "exp", coefs = c(rate = log(0.05)), vcov = nv(0.1, "rate"))
#'     ),
#'     responder = list(
#'       die = list(dist = "exp", coefs = c(rate = log(0.01)), vcov = nv(0.2, "rate"))
#'     )
#'   ),
#'   m_spec, state_names, th, tuota = 30.44 / 365.25, cl = 7 / 365.25
#' )
#' # a decision node at the end of the 3-cycle surgery tunnel: 90% respond,
#' # 10% die, nobody stays. Cycles 1-2 (NA) come from the survival models.
#' n <- 100
#' node <- dirichlet_draws(draw_uniforms(n, c("responder", "die"), seed = 2),
#'                         c(responder = 90, die = 10))
#' tp_draws <- list(surgery = list(
#'   responder = rbind(NA, NA, node[, "responder"]),
#'   die       = rbind(NA, NA, node[, "die"])
#' ))
#' tp_point <- list(surgery = list(responder = c(NA, NA, 0.9), die = c(NA, NA, 0.1)))
#' draws <- draw_psa(spec, n, seed = 1, tp_draws = tp_draws, tp_point = tp_point)
#' draws
#' tp_1 <- tp_source_psa(draws, iter = 1)
#' round(1 - Reduce(`+`, tp_1$surgery), 3) # p_stay in surgery: 0 at cycle 3
draw_psa <- function(
  spec,
  n_draws,
  seed = NULL,
  rng = NULL,
  tp_draws = NULL,
  tp_point = NULL
) {
  assertthat::assert_that(
    inherits(spec, "treatseqr_psa_spec"),
    msg = "'spec' must be the output of specify_psa()"
  )
  assertthat::assert_that(
    is.numeric(n_draws),
    length(n_draws) == 1L,
    n_draws >= 1,
    n_draws == as.integer(n_draws),
    msg = "'n_draws' must be a single positive integer"
  )
  n_draws <- as.integer(n_draws)
  transient <- head(spec$state_names, -1L)
  n_cycles <- spec$n_cycles

  # --- uniform streams and survival parameter draws -------------------------
  has_surv <- length(spec$stream_names) > 0L
  if (has_surv) {
    assertthat::assert_that(
      xor(is.null(seed), is.null(rng)),
      msg = "Provide exactly one of 'seed' or 'rng'"
    )
    if (is.null(rng)) {
      rng <- draw_uniforms(n_draws, spec$stream_names, seed)
    } else {
      missing_streams <- setdiff(spec$stream_names, colnames(rng))
      assertthat::assert_that(
        is.matrix(rng),
        is.numeric(rng),
        nrow(rng) == n_draws,
        length(missing_streams) == 0L,
        msg = sprintf(
          paste0(
            "'rng' must be a numeric matrix with n_draws rows containing every ",
            "stream in spec$stream_names. Missing: %s"
          ),
          paste(sQuote(missing_streams), collapse = ", ")
        )
      )
    }
  } else {
    assertthat::assert_that(
      is.null(seed) && is.null(rng),
      msg = "The specification has no survival models, so 'seed'/'rng' are not used"
    )
    rng <- matrix(numeric(0), nrow = n_draws, ncol = 0L)
  }

  params <- lapply(setNames(names(spec$surv), names(spec$surv)), function(state) {
    dests <- spec$surv[[state]]
    lapply(setNames(names(dests), names(dests)), function(dest) {
      d <- dests[[dest]]
      cols <- paste(state, dest, names(d$coefs), sep = ".")
      est_scale <- mvnorm_draws(rng[, cols, drop = FALSE], d$coefs, vcov = d$vcov)
      # n_parameters x n_draws: one column per iteration (R is column-major)
      list(dist = d$dist, params = t(flexsurv_inv_transform(est_scale, d$dist)))
    })
  })

  # --- injected transition probabilities ------------------------------------
  injected_rows <- lapply(n_cycles, function(n) rep(FALSE, n))

  if (!is.null(tp_draws)) {
    assertthat::assert_that(
      !is.null(tp_point),
      msg = "'tp_point' (the deterministic values of the injected probabilities) is required when 'tp_draws' is given"
    )
    tp_draws <- validate_injected(
      tp_draws, spec, n_draws, what = "tp_draws", is_matrix = TRUE
    )
    tp_point <- validate_injected(
      tp_point, spec, n_draws, what = "tp_point", is_matrix = FALSE
    )
    assertthat::assert_that(
      identical(sort(names(tp_point)), sort(names(tp_draws))),
      msg = "'tp_point' and 'tp_draws' must cover the same states"
    )

    for (state in names(tp_draws)) {
      dests_d <- tp_draws[[state]]
      dests_p <- tp_point[[state]]
      assertthat::assert_that(
        identical(sort(names(dests_d)), sort(names(dests_p))),
        msg = sprintf(
          "tp_point$%s and tp_draws$%s must have the same destinations",
          state, state
        )
      )
      dests_p <- dests_p[names(dests_d)]

      # merge rule: a row is injected if ANY destination is non-NA in it
      row_injected <- Reduce(`|`, lapply(dests_d, function(m) {
        rowSums(!is.na(m)) > 0
      }))
      # within an injected row, every draw of every destination must be given
      # (or NA -> 0); partially-NA rows across draws would be ambiguous
      for (dest in names(dests_d)) {
        m <- dests_d[[dest]]
        partial <- row_injected & (rowSums(is.na(m)) > 0) & (rowSums(is.na(m)) < n_draws)
        assertthat::assert_that(
          !any(partial),
          msg = sprintf(
            paste0(
              "tp_draws$%s$%s: cycle(s) %s mix NA and non-NA values across ",
              "draws. Each cycle must be entirely NA or entirely given."
            ),
            state, dest, paste(which(partial), collapse = ", ")
          )
        )
        m[row_injected & is.na(m)] <- 0
        dests_d[[dest]] <- m

        v <- dests_p[[dest]]
        point_marks <- !is.na(v)
        assertthat::assert_that(
          !any(point_marks & !row_injected),
          msg = sprintf(
            paste0(
              "tp_point$%s$%s has values in cycle(s) %s, which are NA for every ",
              "destination in tp_draws$%s. NA patterns must agree at the row level."
            ),
            state, dest,
            paste(which(point_marks & !row_injected), collapse = ", "),
            state
          )
        )
        v[row_injected & is.na(v)] <- 0
        dests_p[[dest]] <- v
      }

      # row sums <= 1 for every draw and for the point estimates
      row_tot <- Reduce(`+`, dests_d) # n_cycles x n_draws
      bad <- which(row_injected & (apply(row_tot, 1, max) > 1 + 1e-10))
      assertthat::assert_that(
        length(bad) == 0L,
        msg = sprintf(
          "tp_draws$%s: probabilities sum to > 1 in cycle(s) %s for at least one draw",
          state, paste(bad, collapse = ", ")
        )
      )
      point_tot <- Reduce(`+`, dests_p)
      bad <- which(row_injected & point_tot > 1 + 1e-10)
      assertthat::assert_that(
        length(bad) == 0L,
        msg = sprintf(
          "tp_point$%s: probabilities sum to > 1 in cycle(s) %s",
          state, paste(bad, collapse = ", ")
        )
      )

      tp_draws[[state]] <- dests_d
      tp_point[[state]] <- dests_p
      injected_rows[[state]] <- row_injected

      message(sprintf(
        "draw_psa: state '%s' cycle(s) %s injected from tp_draws",
        state, format_cycles(which(row_injected))
      ))
    }
  }

  # --- coverage: every cycle of every transient state has a source ----------
  for (state in transient) {
    surv_rows <- state %in% names(spec$surv)
    uncovered <- which(!injected_rows[[state]] & !surv_rows)
    assertthat::assert_that(
      length(uncovered) == 0L,
      msg = sprintf(
        paste0(
          "State '%s' has no source of transition probabilities for cycle(s) %s: ",
          "add survival models in specify_psa() or inject them via tp_draws"
        ),
        state, format_cycles(uncovered)
      )
    )
  }

  # --- holding-state warning ------------------------------------------------
  for (state in transient) {
    n <- n_cycles[[state]]
    is_fixed_tunnel <- n < spec$th
    last_is_surv <- !injected_rows[[state]][n] && state %in% names(spec$surv)
    if (is_fixed_tunnel && last_is_surv) {
      dists <- vapply(spec$surv[[state]], `[[`, character(1L), "dist")
      non_exp <- names(dists)[dists != "exp"]
      if (length(non_exp) > 0L) {
        warning(sprintf(
          paste0(
            "State '%s' is a fixed-duration tunnel (%d cycles) whose final cycle ",
            "is a holding state powered by non-exponential survival model(s) ",
            "(%s). The cycle-%d hazard increment will be re-applied indefinitely ",
            "to anyone remaining. Use dist = \"exp\" for these transitions, or ",
            "inject a decision node for cycle %d via tp_draws (a dirichlet_draws() ",
            "row with NA in earlier cycles; see vignette(\"probabilistic\"))."
          ),
          state, n,
          paste(sprintf("%s [%s]", non_exp, dists[non_exp]), collapse = ", "),
          n, n
        ), call. = FALSE)
      }
    }
  }

  structure(
    list(
      spec = spec,
      n_draws = n_draws,
      rng = rng,
      params = params,
      tp_draws = tp_draws,
      tp_point = tp_point,
      injected_rows = injected_rows
    ),
    class = "treatseqr_psa_draws"
  )
}


#' Validate the structure of tp_draws / tp_point
#' @noRd
validate_injected <- function(x, spec, n_draws, what, is_matrix) {
  transient <- head(spec$state_names, -1L)
  assertthat::assert_that(
    is.list(x),
    length(x) >= 1L,
    !is.null(names(x)),
    all(nzchar(names(x))),
    !anyDuplicated(names(x)),
    msg = sprintf("'%s' must be a non-empty named list (state -> destination)", what)
  )
  bad_states <- setdiff(names(x), transient)
  assertthat::assert_that(
    length(bad_states) == 0L,
    msg = sprintf(
      "'%s' has states not in state_names (or equal to the death state): %s",
      what, paste(sQuote(bad_states), collapse = ", ")
    )
  )
  lapply(setNames(names(x), names(x)), function(state) {
    dests <- x[[state]]
    n <- spec$n_cycles[[state]]
    assertthat::assert_that(
      is.list(dests),
      length(dests) >= 1L,
      !is.null(names(dests)),
      all(nzchar(names(dests))),
      !anyDuplicated(names(dests)),
      msg = sprintf("%s$%s must be a non-empty named list of destinations", what, state)
    )
    bad_dests <- setdiff(names(dests), setdiff(spec$state_names, state))
    assertthat::assert_that(
      length(bad_dests) == 0L,
      msg = sprintf(
        "%s$%s has destination(s) not in state_names (or equal to itself): %s",
        what, state, paste(sQuote(bad_dests), collapse = ", ")
      )
    )
    lapply(setNames(names(dests), names(dests)), function(dest) {
      v <- dests[[dest]]
      where <- sprintf("%s$%s$%s", what, state, dest)
      if (is_matrix) {
        assertthat::assert_that(
          is.matrix(v),
          is.numeric(v),
          nrow(v) == n,
          ncol(v) == n_draws,
          msg = sprintf(
            "%s must be a numeric %d x %d matrix (n_cycles x n_draws); got %s",
            where, n, n_draws,
            if (is.matrix(v)) paste(dim(v), collapse = " x ") else class(v)[1L]
          )
        )
      } else {
        assertthat::assert_that(
          is.numeric(v),
          is.null(dim(v)),
          length(v) == n,
          msg = sprintf(
            "%s must be a numeric vector of length %d (n_cycles); got length %d",
            where, n, length(v)
          )
        )
      }
      ok <- is.na(v) | (v >= 0 & v <= 1)
      assertthat::assert_that(
        all(ok),
        msg = sprintf("%s has values outside [0, 1]", where)
      )
      v
    })
  })
}


#' Compact "1-2, 5" style formatting of cycle indices
#' @noRd
format_cycles <- function(idx) {
  if (length(idx) == 0L) return("")
  runs <- split(idx, cumsum(c(1L, diff(idx) != 1L)))
  paste(vapply(runs, function(r) {
    if (length(r) == 1L) as.character(r) else sprintf("%d-%d", r[1L], r[length(r)])
  }, character(1L)), collapse = ", ")
}


#' Build the transition probability list for one PSA iteration
#'
#' Combines the survival-model parameters drawn by \code{\link{draw_psa}} for
#' iteration \code{iter} (converted to competing-risk transition probabilities
#' on the hazard scale) with any injected probabilities for that iteration,
#' and returns a \code{tp_source} list ready for \code{\link{generate_m_list}}.
#' With \code{iter = NULL} the base case is returned: survival models at their
#' point estimates and injected rows from \code{tp_point}.
#'
#' @param draws A \code{treatseqr_psa_draws} object from \code{\link{draw_psa}}.
#' @param iter Integer in \code{1:n_draws}, or \code{NULL} for the base case.
#'
#' @return A named list in the package's \code{state -> destination -> vector}
#'   structure (transient states in \code{state_names} order, destinations in
#'   canonical order among those present), suitable for
#'   \code{\link{validate_tp_source}} / \code{\link{generate_m_list}}.
#'
#' @details
#' Survival models are evaluated on the state clock: for a state with
#' \eqn{n} cycles of transition probabilities, the cumulative hazards are
#' evaluated at \eqn{t = 0, 1, \dots, n} cycles since entry, converted to the
#' survival model's time unit with \code{cl / tuota}, and differenced per cycle
#' (see \code{\link{cum_haz_to_tps}}).
#'
#' @seealso \code{\link{draw_psa}}
#' @export
#'
#' @examplesIf requireNamespace("flexsurv", quietly = TRUE)
#' # see ?draw_psa for a full example
tp_source_psa <- function(draws, iter = NULL) {
  assertthat::assert_that(
    inherits(draws, "treatseqr_psa_draws"),
    msg = "'draws' must be the output of draw_psa()"
  )
  spec <- draws$spec
  base_case <- is.null(iter)
  if (!base_case) {
    assertthat::assert_that(
      is.numeric(iter),
      length(iter) == 1L,
      iter == as.integer(iter),
      iter >= 1,
      iter <= draws$n_draws,
      msg = sprintf("'iter' must be a single integer in 1:%d, or NULL for the base case", draws$n_draws)
    )
    iter <- as.integer(iter)
  }

  transient <- head(spec$state_names, -1L)
  out <- lapply(setNames(transient, transient), function(state) {
    n <- spec$n_cycles[[state]]
    injected <- draws$injected_rows[[state]]

    # survival-powered rows
    tp_surv <- NULL
    if (state %in% names(spec$surv)) {
      times <- seq.int(0L, n) * spec$cl / spec$tuota
      H_list <- lapply(spec$surv[[state]], function(d) NULL)
      for (dest in names(spec$surv[[state]])) {
        d <- spec$surv[[state]][[dest]]
        pars <- if (base_case) {
          coef_mat <- matrix(d$coefs, nrow = 1L, dimnames = list(NULL, names(d$coefs)))
          drop(flexsurv_inv_transform(coef_mat, d$dist))
        } else {
          draws$params[[state]][[dest]]$params[, iter]
        }
        H_list[[dest]] <- do.call(
          flexsurv_cum_haz(d$dist),
          c(list(times), as.list(pars))
        )
      }
      tp_surv <- cum_haz_to_tps(H_list)
    }

    # injected rows
    inj <- if (base_case) draws$tp_point[[state]] else {
      lapply(draws$tp_draws[[state]], function(m) m[, iter])
    }

    dests <- setdiff(spec$state_names, state)
    dests <- dests[dests %in% union(colnames(tp_surv), names(inj))]
    lapply(setNames(dests, dests), function(dest) {
      v <- numeric(n)
      if (!is.null(tp_surv) && dest %in% colnames(tp_surv)) {
        v[!injected] <- tp_surv[!injected, dest]
      }
      if (!is.null(inj) && dest %in% names(inj)) {
        v[injected] <- inj[[dest]][injected]
      }
      v
    })
  })
  out
}


#' @export
print.treatseqr_psa_draws <- function(x, ...) {
  cat("treatseqr PSA draws\n")
  cat(sprintf(
    "  %d iteration(s); %d uniform stream(s) used\n",
    x$n_draws, ncol(x$rng)
  ))
  for (state in names(x$spec$n_cycles)) {
    n <- x$spec$n_cycles[[state]]
    inj <- which(x$injected_rows[[state]])
    src <- if (length(inj) == 0L) {
      "survival"
    } else if (length(inj) == n) {
      "injected"
    } else {
      sprintf("survival; injected at cycle(s) %s", format_cycles(inj))
    }
    cat(sprintf("    %-14s %4d cycles  %s\n", state, n, src))
  }
  invisible(x)
}
