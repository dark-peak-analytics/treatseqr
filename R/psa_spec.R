# Specification of a probabilistic sensitivity analysis (PSA) for a treatseqr
# model.
#
# The deterministic workflow is:
# \itemize{
#   \item \code{tp_source}
#   \item \code{validate_tp_source()}
#   \item \code{generate_m_list()}
#   \item \code{extrapolate_treatseqr()}
# }
#
# The PSA setp layers on top. One tp_source is made per iteration. Transition
# probabilities can reach it two ways:
#
#   1. survival-powered (the default): parametric survival models per
#      destination, using flexsurv and converted to competing-risk TPs on the
#      hazard scale.
#   2. injection: already-drawn TPs supplied by the modeller to draw_psa()
#
# specify_psa() validates the structure and the survival models once rather than
# every PSA iteration.

#' Specify a probabilistic sensitivity analysis for a treatseqr model
#'
#' Validates the model structure. Returns a specification object for
#' \code{\link{draw_psa}}. Transitions not covered by a survival model are
#' expected to be injected as pre-drawn probabilities via the \code{tp_draws}
#' argument of \code{\link{draw_psa}}.
#'
#' @param surv Named list, one element per state that has at least
#'   one survival-powered transition, in the same structure as a deterministic
#'   \code{tp_source}. Destinations are:
#'   \itemize{
#'     \item a fitted \code{flexsurv::flexsurvreg} object with no covariates.
#'     \item a list with elements \code{dist}, \code{coefs}, and \code{vcov}.
#'   }
#'   May be \code{NULL} if every transition is injected directly.
#' @param m_spec Model specification from \code{\link{specify_m}}.
#' @param state_names Character vector of all state names in the right order,
#' with the absorbing (death) state last.
#' @param th Integer time horizon in model cycles.
#' @param tuota Time unit of the survival analysis, in years.
#' @param cl Cycle length of the model, in years.
#'
#' @return A \code{treatseqr_psa_spec} object. List with elements \code{surv},
#' \code{state_names}, \code{death_name}, \code{n_cycles}, \code{th},
#' \code{tuota}, \code{cl} and \code{stream_names}.
#'
#' @details
#' Clock-reset assumption is used throughout. Cumulative hazards for a
#' transition out of state \eqn{s} are evaluated at time since entry. The
#' survival model for each transition must therefore have been fitted on the
#' same time scale, and must censor for other destinations out of \eqn{s}.
#'
#' Covariates are not supported. You could potentially pre-evaluate the model
#' for the population of interest first using bootstrapping, for instance
#' \code{flexsurv::normboot.flexsurvreg} and then passing in
#' \code{dist}/\code{coefs}/\code{vcov} instead of the \code{flexsurvreg}
#' object.
#'
#' @seealso \code{\link{draw_psa}}, \code{\link{tp_source_psa}}
#' @importFrom utils head tail
#' @importFrom stats setNames
#' @export
#'
#' @examplesIf requireNamespace("flexsurv", quietly = TRUE)
#' th <- 20L
#' m_spec <- specify_m(tunnel_lengths = c(3L, th), pre_tunnel_states = 1L)
#' state_names <- c("on_trt", "surgery", "responder", "die")
#' named_vcov <- function(se, nms, rho = 0) {
#'   m <- outer(se, se) * (diag(length(nms)) * (1 - rho) + rho)
#'   dimnames(m) <- list(nms, nms)
#'   m
#' }
#' spec <- specify_psa(
#'   surv = list(
#'     on_trt = list(
#'       surgery = list(
#'         dist = "exp", coefs = c(rate = log(0.05)),
#'         vcov = named_vcov(0.1, "rate")
#'       ),
#'       die = list(
#'         dist = "exp", coefs = c(rate = log(0.01)),
#'         vcov = named_vcov(0.2, "rate")
#'       )
#'     ),
#'     surgery = list(
#'       responder = list(
#'         dist = "weibull",
#'         coefs = c(shape = log(1.2), scale = log(5)),
#'         vcov = named_vcov(c(0.1, 0.1), c("shape", "scale"), -0.3)
#'       ),
#'       die = list(
#'         dist = "exp", coefs = c(rate = log(0.05)),
#'         vcov = named_vcov(0.1, "rate")
#'       )
#'     ),
#'     responder = list(
#'       die = list(
#'         dist = "exp", coefs = c(rate = log(0.01)),
#'         vcov = named_vcov(0.2, "rate")
#'       )
#'     )
#'   ),
#'   m_spec = m_spec, state_names = state_names, th = th,
#'   tuota = 30.44 / 365.25, cl = 7 / 365.25
#' )
#' spec
specify_psa <- function(surv, m_spec, state_names, th, tuota, cl) {
  assertthat::assert_that(
    is.list(m_spec),
    !is.null(m_spec$pre_tunnels),
    !is.null(m_spec$tunnel_lengths),
    msg = "'m_spec' must be the output of specify_m()"
  )
  n_pre <- m_spec$pre_tunnels
  tunnel_lengths <- m_spec$tunnel_lengths

  assertthat::assert_that(
    is.character(state_names),
    !anyDuplicated(state_names),
    all(nzchar(state_names)),
    msg = "'state_names' must be a character vector of unique, non-empty names"
  )
  assertthat::assert_that(
    length(state_names) == n_pre + length(tunnel_lengths) + 1L,
    msg = sprintf(
      paste0(
        "length(state_names) must equal pre_tunnels + length(tunnel_lengths) + 1 ",
        "(death). Got %d names for %d pre-tunnel state(s) and %d tunnel(s)."
      ),
      length(state_names), n_pre, length(tunnel_lengths)
    )
  )
  assertthat::assert_that(
    is.numeric(th),
    length(th) == 1L,
    th == as.integer(th),
    th >= 2,
    msg = "'th' must be a single integer >= 2"
  )
  assertthat::assert_that(
    is.numeric(tuota), length(tuota) == 1L, tuota > 0,
    is.numeric(cl), length(cl) == 1L, cl > 0,
    msg = "'tuota' and 'cl' must be single positive numbers (years)"
  )

  death_name <- tail(state_names, 1L)
  transient <- head(state_names, -1L)
  n_cycles <- setNames(
    as.integer(c(rep(th, n_pre), tunnel_lengths)),
    transient
  )

  # --- survival models ------------------------------------------------------
  if (is.null(surv)) surv <- list()
  assertthat::assert_that(
    is.list(surv),
    length(surv) == 0L || !is.null(names(surv)),
    msg = "'surv' must be a named list (state -> destination -> model), or NULL"
  )
  if (length(surv) > 0L) {
    bad_states <- setdiff(names(surv), transient)
    assertthat::assert_that(
      length(bad_states) == 0L,
      msg = sprintf(
        "'surv' has states not in state_names (or equal to the death state): %s",
        paste(sQuote(bad_states), collapse = ", ")
      )
    )
    assertthat::assert_that(
      !anyDuplicated(names(surv)),
      msg = "'surv' has duplicated state names"
    )
  }

  surv <- lapply(setNames(names(surv), names(surv)), function(state) {
    dests <- surv[[state]]
    assertthat::assert_that(
      is.list(dests),
      length(dests) >= 1L,
      !is.null(names(dests)),
      all(nzchar(names(dests))),
      !anyDuplicated(names(dests)),
      msg = sprintf(
        "surv$%s must be a non-empty named list of destinations",
        state
      )
    )
    allowed <- setdiff(state_names, state)
    bad_dests <- setdiff(names(dests), allowed)
    assertthat::assert_that(
      length(bad_dests) == 0L,
      msg = sprintf(
        "surv$%s has destination(s) not in state_names (or equal to itself): %s",
        state,
        paste(sQuote(bad_dests), collapse = ", ")
      )
    )
    lapply(setNames(names(dests), names(dests)), function(dest) {
      normalise_surv_spec(dests[[dest]], state, dest)
    })
  })

  stream_names <- unlist(lapply(names(surv), function(state) {
    lapply(names(surv[[state]]), function(dest) {
      paste(state, dest, names(surv[[state]][[dest]]$coefs), sep = ".")
    })
  }), use.names = FALSE)
  if (is.null(stream_names)) stream_names <- character(0)

  structure(
    list(
      surv = surv,
      state_names = state_names,
      death_name = death_name,
      n_cycles = n_cycles,
      th = as.integer(th),
      tuota = tuota,
      cl = cl,
      stream_names = stream_names
    ),
    class = "treatseqr_psa_spec"
  )
}


#' Normalise one survival specification to list(dist, coefs, vcov)
#'
#' Accepts a fitted \code{flexsurvreg} object (no covariates) or a list with
#' \code{dist}, \code{coefs} and \code{vcov}, and returns the latter form with
#' the distribution name canonicalised and everything checked against
#' \code{flexsurv.dists}.
#'
#' @param x The user-supplied specification.
#' @param state,dest Names, used in error messages only.
#' @noRd
normalise_surv_spec <- function(x, state, dest) {
  where <- sprintf("surv$%s$%s", state, dest)

  if (inherits(x, "flexsurvreg")) {
    assertthat::assert_that(
      isTRUE(x$ncovs == 0),
      msg = sprintf(
        paste0(
          "%s: flexsurvreg fits with covariates are not supported. Evaluate ",
          "the model for your population and pass dist/coefs/vcov instead."
        ),
        where
      )
    )
    # res.t always carries parameter names (coef()/vcov() drop them for
    # single-parameter fits). Internal names like "weibull.quiet" map to the
    # public distribution name.
    dist <- sub("\\.quiet$", "", x$dlist$name)
    coefs <- setNames(x$res.t[, "est"], rownames(x$res.t))
    vcov <- as.matrix(stats::vcov(x))
    dimnames(vcov) <- list(names(coefs), names(coefs))
    x <- list(dist = dist, coefs = coefs, vcov = vcov)
  }

  assertthat::assert_that(
    is.list(x),
    all(c("dist", "coefs", "vcov") %in% names(x)),
    msg = sprintf(
      paste0(
        "%s must be a flexsurvreg fit or a list with elements ",
        "'dist', 'coefs' and 'vcov'"
      ),
      where
    )
  )
  extra <- setdiff(names(x), c("dist", "coefs", "vcov"))
  assertthat::assert_that(
    length(extra) == 0L,
    msg = sprintf(
      "%s has unexpected element(s): %s",
      where,
      paste(sQuote(extra), collapse = ", ")
    )
  )

  info <- flexsurv_dist_info(x$dist)
  coefs <- x$coefs
  assertthat::assert_that(
    is.numeric(coefs),
    all(is.finite(coefs)),
    identical(names(coefs), info$pars),
    msg = sprintf(
      "%s: coefs must be a finite numeric vector named exactly c(%s) for dist '%s'",
      where,
      paste(sQuote(info$pars), collapse = ", "),
      info$name
    )
  )
  vcov <- x$vcov
  assertthat::assert_that(
    is.matrix(vcov),
    is.numeric(vcov),
    all(is.finite(vcov)),
    nrow(vcov) == length(coefs),
    ncol(vcov) == length(coefs),
    identical(rownames(vcov), names(coefs)),
    identical(colnames(vcov), names(coefs)),
    msg = sprintf(
      paste0(
        "%s: vcov must be a finite %d x %d numeric matrix with row and column ",
        "names identical to names(coefs)"
      ),
      where, length(coefs), length(coefs)
    )
  )
  assertthat::assert_that(
    isSymmetric(unname(vcov), tol = 1e-8),
    all(diag(vcov) >= 0),
    msg = sprintf("%s: vcov must be symmetric with non-negative diagonal", where)
  )

  list(dist = info$name, coefs = coefs, vcov = vcov)
}


#' @export
print.treatseqr_psa_spec <- function(x, ...) {
  cat("treatseqr PSA specification\n")
  cat(sprintf(
    "  %d transient state(s) + death ('%s'); th = %d cycles\n",
    length(x$n_cycles), x$death_name, x$th
  ))
  cat(sprintf(
    "  tuota = %.4f years, cl = %.4f years (%.3f model cycles per TUoTA)\n",
    x$tuota, x$cl, x$tuota / x$cl
  ))
  cat("  cycles of TPs per state:\n")
  for (s in names(x$n_cycles)) {
    dests <- x$surv[[s]]
    desc <- if (is.null(dests)) {
      "(no survival models: inject via tp_draws)"
    } else {
      paste(
        vapply(names(dests), function(d) {
          sprintf("%s [%s]", d, dests[[d]]$dist)
        }, character(1L)),
        collapse = ", "
      )
    }
    cat(sprintf("    %-14s %4d  %s\n", s, x$n_cycles[[s]], desc))
  }
  cat(sprintf("  %d survival parameter stream(s)\n", length(x$stream_names)))
  invisible(x)
}
