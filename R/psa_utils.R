# Low-level utilities for probabilistic sensitivity analysis (PSA).
#
# These are the building blocks that `specify_psa()`, `draw_psa()` and the
# `draw_tp_*()` helpers are built from. They are exported because modellers
# who draw their own transition probabilities (and inject them via the
# `tp_draws` argument of `draw_psa()`) will usually want the same pieces:
# named uniform streams, multivariate-normal and Dirichlet draws, and the
# competing-risks conversion from cumulative hazards to transition
# probabilities.
#
# See `vignette("probabilistic", package = "treatseqr")` for the reasoning
# behind each of these.

#' Draw named uniform random-number streams without disturbing the global RNG
#'
#' Draws an \code{n_draws x length(stream_names)} matrix of \code{runif()}
#' values, one named column ("stream") per parameter that will be sampled in
#' the PSA. The draw happens inside a local RNG scope: \code{.Random.seed} is
#' restored on exit, so \code{set.seed()} calls elsewhere in your model behave
#' identically whether or not this function has run.
#'
#' Naming the streams (for example \code{"surgery.die.rate"}) and looking them
#' up by name rather than by column position means that adding or re-ordering
#' a transition later does not shift the random numbers feeding every other
#' transition, which preserves common random numbers between scenarios.
#'
#' @param n_draws Integer. Number of PSA iterations (rows).
#' @param stream_names Character. rng stream id (origin.destination.parameter).
#' @param seed Integer seed passed to \code{set.seed()} inside the local scope.
#'
#' @return A numeric matrix with \code{n_draws} rows and
#'   \code{length(stream_names)} named columns of values in (0, 1).
#' @importFrom stats runif
#' @export
#'
#' @examples
#' rng <- draw_uniforms(
#'   5,
#'   c("on_trt.die.rate", "on_trt.surgery.rate"),
#'   seed = 1
#' )
#' dim(rng)
#' colnames(rng)
draw_uniforms <- function(n_draws, stream_names, seed) {
  assertthat::assert_that(
    is.numeric(n_draws),
    length(n_draws) == 1L,
    n_draws >= 1,
    msg = "'n_draws' must be a single positive integer"
  )
  assertthat::assert_that(
    is.character(stream_names),
    length(stream_names) >= 1L,
    !anyDuplicated(stream_names),
    all(nzchar(stream_names)),
    msg = paste0(
      "'stream_names' must be a non-empty character vector of unique,",
      "  non-empty names"
    )
  )
  assertthat::assert_that(
    is.numeric(seed),
    length(seed) == 1L,
    msg = "'seed' must be a single number"
  )

  # Save and restore the global RNG state so this draw is side-effect free.
  has_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (has_seed) get(".Random.seed", envir = globalenv())
  on.exit({
    if (has_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else {
      rm(".Random.seed", envir = globalenv())
    }
  })

  set.seed(seed)
  matrix(
    runif(n_draws * length(stream_names)),
    nrow = n_draws,
    dimnames = list(NULL, stream_names)
  )
}


#' Multivariate normal draws from uniform streams (Cholesky method)
#'
#' Inverse-transform sampling: \code{qnorm()} turns the uniform streams into
#' independent standard normals, which are then correlated with the Cholesky
#' factor of \code{vcov} and shifted by \code{params}. This is what
#' \code{mvtnorm::rmvnorm()} does internally, but driven by pre-drawn uniforms
#' so that the random-number stream is under your control.
#'
#' Provide exactly one of \code{vcov} or \code{R}. Pre-computing
#' \code{R <- chol(vcov, pivot = TRUE)} (and un-pivoting it) once per PSA saves
#' time when the same covariance is used for many draws.
#'
#' @param rng Numeric matrix of uniform draws with one column per parameter and
#'   one row per PSA iteration (see \code{\link{draw_uniforms}}).
#' @param params Named numeric vector of means (for example regression
#'   coefficients on the estimation scale).
#' @param vcov Variance-covariance matrix for \code{params}. Its row and column
#'   names must be identical to \code{names(params)} so that nothing can be
#'   silently misaligned.
#' @param R Optional pre-computed, un-pivoted upper Cholesky factor of
#'   \code{vcov} (\code{ncol(R) == length(params)}).
#'
#' @return Numeric matrix with \code{nrow(rng)} rows and one named column per
#'   parameter.
#' @importFrom stats qnorm
#' @export
#'
#' @examples
#' rng <- draw_uniforms(1000, c("shape", "scale"), seed = 1)
#' vc <- matrix(c(0.01, 0.003, 0.003, 0.01), 2, dimnames = list(
#'   c("shape", "scale"), c("shape", "scale")
#' ))
#' draws <- mvnorm_draws(rng, c(shape = log(1.5), scale = log(3)), vcov = vc)
#' round(colMeans(draws), 2)
mvnorm_draws <- function(rng, params, vcov = NULL, R = NULL) {
  assertthat::assert_that(
    xor(is.null(vcov), is.null(R)),
    msg = "Provide exactly one of 'vcov' or 'R'"
  )
  assertthat::assert_that(
    is.numeric(params),
    !is.null(names(params)),
    all(nzchar(names(params))),
    msg = "'params' must be a named numeric vector"
  )
  assertthat::assert_that(
    is.matrix(rng),
    is.numeric(rng),
    nrow(rng) >= 1L,
    ncol(rng) == length(params),
    msg = "'rng' must be a numeric matrix with one column per element of 'params'"
  )

  if (is.null(R)) {
    # NOTE: `colnames(vcov) == names(params)` would be logical(0) when vcov has
    # no dimnames, and assert_that(logical(0)) passes. identical() does not.
    assertthat::assert_that(
      is.matrix(vcov),
      is.numeric(vcov),
      identical(colnames(vcov), names(params)),
      identical(rownames(vcov), names(params)),
      msg = paste0(
        "'vcov' must be a numeric matrix whose row and column names are ",
        "identical to names(params)"
      )
    )
    # pivoting copes with near-singular matrices (very high correlation)
    R <- chol(vcov, pivot = TRUE)
    R <- R[, order(attr(R, "pivot")), drop = FALSE]
  }

  assertthat::assert_that(
    is.matrix(R),
    ncol(R) == length(params),
    msg = "'R' must be a matrix with ncol(R) == length(params)"
  )

  draws <- qnorm(rng) %*% R
  out <- sweep(draws, 2, params, "+")
  colnames(out) <- names(params)
  out
}


#' Dirichlet draws from uniform streams
#'
#' A Dirichlet(\eqn{\alpha}) draw is a set of independent
#' Gamma(\eqn{\alpha_k}, 1) draws normalised to sum to 1. Using
#' \code{qgamma()} on pre-drawn uniform streams keeps the draw reproducible and
#' under the control of the named random-number streams used throughout the
#' PSA. Each row of the result is one probability vector summing to 1, so it
#' is the natural way to draw a decision-tree split or one row of a transition
#' probability matrix from pseudo-counts.
#'
#' @param rng Numeric matrix of uniform draws with one column per category and
#'   one row per PSA iteration.
#' @param alpha Positive numeric vector of Dirichlet parameters
#'   (pseudo-counts), one per column of \code{rng}. If named, the names are
#'   used as column names of the result.
#'
#' @return Numeric matrix with \code{nrow(rng)} rows, each summing to 1, and
#'   \code{length(alpha)} columns. Column means are approximately
#'   \code{alpha / sum(alpha)}.
#' @importFrom stats qgamma
#' @export
#'
#' @examples
#' rng <- draw_uniforms(1000, c("a", "b", "c"), seed = 2)
#' d <- dirichlet_draws(rng, c(a = 35, b = 55, c = 10))
#' round(colMeans(d), 2)
#' all(abs(rowSums(d) - 1) < 1e-12)
dirichlet_draws <- function(rng, alpha) {
  assertthat::assert_that(
    is.numeric(alpha),
    length(alpha) >= 1L,
    all(alpha > 0),
    msg = "'alpha' must be a numeric vector of positive values"
  )
  assertthat::assert_that(
    is.matrix(rng),
    is.numeric(rng),
    nrow(rng) >= 1L,
    ncol(rng) == length(alpha),
    msg = "'rng' must be a numeric matrix with one column per element of 'alpha'"
  )

  g <- qgamma(rng, shape = rep(alpha, each = nrow(rng)), rate = 1)
  out <- g / rowSums(g)
  colnames(out) <- names(alpha)
  out
}


#' Convert cumulative hazards to competing-risk transition probabilities
#'
#' Takes the cumulative hazard \eqn{H_k(t)} for each destination \eqn{k},
#' evaluated at the cycle boundaries \eqn{t = 0, 1, \dots, n}, and returns the
#' per-cycle probability of moving to each destination.
#'
#' For each cycle the \eqn{\Delta H_k} are summed across destinations to give
#' the total hazard. Exit probability is \eqn{1 - \exp(-\sum_k \Delta H_k)}, and
#' that is shared between destinations in proportion to their increments.
#' Converting each destination via \eqn{1 - \exp(-\Delta H_k)} and summing is
#' wrong. Conversion is exact when the hazards are proportional within a cycle,
#' and a very close approximation otherwise.
#'
#' \code{p_stay} (\eqn{1 - p_{exit}}) is not returned. \code{treatseqr} computes
#' it as the complement when it builds the transition matrices.
#'
#' @param H_list Named list with one element per destination. Each element is a
#' numeric vector of cumulative hazards of length \code{n + 1}, evaluated each
#' cycle starting from \eqn{t = 0} (where \eqn{H = 0}).
#'
#' @return Numeric \code{n x K} matrix of transition probabilities: one row per
#'   cycle, one named column per destination. Row sums are exit probabilities.
#' @export
#'
#' @examples
#' # exponential hazard, rate 0.1 per time unit, 5 cycles of length 1
#' H <- list(die = 0.1 * 0:5, progress = 0.05 * 0:5)
#' p <- cum_haz_to_tps(H)
#' p
#' # total exit probability each cycle is 1 - exp(-0.15)
#' all.equal(rowSums(p), rep(1 - exp(-0.15), 5))
cum_haz_to_tps <- function(H_list) {
  assertthat::assert_that(
    is.list(H_list),
    length(H_list) >= 1L,
    !is.null(names(H_list)),
    all(nzchar(names(H_list))),
    msg = "'H_list' must be a non-empty named list of cumulative hazard vectors"
  )
  lens <- lengths(H_list)
  assertthat::assert_that(
    all(lens == lens[1L]),
    lens[1L] >= 2L,
    msg = "All elements of 'H_list' must be numeric vectors of the same length (>= 2)"
  )
  assertthat::assert_that(
    all(vapply(H_list, is.numeric, logical(1L))),
    all(vapply(H_list, function(h) all(is.finite(h)), logical(1L))),
    all(vapply(H_list, function(h) all(diff(h) >= 0), logical(1L))),
    msg = "Cumulative hazards must be finite and non-decreasing"
  )

  # (n + 1) x K, keeps destination names
  H_mat <- do.call(cbind, H_list)
  dH <- diff(H_mat)
  total <- rowSums(dH)
  p_exit <- 1 - exp(-total)

  # Share p_exit between destinations in proportion to dH.
  dH * ifelse(total > 0, p_exit / total, 0)
}


# ---------------------------------------------------------------------------
# flexsurv interface
# ---------------------------------------------------------------------------

#' Check that flexsurv is available
#' @noRd
check_flexsurv <- function() {
  assertthat::assert_that(
    requireNamespace("flexsurv", quietly = TRUE),
    msg = paste0(
      "Package 'flexsurv' is required for survival-powered PSA draws. ",
      "Install it with install.packages(\"flexsurv\")."
    )
  )
  invisible(TRUE)
}


#' Look up a distribution in flexsurv.dists (case-insensitive)
#'
#' Returns the \code{flexsurv::flexsurv.dists} entry for \code{dist}, which
#' records the parameter names (in order), the transforms used for estimation
#' and their inverses. \code{"exponential"} is accepted as an alias for
#' \code{"exp"}.
#'
#' @param dist Character scalar distribution name as used by
#'   \code{flexsurv::flexsurvreg()} (e.g. \code{"weibull"}, \code{"weibullPH"},
#'   \code{"gompertz"}, \code{"lnorm"}, \code{"llogis"}, \code{"gamma"},
#'   \code{"gengamma"}, \code{"genf"}, \code{"exp"}).
#' @return A list with elements \code{name}, \code{pars}, \code{transforms} and
#'   \code{inv.transforms}.
#' @noRd
flexsurv_dist_info <- function(dist) {
  check_flexsurv()
  assertthat::assert_that(
    is.character(dist),
    length(dist) == 1L,
    nzchar(dist),
    msg = "'dist' must be a single distribution name"
  )
  dists <- flexsurv::flexsurv.dists
  key <- tolower(dist)
  if (key == "exponential") key <- "exp"
  idx <- match(key, tolower(names(dists)))
  assertthat::assert_that(
    !is.na(idx),
    msg = sprintf(
      "Unsupported distribution '%s'. Must be one of: %s",
      dist,
      paste(names(dists), collapse = ", ")
    )
  )
  info <- dists[[idx]]
  info$name <- names(dists)[idx]
  info
}


#' Map flexsurv parameter draws back to the natural scale
#'
#' \code{flexsurv} produces transformed vcovs (for example \code{log(rate)},
#' \code{log(shape)}, \code{log(scale)}), so mvnorm draws reflect that. This
#' function applies the distribution's inverse transforms column by column,
#' which then respects the scale of the estimation from \code{flexsurv}. when
#' providing \code{flexsurvreg} results that you have pulled out yourself, be
#' very careful that you provide them as they come out of the package, and don't
#' tranform them yourself. It is simpler to use the \code{flexsurvreg} objects
#' themselves if you are in doubt, which will bypass the use of this function.
#'
#' @param draws Numeric matrix of parameter draws, one row per PSA iteration.
#' Column names must equal the distribution's parameter names in
#' \code{flexsurv}'s order (as returned by \code{names(coef(fit))}).
#' @param dist Distribution name (see \code{\link{flexsurv_cum_haz}}).
#'
#' @return \code{draws} with each column mapped to the natural scale.
#' @export
#'
#' @examplesIf requireNamespace("flexsurv", quietly = TRUE)
#' fit <- flexsurv::flexsurvreg(
#'   survival::Surv(recyrs, censrec) ~ 1,
#'   data = flexsurv::bc, dist = "weibull"
#' )
#' est <- flexsurv_inv_transform(t(coef(fit)), "weibull")
#' all.equal(as.numeric(est), unname(fit$res[, "est"]))
flexsurv_inv_transform <- function(draws, dist) {
  info <- flexsurv_dist_info(dist)
  assertthat::assert_that(
    is.matrix(draws),
    is.numeric(draws),
    identical(colnames(draws), info$pars),
    msg = sprintf(
      "colnames(draws) must be exactly c(%s) for dist '%s'",
      paste(sQuote(info$pars), collapse = ", "),
      info$name
    )
  )
  for (j in seq_along(info$pars)) {
    draws[, j] <- info$inv.transforms[[j]](draws[, j])
  }
  draws
}


#' Get the cumulative hazard function for a flexsurv distribution
#'
#' Returns the \code{flexsurv::H<dist>()} function (for example
#' \code{flexsurv::Hweibull}) for use in building competing-risk transition
#' probabilities. The returned function takes a vector of times as its first
#' argument followed by the distribution's natural-scale parameters by name.
#'
#' @param dist Distribution name as used by \code{flexsurv::flexsurvreg()}.
#'   Case-insensitive; \code{"exponential"} is accepted for \code{"exp"}.
#'
#' @return A function \code{H(x, ...)} returning cumulative hazards.
#' @export
#'
#' @examplesIf requireNamespace("flexsurv", quietly = TRUE)
#' H <- flexsurv_cum_haz("weibull")
#' H(0:5, shape = 1.5, scale = 3)
flexsurv_cum_haz <- function(dist) {
  info <- flexsurv_dist_info(dist)
  fname <- paste0("H", info$name)
  ns <- asNamespace("flexsurv")
  assertthat::assert_that(
    exists(fname, envir = ns, mode = "function", inherits = FALSE),
    msg = sprintf(
      "flexsurv does not provide a cumulative hazard function '%s'",
      fname
    )
  )
  get(fname, envir = ns, mode = "function", inherits = FALSE)
}
