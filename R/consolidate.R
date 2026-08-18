# Functions to consolidate `trace` after it is calculated using
# `extrapolate_treatseqr`

#' Consolidate a full trace matrix into summarized tunnel and pre-tunnel blocks
#'
#' This function takes a full trace matrix from a health economic state
#' transition model and consolidates it into summarized blocks for pre-tunnel,
#' tunnel, and dead states. It returns two matrices: one with time-based (t) and
#' one with difference-based (d) summaries.
#'
#' @param full_trace A matrix representing the full trace of state occupancies
#' over time.
#' @param spec A list containing model specifications, including:
#'   - pre_tunnels: Number of pre-tunnel states
#'   - tunnel_lengths: Integer vector of tunnel block lengths
#' @param state_names A character vector of state names for labeling output
#' matrices.
#' @param m Optional list as returned by \code{generate_m_list()}. When
#' supplied, the 3d array of pre-tunnel transition probabilities are used to
#' compute true time-in-state (sojourn) curves for the pre-tunnel columns of
#' \code{d}, weighting each entry cohort by its size. When \code{NULL}, the
#' pre-tunnel columns of \code{d} fall back to wall-time occupancy, which only
#' coincides with a sojourn curve for a single pre-tunnel start state; a
#' warning is issued if \code{pre_tunnels > 1}.
#'
#' @return A list with two elements:
#'   - t: Matrix of time-based state occupancies (rows: cycles, columns: states)
#'   - d: Matrix of difference-based state occupancies (rows: cycles, columns:
#'   states)
#'
#' @details
#' - Pre-tunnel blocks are summarized as single rows.
#' - Tunnel blocks are summarized by summing across the appropriate rows for
#' each tunnel.
#' - The dead block is taken as the last row of the trace.
#' - The function is useful for cost-effectiveness modeling where tunnel states
#' are present.
#'
#' The \code{d} matrix is indexed by cycles since state entry: tunnel columns
#' are person-cycles at each tunnel position normalized by total traffic into
#' the tunnel, pre-tunnel columns are entry-cohort-weighted survival-in-state
#' curves (requires \code{m}), and the dead column is overall survival. Every
#' column of \code{d} therefore starts at 1 (or 0 for states never entered).
#' As with tunnels, entry cohorts are aggregated across wall time, so curves
#' mix cohorts entering at different model cycles.
#'
#' A tunnel column of \code{d} is only a true time-in-state curve while
#' patients cannot accumulate at a single \code{d}. If the transition
#' probabilities leave a non-zero \code{p_stay} at a tunnel's final cycle (a
#' "holding state"), patients pile up in that last position and \code{d} stops
#' measuring time since entry for them. To keep the interpretation clean, close
#' a fixed-duration tunnel with a decision-tree node instead: set the final
#' cycle's destination probabilities to sum to 1 so \code{p_stay} is 0. See the
#' vignette \code{fixed-duration-tunnels.Rmd} for a worked comparison of the two
#' designs.
#'
#' @export
consolidate_treatseqr_trace <- function(
  full_trace,
  spec,
  state_names,
  m = NULL
) {
  pre_tun <- spec$pre_tunnels
  n_states <- pre_tun + length(spec$tunnel_lengths) + 1L
  tunnel_starts <- pre_tun +
    cumsum(c(1L, utils::head(spec$tunnel_lengths, -1L)))
  tun_blocks <- lapply(seq_along(tunnel_starts), function(k) {
    c(tunnel_starts[k], tunnel_starts[k] + spec$tunnel_lengths[k] - 1L)
  })

  # th is trace columns less 1 for dead
  th <- ncol(full_trace) - 1

  # pre-tunnel blocks are just one row each, so just a sequence is fine:
  pre_blocks <- seq_len(pre_tun)

  # compute trace values and give names:
  trace_pre <- t(full_trace[pre_blocks, , drop = FALSE])

  # Dead block is just the last row:
  trace_dead <- t(full_trace[nrow(full_trace), , drop = FALSE])

  # tunnel populations:
  trace_tun_t <- lapply(tun_blocks, function(tun_idx) {
    index <- seq(tun_idx[1], tun_idx[2])
    colSums(full_trace[index, , drop = FALSE])
  })
  trace_tun_d <- lapply(tun_blocks, function(tun_idx) {
    index <- seq(tun_idx[1], tun_idx[2])
    total <- sum(full_trace[index[1], ])
    if (total == 0) {
      d_vec <- rep(0, length(index))
    } else {
      d_vec <- rowSums(full_trace[index, , drop = FALSE]) / total
    }
    c(d_vec, rep(0, th - length(d_vec)))
  })

  # pre-tunnel d columns. Tunnels encode time-in-state positionally in the
  # trace, but a pre-tunnel state is a single row, so its sojourn curve needs
  # the per-cycle staying probabilities from m$m1.
  if (is.null(m)) {
    if (pre_tun > 1) {
      warning(
        "With multiple pre-tunnel states, pre-tunnel columns of `d` are ",
        "wall-time occupancy, not time-in-state curves. Supply `m` to ",
        "compute true sojourn curves."
      )
    }
    trace_pre_d <- trace_pre[-(th + 1), , drop = FALSE]
  } else {
    assertthat::assert_that(
      dim(m$m1)[1] == th,
      msg = "dim(m$m1)[1] must match the trace time horizon (ncol - 1)"
    )
    # p_stay[s, u]: probability that pre-tunnel state s retains its occupants
    # across the transition from time u - 1 to time u
    p_stay <- t(vapply(pre_blocks, function(s) m$m1[, s, s], numeric(th)))

    trace_pre_d <- vapply(
      pre_blocks,
      FUN.VALUE = numeric(th),
      FUN = function(s) {
        occ <- full_trace[s, ]
        ps <- p_stay[s, ]

        # entrants at each time: initial occupancy, then any occupancy not
        # explained by last cycle's occupants staying put (guard tiny negative
        # values from floating point)
        inflow <- c(occ[1], occ[-1] - occ[-(th + 1)] * ps)
        inflow[inflow < 0] <- 0
        total <- sum(inflow)
        if (total == 0) {
          return(rep(0, th))
        }

        # aggregate each entry cohort's survival-in-state, mirroring the
        # tunnel normalization: numerator only includes cohorts observable
        # within the horizon, denominator is total traffic into the state
        d_vec <- numeric(th)
        d_vec[1] <- 1
        surv <- rep(1, th + 1)
        for (k in seq_len(th)[-1]) {
          cohort <- seq_len(th - k + 2)
          surv[cohort] <- surv[cohort] * ps[cohort + k - 2]
          d_vec[k] <- sum(inflow[cohort] * surv[cohort]) / total
        }
        d_vec
      }
    )
    trace_pre_d <- matrix(trace_pre_d, nrow = th)
  }

  # return list with 2 matrices, one for t-based, one for d based:
  t_mat <- matrix(
    c(trace_pre, unlist(trace_tun_t, use.names = FALSE), trace_dead),
    ncol = n_states,
    nrow = th + 1,
    dimnames = list(
      NULL,
      state_names
    )
  )

  # Note that this is flexible to fixed-duration tunnels

  d_mat <- matrix(
    c(
      trace_pre_d,
      unlist(trace_tun_d, use.names = FALSE),
      1 - trace_dead[-(th + 1)]
    ),
    ncol = n_states,
    nrow = th,
    dimnames = list(
      NULL,
      state_names
    )
  )

  # return these traces, as they are both very useful for CE modelling.
  list(
    t = t_mat,
    d = d_mat
  )
}
