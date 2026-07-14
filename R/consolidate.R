#' Functions to consolidate `trace` after it is calculated using
#' `extrapolate_treatseqr`

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
#'   - m1_ijx: Matrix with state indices and bounds
#' @param state_names A character vector of state names for labeling output
#' matrices.
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
#' @export
consolidate_treatseqr_trace <- function(full_trace, spec, state_names) {
  pre_tun <- spec$pre_tunnels
  bounds <- spec$m1_ijx[, "j"]

  # th is trace columns less 1 for dead
  th <- ncol(full_trace) - 1
  n_states <- nrow(spec$m1_ijx)

  # pre-tunnel blocks are just one row each, so just a sequence is fine:
  pre_blocks <- seq_len(pre_tun)

  # within-tunnel blocks are pairings of first row and last row, which go from
  # element pre_tun + 1 to length(bounds). However, each 2nd item needs reducing
  # by 1 because that is the starting point of the next tunnel not the end of
  # this one.
  tun_indices <- (pre_tun + 1):(length(bounds) - 1)

  # Generate the coordinate bounds for each tunnel block. This is the same when
  # summing across rows, or in columns restricting to rows.
  tun_blocks <- lapply(
    seq_along(tun_indices),
    function(k) {
      tunnel_start <- bounds[tun_indices[k]]
      c(tunnel_start, tunnel_start + spec$tunnel_lengths[k] - 1)
    }
  )

  # compute trace values and give names:
  trace_pre <- t(full_trace[pre_blocks, ])

  # Dead block is just the last row:
  trace_dead <- t(full_trace[nrow(full_trace), ])

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
      trace_pre[-(th + 1)],
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
