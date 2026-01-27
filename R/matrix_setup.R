#' Calculate Matrix Size and Top Row Indices for Sequence Model
#'
#' Computes the indices of the top row (`m1`) of the main transition matrix and
#' the overall matrix size, given the number of cycles and tunnel states.
#'
#' @param n_cycles Integer. Number of cycles in the model (must be > 0).
#' @param n_tunnels Integer. Number of tunnel states in the model (must be > 0).
#'
#' @return A list with:
#'   \describe{
#'     \item{m1}{Integer vector of column indices for the top row of the
#'     matrix.}
#'     \item{matrix_size}{Integer. The total number of columns (states) in the
#'     matrix.}
#'   }
#' @examples
#' m_size(5, 2)
#' @export
m_size <- function(n_cycles, n_tunnels, pre_tunnel_states = 1) {
  assertthat::assert_that(
    is.numeric(n_cycles),
    length(n_cycles) == 1,
    n_cycles > 0,
    n_cycles == as.integer(n_cycles),
    msg = "n_cycles must be a single positive integer"
  )
  assertthat::assert_that(
    is.numeric(n_tunnels),
    length(n_tunnels) == 1,
    n_tunnels > 0,
    n_tunnels == as.integer(n_tunnels),
    msg = "n_tunnels must be a single positive integer"
  )
  assertthat::assert_that(
    length(pre_tunnel_states) == 1,
    is.numeric(pre_tunnel_states),
    pre_tunnel_states == floor(pre_tunnel_states),
    pre_tunnel_states >= 1,
    pre_tunnel_states <= 10,
    msg = "'pre_tunnel_states' must be a single positive integer. Limit of 10"
  )

  # pre-tunnel row/column indices. These are diagonal elements as probability of
  # staying in pre-tunnel states:
  pt_i <- seq_len(pre_tunnel_states)

  # starting indices for first tunnel's top left element. remember that this
  # element is empty as probability of stayin in tunnel is superdiagonal (1)
  tun_start <- c(pre_tunnel_states + 1, pre_tunnel_states + 1)

  # number of destinations from the first health state (pre-tunnel)
  n_dest_first <- n_tunnels + pre_tunnel_states + 1
  n_dest_last_non_tunnel <- n_dest_first - pre_tunnel_states + 1
  n_dest <- n_dest_first:n_dest_last_non_tunnel

  # indices for the top left of each tunnel block are the same for all
  # pre-tunnel states, so work them out just once:
  tunnel_starts <- seq(
    from = tun_start[2],
    to = tun_start[2] + (n_tunnels - 1) * n_cycles,
    by = n_cycles
  )

  # add one row/column at the extreme right and bottom for the dead state:
  dead_rowcol <- max(tunnel_starts) + n_cycles + 1

  tun_j <- c(tunnel_starts, dead_rowcol)

  # rows for m1. Note each successive pre-tunnel has one less "place to go".
  # This is relatively simple to derive:
  m1_i <- rep(
    seq_len(pre_tunnel_states),
    times = n_dest
  )

  # The columns are more involved. essentially columns to move forward to
  # another state, and then the set of tunnel "beginnings", plus dead for each
  # pre-tunnel.
  m1_ijx <- Reduce(
    x = pt_i,
    init = matrix(
      c(m1_i, rep(rep(0, length(m1_i)), 2)),
      ncol = 3,
      dimnames = list(NULL, c("i", "j", "x"))
    ),
    function(coord_mat, pre_tun_number) {
      # pull out the relevant rows of the coordinate matrix for this pre-tunnel
      relevant_rows <- which(coord_mat[, "i"] == pre_tun_number)

      # identify the movements patients can make for the pre-tunnel states
      pre_tun_moves <- pre_tun_number:pre_tunnel_states
      j_vals <- c(pre_tun_moves, tun_j)
      coord_mat[relevant_rows, "j"] <- j_vals
      coord_mat
    }
  )

  # return a list of a coordinate matrix for m1, and bounds for m2
  list(
    m1_ijx = m1_ijx,
    matrix_size = dead_rowcol
  )
}
