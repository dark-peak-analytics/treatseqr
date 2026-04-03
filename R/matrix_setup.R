#' Specify Matrix Indices and Size for Sequence Model with Tunnel States
#'
#' Computes the coordinate matrices for the main transition matrix of a sequence
#' model with tunnel states, as well as the overall matrix size, given the
#' number of cycles, tunnel states, and pre-tunnel states.
#'
#' @param n_cycles Integer. Number of cycles in each tunnel (must be > 0).
#' @param n_tunnels Integer. Number of tunnel state blocks (must be > 0).
#' @param pre_tunnel_states Integer. Number of pre-tunnel states before entering
#' the first tunnel (default is 1, must be between 1 and 10).
#' @param tunnel_lengths the lengths of the individual tunnels. If left blank
#' the function assumes that all the tunnels are `n_cycles` in length
#'
#' @return A list with:
#'   \describe{
#'     \item{m1_ijx}{A matrix of coordinates (i, j, x) for the pre-tunnel state
#'     transitions.}
#'     \item{m2_ijx}{A matrix of coordinates (i, j, x) for the tunnel state
#'     transitions and death state.}
#'     \item{matrix_size}{Integer. The total number of states (rows/columns) in
#'     the transition matrix.}
#'   }
#' @details
#' The function generates coordinate matrices for efficiently populating a
#' sparse transition matrix for a sequence model with tunnel states. The
#' \code{m1_ijx} matrix contains the coordinates for transitions from pre-tunnel
#' states, while \code{m2_ijx} contains the coordinates for transitions within
#' tunnel states and to the absorbing (death) state.
#'
#' @examples
#' specify_m(5, 2)
#' specify_m(3, 3, pre_tunnel_states = 2)
#' @export
specify_m <- function(
  n_cycles,
  n_tunnels,
  tunnel_lengths = NULL,
  pre_tunnel_states = 1
) {
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

  if (!is.null(tunnel_lengths)) {
    # if variable tunnel lengths, use the lengths provided
    tunnel_starts <- cumsum(tunnel_lengths) + tun_start[1]
    dead_rowcol <- sum(
      max(tunnel_starts),
      tunnel_lengths[length(tunnel_lengths)],
      1
    )
  } else {
    # indices for the top left of each tunnel block are the same for all
    # pre-tunnel states, so work them out just once:
    tunnel_starts <- seq(
      from = tun_start[2],
      to = tun_start[2] + (n_tunnels - 1) * n_cycles,
      by = n_cycles
    )
    # add one row/column at the extreme right and bottom for the dead state:
    dead_rowcol <- max(tunnel_starts) + n_cycles + 1
  }

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

  # now m2, which is for the successive tunnel states.

  # the top left column of each tunnel block, with the column index for the
  # vertical transitions to subsequent tunnel states, plus death
  tun_topleft <- lapply(seq_along(tun_j), function(i) {
    tun_j[i:length(tun_j)]
  })

  # we can cycle through the above to generate the i,j,x for m2:
  l_m2_ijx <- lapply(
    X = tun_topleft,
    FUN = function(tun_start_cols) {
      if (length(tun_start_cols) == 1) {
        # This is the death state. coordinates are bottom right element:
        list(i = dead_rowcol, j = dead_rowcol, x = 1)
      } else {
        # this is a tunnel block with n_cycles states. The first set are
        # superdiagonal (1), whilst the others are veritcally arranged

        # first element's row is also the leftmost column of this tunnel block
        tun_start_row <- tun_start_cols[1]

        # rows for superdiagonal elements go from that point. cols are +1
        tun_sdiag_i <- tun_start_row + (seq_len(n_cycles) - 1)
        tun_sdiag_j <- tun_sdiag_i + 1

        # rest of the transitions are vertically arranged:
        vertical_strips <- tun_start_cols[-1]

        # row indices are the same as for the superdiagonals:
        tun_vert_i <- rep(tun_sdiag_i, length(vertical_strips))

        # column indices just repeat element of vertical_strips n_cycles times
        tun_vert_j <- rep(vertical_strips, each = n_cycles)

        list(
          i = c(tun_sdiag_i, tun_vert_i),
          j = c(tun_sdiag_j, tun_vert_j),
          x = rep(0, length(tun_sdiag_i) + length(tun_vert_i))
        )
      }
    }
  )

  # collapse m2_ijx list into a single coordinate matrix. one can populate m2
  # with data simply using m2_ijx[, x] <- transition_probaibilites later on
  m2_ijx <- matrix(
    c(
      lapply(l_m2_ijx, function(x) .subset2(x, "i")) |> unlist(),
      lapply(l_m2_ijx, function(x) .subset2(x, "j")) |> unlist(),
      lapply(l_m2_ijx, function(x) .subset2(x, "x")) |> unlist()
    ),
    ncol = 3,
    dimnames = list(NULL, c("i", "j", "x"))
  )

  # return a list of a coordinate matrix for m1, and bounds for m2
  list(
    m1_ijx = m1_ijx,
    m2_ijx = m2_ijx,
    matrix_size = dead_rowcol,
    n_cycles = n_cycles,
    pre_tunnels = pre_tunnel_states
  )
}
