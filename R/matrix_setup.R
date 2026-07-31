#' Specify Matrix Indices and Size for Sequence Model with Tunnel States
#'
#' Computes the coordinate matrices for the main transition matrix of a sequence
#' model with tunnel states, as well as the overall matrix size, given a vector
#' of tunnel lengths and the number of pre-tunnel states.
#'
#' @param tunnel_lengths Integer vector. Length of each tunnel state block (one
#'   element per tunnel; all values must be positive integers). The number of
#'   tunnels is inferred from the length of this vector.
#' @param pre_tunnel_states Integer. Number of pre-tunnel states before entering
#' the first tunnel (default is 1, must be between 1 and 10).
#'
#' @return A list with:
#'   \describe{
#'     \item{m1_ijx}{A matrix of coordinates (i, j, x) for the pre-tunnel state
#'     transitions.}
#'     \item{m2_ijx}{A matrix of coordinates (i, j, x) for the tunnel state
#'     transitions and death state.}
#'     \item{matrix_size}{Integer. The total number of states (rows/columns) in
#'     the transition matrix.}
#'     \item{tunnel_lengths}{The \code{tunnel_lengths} vector as supplied.}
#'     \item{pre_tunnels}{The \code{pre_tunnel_states} value as supplied.}
#'   }
#' @details
#' The function generates coordinate matrices for efficiently populating a
#' sparse transition matrix for a sequence model with tunnel states. The
#' \code{m1_ijx} matrix contains the coordinates for transitions from pre-tunnel
#' states, while \code{m2_ijx} contains the coordinates for transitions within
#' tunnel states and to the absorbing (death) state.
#'
#' Tunnels may have different lengths (variable-duration tunnels). The total
#' matrix size is \code{pre_tunnel_states + sum(tunnel_lengths) + 1}.
#'
#' @examples
#' # Two tunnels of equal length (5 cycles each)
#' specify_m(rep(5, 2))
#'
#' # Three tunnels of different lengths
#' specify_m(c(5, 10, 15), pre_tunnel_states = 2)
#'
#' # Multiple pre-tunnel states
#' specify_m(rep(3, 3), pre_tunnel_states = 2)
#' @importFrom utils head
#' @export
specify_m <- function(tunnel_lengths, pre_tunnel_states = 1) {
  assertthat::assert_that(
    is.numeric(tunnel_lengths),
    length(tunnel_lengths) >= 1,
    all(tunnel_lengths == as.integer(tunnel_lengths)),
    all(tunnel_lengths > 0),
    msg = "tunnel_lengths must be a vector of positive integers"
  )
  assertthat::assert_that(
    length(pre_tunnel_states) == 1,
    is.numeric(pre_tunnel_states),
    pre_tunnel_states == floor(pre_tunnel_states),
    pre_tunnel_states >= 1,
    pre_tunnel_states <= 10,
    msg = "'pre_tunnel_states' must be a single positive integer. Limit of 10"
  )

  # tunnel lengths shows us how many tunnels:
  n_tunnels <- length(tunnel_lengths)

  # pre-tunnel row/column indices. These are diagonal elements as probability of
  # staying in pre-tunnel states:
  pt_i <- seq_len(pre_tunnel_states)

  # number of destinations from each pre-tunnel state (uniform: all pre-tunnel
  # states have the same destination set — self, all other pre-tunnel states,
  # all tunnel starts, and death):
  n_dest_first <- n_tunnels + pre_tunnel_states + 1

  # the tunnels are different sizes, so compute the start points for each.
  # this can be done by taking n(pre-tunnel), then adding up the sizing of all
  # tunnels except the last one:
  tunnel_starts <- pre_tunnel_states + cumsum(c(1L, head(tunnel_lengths, -1L)))

  # add one row/column at the extreme right and bottom for the dead state:
  dead_rowcol <- max(tunnel_starts) + tunnel_lengths[length(tunnel_lengths)]

  tun_j <- c(tunnel_starts, dead_rowcol)

  # rows for m1. Every pre-tunnel state has the same number of destinations
  # (uniform). Backwards pre-tunnel coordinates are included here; they default
  # to 0 and are filtered out at the matrix population stage if not used.
  m1_i <- rep(
    seq_len(pre_tunnel_states),
    each = n_dest_first
  )

  # The columns: self first (maps to p_stay in populate_m), then all other
  # pre-tunnel states in index order (backwards then forwards), then all tunnel
  # starts, then death.
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

      # self first, then all other pre-tunnel states in index order, then tunnels
      pre_tun_moves <- c(
        pre_tun_number,
        setdiff(seq_len(pre_tunnel_states), pre_tun_number)
      )
      j_vals <- c(pre_tun_moves, tun_j)
      coord_mat[relevant_rows, "j"] <- j_vals
      coord_mat
    }
  )

  # now m2, which is for the successive tunnel states.

  # the column indices for the vertical strips in each tunnel block.
  # Each tunnel block gets strips for ALL other tunnels + death (not just
  # forward). Backwards strips default to 0 and are filtered at population
  # stage if unused. The death block entry is just itself (absorbing state).
  tun_topleft <- lapply(seq_along(tun_j), function(i) {
    if (i == length(tun_j)) {
      # death: just itself
      tun_j[i]
    } else {
      tun_j[-i]
    } # all other tunnel starts + death
  })

  # we can cycle through the above to generate the i,j,x for m2. Note that in
  # each tunnel, the last transition is placed on the diagonal to give users the
  # flexibility to "hold" patients at the end of a tunnel, or set the tp_source
  # to have full probability of exit at the end.
  l_m2_ijx <- lapply(
    X = seq_along(tun_topleft),
    FUN = function(i_tunnel_block) {
      if (i_tunnel_block == length(tun_topleft)) {
        # This is the death state. Coordinates are bottom right element:
        list(i = dead_rowcol, j = dead_rowcol, x = 1)
      } else {
        # this is a tunnel block with user-defined size. The first set are
        # superdiagonal (1), whilst the others are vertically arranged. The last
        # transition in the tunnel goes on the diagonal to allow users to
        # specify "holding" at the end of the tunnel if they wish

        # pull the block for readability:
        block <- tun_topleft[[i_tunnel_block]]

        # first element's row is also the leftmost column of this tunnel block
        tun_start_row <- tunnel_starts[i_tunnel_block]
        tun_len <- tunnel_lengths[i_tunnel_block]

        # rows for superdiagonal elements go from that point. cols are +1
        tun_sdiag_i <- tun_start_row + (seq_len(tun_len) - 1)
        tun_sdiag_j <- tun_sdiag_i + 1

        # each block self-loops at the end of its extent, instead of imposing
        # any transitions. This means that the user can specify probability of
        # exit as 1 themselves instead of it being forced by the function.
        tun_sdiag_j[tun_len] <- tun_sdiag_i[tun_len]

        # vertical strips for all other tunnels + death. tun_topleft already
        # excludes self, so no [-1] needed here.
        vertical_strips <- block

        # row indices are the same as for the superdiagonals:
        tun_vert_i <- rep(tun_sdiag_i, length(vertical_strips))

        # column indices: j = repeat horizontal location tun_len times
        tun_vert_j <- rep(vertical_strips, each = tun_len)

        # return a list for consolidation later.
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
    tunnel_lengths = tunnel_lengths,
    pre_tunnels = pre_tunnel_states
  )
}
