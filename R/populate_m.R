# These functions are dedicated to using the validated version of the
# transition probabilities which the functions in `validate_tp.R` produce to
# generate the sparce matrices required to extrapolate the sequencing model
#
# See the vignette `sparse-matrix-functions.Rmd` for more details on how these
# work.

#' Generate Markov Transition Matrices (m1 and m2) Based on Model Specification
#'
#' This function generates two transition matrices (m1 and m2) for a Markov
#' model, using the provided model specification and a validated list of
#' transition probabilities.
#' - m1 corresponds to pre-tunnel states for each model cycle.
#' - m2 corresponds to tunnel states after entering the sequence.
#'
#' @param m_specification A list containing model parameters, as returned by
#'   \code{specify_m()}. Required fields:
#'   - tunnel_lengths: Integer vector of per-tunnel cycle counts.
#'   - pre_tunnels: Number of pre-tunnel states.
#'   - matrix_size: Size of the transition matrix.
#'   - m1_ijx: Coordinate matrix for m1.
#'   - m2_ijx: Coordinate matrix for m2.
#' @param transition_prob_list A list of transition probabilities for each
#' state.
#' @param state_names A character vector used to validate
#' \code{transition_prob_list}. It must equal \code{names(transition_prob_list)}
#' in the same order, with the absorbing (death) state appended as the last
#' element. It cannot be used to reorder states: the matrix layout produced by
#' \code{specify_m()} is positional, so a vector in any other order is an
#' error rather than a request to permute the columns of M.
#'
#' @return A list with four elements:
#'   - m1: A numeric array of dim c(th, n_pre, n_dest) holding the pre-tunnel
#'     transition probabilities, indexed
#'     `[cycle, pre-tunnel state, destination]`.
#'     Destinations are in canonical order (all states, then death).
#'   - m1_dest: Integer vector giving the column of the full matrix M that each
#'     destination slot of m1 corresponds to.
#'   - m2: A sparse matrix for tunnel states.
#'   - state_names: Character vector of state names, death last.
#'
#' @details
#' The function validates the transition probability list, splits it into
#' pre-tunnel and tunnel components, and constructs the corresponding sparse
#' matrices. It ensures that row sums are appropriate and asserts that the
#' tunnel matrix rows sum to 1. Note that it filters out zero probabilities for
#' m2 before entering them into the sparse matrix. This avoids those values
#' being put into the matrix and then used in the matrix multiplications, which
#' would ultimately always result in values of 0 being propagated. This then has
#' a computational gain for no cost.
#'
#' The time horizon is taken as the longest transition probability vector
#' supplied for the first state, and must be at least 2 cycles. A single-cycle
#' horizon is rejected, as it leaves no transition for the engine to apply.
#'
#' @importFrom Matrix sparseMatrix rowSums
#' @importFrom assertthat assert_that
#' @export
#'
generate_m_list <- function(
  m_specification,
  transition_prob_list,
  state_names
) {
  # validate the transition probability source list.
  # This will return a valid list or an error.
  valid_tp <- validate_tp_source(
    tp_source = transition_prob_list,
    tunnel_lengths = m_specification$tunnel_lengths,
    state_names = state_names
  )

  # go through valid_tp making sure all of state_names except self are inside
  # each element:
  for (o_state in state_names[-length(state_names)]) {
    missing_states <- unique(
      state_names[
        state_names != o_state & !state_names %in% names(valid_tp[[o_state]])
      ]
    )
    assertthat::assert_that(
      length(missing_states) == 0,
      msg = sprintf(
        paste0(
          "State(s) %s are not present in the validated",
          " transition probability list for '%s'"
        ),
        paste(shQuote(missing_states), collapse = ", "),
        o_state
      )
    )

    # check ordering (excluding the origin state itself)
    expected_order <- state_names[state_names != o_state]
    actual_order <- names(valid_tp[[o_state]])
    # ensure origin is removed from the actual order if present
    actual_order <- actual_order[actual_order != o_state]
    if (!identical(expected_order, actual_order)) {
      msg <- sprintf(
        "Ordering mismatch for '%s'. Expected (excluding origin): %s; got: %s",
        o_state,
        paste(shQuote(expected_order), collapse = ", "),
        paste(shQuote(actual_order), collapse = ", ")
      )
      assertthat::assert_that(FALSE, msg = msg)
    }
  }

  # time horizon must be the length of the pre-tunnels. At least 1 pre-tunnel:
  th <- max(unlist(lapply(transition_prob_list[[1]], length)))

  # time horizon must be at least 2
  assertthat::assert_that(
    th >= 2,
    msg = sprintf("Time horizon (%d) must be at least 2 cycles", th)
  )

  # Other important inputs
  pre_tun_states <- m_specification$pre_tunnels
  n_states <- length(valid_tp)
  matrix_size <- m_specification$matrix_size
  pre_tun_names <- names(valid_tp)[seq_len(pre_tun_states)]

  # coordinates for m1 and m2:
  m1_ijx <- m_specification$m1_ijx
  m2_ijx <- m_specification$m2_ijx

  # for pre-tunnels, we just need all the elements from the first pre_tun_states
  # bits of transition_prob_list. These go into m1
  tp_pre_tun <- valid_tp[1:pre_tun_states]

  # all states are tunnels after entering the sequence. These all go into m2
  tp_tun <- valid_tp[(pre_tun_states + 1):(n_states)]

  # generate a 3d array with 2nd dimension pre-tunnel state, 3rd dimesnion
  # destination state (1st is model cycle)
  dest_names <- c(names(valid_tp), tail(state_names, 1))
  n_dest <- length(dest_names)
  m1_dest <- sort(unique(m1_ijx[, "j"]))

  # make the 3d array for all pre-tunnel states:
  m1_array <- array(
    data = 0,
    dim = c(th, pre_tun_states, n_dest)
  )

  # populate it by pulling the numbers in:
  for (pre_tun in seq_len(pre_tun_states)) {
    tp_move <- .subset2(tp_pre_tun, pre_tun)
    cycle_mat <- matrix(0, nrow = th, ncol = n_dest)
    cycle_mat[, match(names(tp_move), dest_names)] <- do.call(cbind, tp_move)

    # p_stay is complement (1 - sum):
    cycle_mat[, pre_tun] <- 1 - rowSums(cycle_mat)
    m1_array[, pre_tun, ] <- cycle_mat
  }

  # Now that m1 is populated, let's move on to m2. It is quite simple because we
  # have already validated the transition probability list
  x_m2 <- unlist(
    use.names = FALSE,
    lapply(seq_along(tp_tun), function(tunnel) {
      self_name <- names(tp_tun)[tunnel]
      state_dests <- tp_tun[[tunnel]]

      # validation
      dropped <- state_dests[names(state_dests) %in% pre_tun_names]
      offending <- names(dropped)[vapply(
        dropped,
        function(v) any(v != 0),
        logical(1L)
      )]
      assertthat::assert_that(
        length(offending) == 0L,
        msg = sprintf(
          paste0(
            "State '%s' specifies non-zero transitions to pre-tunnel state(s) %s, ",
            "which are not supported: tunnel states cannot transition back to ",
            "pre-tunnel states."
          ),
          self_name,
          paste(sQuote(offending), collapse = ", ")
        )
      )

      # Exclude self (p_stay computed as complement) and pre-tunnel destinations
      # (no m2 coordinate slots exist for pre-tunnel columns)
      keep_names <- names(state_dests)[
        !names(state_dests) %in% c(pre_tun_names, self_name)
      ]
      ordered_dests <- state_dests[keep_names]
      p_stay <- 1 - Reduce(`+`, ordered_dests)
      unlist(c(list(p_stay = p_stay), ordered_dests), use.names = FALSE)
    })
  )

  m2_ijx[seq_len(nrow(m2_ijx) - 1), "x"] <- x_m2

  m2_ijx_nz <- m2_ijx[m2_ijx[, "x"] != 0, , drop = FALSE]
  sm_m2 <- Matrix::sparseMatrix(
    i = m2_ijx_nz[, "i"],
    j = m2_ijx_nz[, "j"],
    x = m2_ijx_nz[, "x"],
    dims = c(matrix_size, matrix_size)
  )

  # quick assert check:
  assertthat::assert_that(
    all(
      abs(Matrix::rowSums(sm_m2[(pre_tun_states + 1):nrow(sm_m2), ]) - 1) <
        1e-10
    ),
    msg = "Row sums of m2 do not equal 1 after populating transition probabilities"
  )

  list(
    m1 = m1_array,
    m1_dest = m1_dest,
    m2 = sm_m2,
    state_names = state_names
  )
}
