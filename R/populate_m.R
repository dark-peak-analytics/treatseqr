#' These functions are dedicated to using the validated version of the
#' transition probabilities which the functions in `validate_tp.R` produce to
#' generate the sparce matrices required to extrapolate the sequencing model
#'
#' See the vignette `sparse-matrix-functions.Rmd` for more details on how these
#' work.

#' Generate Markov Transition Matrices (m1 and m2) Based on Model Specification
#'
#' This function generates two transition matrices (m1 and m2) for a Markov
#' model, using the provided model specification and a validated list of
#' transition probabilities.
#' - m1 corresponds to pre-tunnel states for each model cycle.
#' - m2 corresponds to tunnel states after entering the sequence.
#'
#' @param m_specification A list containing model parameters, including:
#'   - n_cycles: Number of model cycles.
#'   - pre_tunnels: Number of pre-tunnel states.
#'   - matrix_size: Size of the transition matrix.
#'   - m1_ijx: Coordinate matrix for m1.
#'   - m2_ijx: Coordinate matrix for m2.
#' @param transition_prob_list A list of transition probabilities for each
#' state.
#'
#' @return A list with two elements:
#'   - m1: A list of sparse matrices (one per cycle) for pre-tunnel states.
#'   - m2: A sparse matrix for tunnel states.
#'
#' @details
#' The function validates the transition probability list, splits it into
#' pre-tunnel and tunnel components, and constructs the corresponding sparse
#' matrices. It ensures that row sums are appropriate and asserts that the
#' tunnel matrix rows sum to 1.
#'
#' @importFrom Matrix sparseMatrix rowSums
#' @importFrom assertthat assert_that
#'
generate_m_list <- function(m_specification, transition_prob_list) {
  # validate the transition probability source list.
  # This will return a valid list or an error.
  valid_tp <- validate_tp_source(transition_prob_list)

  # some useful variables to use throughout the funciton
  th <- m_specification$n_cycles
  pre_tun_states <- m_specification$pre_tunnels
  n_states <- length(valid_tp)
  tunnels <- n_states - pre_tun_states - 1
  matrix_size <- m_specification$matrix_size

  # coordinates for m1 and m2:
  m1_ijx <- m_specification$m1_ijx
  m2_ijx <- m_specification$m2_ijx

  # for pre-tunnels, we just need all the elements from the first pre_tun_states
  # bits of transition_prob_list. These go into m1
  tp_pre_tun <- valid_tp[1:pre_tun_states]

  # all states are tunnels after enteirng the sequence. these all go into m2
  tp_tun <- valid_tp[(pre_tun_states + 1):(n_states)]

  # Now that we have done this, we know that the list of transition
  # probabilities that have been supplied will fit into m1 and m2 correctly, and
  # that the rowSums of the resulting matrix are appropriate prior to
  # calculating the super-diagonal (1) elements for tunnel states

  m1_list <- lapply(seq_len(th), function(model_cycle) {
    coord <- m1_ijx
    x_values <- lapply(seq_len(pre_tun_states), function(state_index) {
      tp_move <- as.numeric(vapply(
        .subset2(tp_pre_tun, state_index),
        FUN.VALUE = numeric(1),
        FUN = function(x) x[model_cycle]
      ))
      tp_stay <- 1 - sum(tp_move, na.rm = TRUE)
      c(tp_stay, tp_move)
    })

    coord[, "x"] <- unlist(x_values)
    coord
    Matrix::sparseMatrix(
      i = coord[, "i"],
      j = coord[, "j"],
      x = coord[, "x"],
      dims = c(pre_tun_states, matrix_size)
    )
  })

  # Now that m1 is populated, let's move on to m2. It is quite simple because we
  # have already validated the transition probability list
  x_m2 <- unlist(
    use.names = FALSE,
    lapply(tp_tun, function(tunnel_state) {
      # in each tunnel state, add the prob of staying first and then return
      # collapsed list:
      p_stay <- 1 - Reduce(`+`, tunnel_state)
      tunnel_state <- c(list(p_stay = p_stay), tunnel_state)
      unlist(tunnel_state, use.names = FALSE)
    })
  )

  m2_ijx[seq_len(nrow(m2_ijx) - 1), "x"] <- x_m2

  sm_m2 <- Matrix::sparseMatrix(
    i = m2_ijx[, "i"],
    j = m2_ijx[, "j"],
    x = m2_ijx[, "x"],
    dims = c(matrix_size, matrix_size)
  )
  sm_m2[nrow(sm_m2) - 1, ncol(sm_m2)] <- 1

  # quick assert check:
  assertthat::assert_that(
    all(abs(Matrix::rowSums(sm_m2[2:nrow(sm_m2), ]) - 1) < 1e-10),
    msg = "Row sums of m2 do not equal 1 after populating transition probabilities"
  )

  list(
    m1 = m1_list,
    m2 = sm_m2
  )
}
