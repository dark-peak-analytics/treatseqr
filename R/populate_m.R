#' These functions are dedicated to using the validated version of the
#' transition probabilities which the functions in `validate_tp.R` produce to
#' generate the sparce matrices required to extrapolate the sequencing model
#'
#' See the vignette `sparse-matrix-functions.Rmd` for more details on how these
#' work.

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
  m2 <- xxx
}
