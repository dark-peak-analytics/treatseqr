#' Extrapolate Health Economic Sequence Model
#'
#' This function extrapolates a health economic sequence model using sparse
#' matrices. It generates a population trace matrix showing state occupancy
#' over time.
#'
#' @param m A list containing the transition matrices:
#'   - m1: A numeric array of dim c(th, n_pre, n_dest) holding the pre-tunnel
#'     transition probabilities, indexed
#'     `[cycle, pre-tunnel state, destination]`.
#'     Destinations are in canonical order (all states, then death).
#'   - m1_dest: Integer vector giving the column of the full matrix M that each
#'     destination slot of m1 corresponds to.
#'   - m2: A sparse matrix for tunnel states.
#'   - state_names: Character vector of state names, death last.
#' @param spec A model specification list containing:
#'   - n_cycles: Number of model cycles.
#'   - matrix_size: Size of the transition matrix.
#'
#' @return A matrix where rows represent states and columns represent time
#'   points (cycles). The matrix shows the proportion of the cohort in each
#'   state at each time point.
#'
#' @export
extrapolate_treatseqr <- function(m, spec) {
  # output population matrix:
  th <- dim(m$m1)[1]
  pop <- matrix(0, nrow = spec$matrix_size, ncol = th + 1)
  pop[1, 1] <- 1

  # get the TPs
  m1 <- m$m1
  m2t <- Matrix::t(m$m2)
  pre_pos <- seq_len(spec$pre_tunnels)
  pre_dest <- m$m1_dest

  # Extrapolate. Scatter-add for m1, matmult for m2
  for (cyc in seq_len(th) + 1) {
    pop[, cyc] <- as.numeric(m2t %*% pop[, cyc - 1, drop = FALSE])
    pop[pre_dest, cyc] <- pop[pre_dest, cyc] +
      as.numeric(pop[pre_pos, cyc - 1] %*% m1[cyc - 1, , ])
  }
  pop
}
