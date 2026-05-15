#' Extrapolate Health Economic Sequence Model
#'
#' This function extrapolates a health economic sequence model using sparse
#' matrices. It generates a population trace matrix showing state occupancy
#' over time.
#'
#' @param m A list containing the transition matrices:
#'   - m1: A list of sparse matrices (one per cycle) for pre-tunnel states.
#'   - m2: A sparse matrix for tunnel states.
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
  th <- spec$n_cycles
  pop <- matrix(0, nrow = spec$matrix_size, ncol = th + 1)
  pop[1, 1] <- 1

  # get the TPs
  m1 <- m$m1
  m2 <- m$m2

  # extrapolate the model:
  for (cyc in seq_len(th) + 1) {
    pop[, cyc] <-
      as.numeric(pop[1, cyc - 1] %*% m1[[cyc - 1]]) +
      as.numeric(pop[, cyc - 1] %*% m2)
  }
  pop
}
