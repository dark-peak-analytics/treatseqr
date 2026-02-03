extrapolate_heseqmat <- function(m, spec) {
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
