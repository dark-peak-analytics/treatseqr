test_that("consolidate_treatseqr_trace returns correct structure", {
  m_spec <- specify_m(tunnel_lengths = rep(5, 2), pre_tunnel_states = 1)

  tp_list <- list(
    state1 = list(
      state2 = rep(0.2, 5),
      state3 = rep(0.1, 5),
      die = rep(0.05, 5)
    ),
    state2 = list(
      state3 = rep(0.3, 5),
      die = rep(0.1, 5)
    ),
    state3 = list(
      die = rep(0.2, 5)
    )
  )

  m <- generate_m_list(m_spec, tp_list)
  full_trace <- extrapolate_treatseqr(m, m_spec)

  state_names <- c("state1", "state2", "state3", "dead")

  result <- consolidate_treatseqr_trace(full_trace, m_spec, state_names)

  expect_type(result, "list")
  expect_named(result, c("t", "d"))
  expect_true(is.matrix(result$t))
  expect_true(is.matrix(result$d))
})

test_that("consolidate_treatseqr_trace has correct dimensions", {
  m_spec <- specify_m(tunnel_lengths = rep(10, 3), pre_tunnel_states = 1)

  tp_list <- list(
    s1 = list(
      s2 = rep(0.1, 10),
      s3 = rep(0.1, 10),
      s4 = rep(0.1, 10),
      die = rep(0.05, 10)
    ),
    s2 = list(
      s3 = rep(0.2, 10),
      s4 = rep(0.1, 10),
      die = rep(0.1, 10)
    ),
    s3 = list(
      s4 = rep(0.15, 10),
      die = rep(0.15, 10)
    ),
    s4 = list(
      die = rep(0.25, 10)
    )
  )

  m <- generate_m_list(m_spec, tp_list)
  full_trace <- extrapolate_treatseqr(m, m_spec)
  state_names <- c("s1", "s2", "s3", "s4", "dead")

  result <- consolidate_treatseqr_trace(full_trace, m_spec, state_names)

  # th = ncol(full_trace) - 1 = 10
  th <- ncol(full_trace) - 1L
  # t matrix should have th + 1 rows (includes time 0)
  expect_equal(nrow(result$t), th + 1)
  # d matrix should have th rows (no time 0)
  expect_equal(nrow(result$d), th)
  # Both should have correct number of columns (states)
  expect_equal(ncol(result$t), length(state_names))
  expect_equal(ncol(result$d), length(state_names))
})

test_that("consolidate_treatseqr_trace column names are correct", {
  m_spec <- specify_m(tunnel_lengths = c(3), pre_tunnel_states = 1)

  tp_list <- list(
    healthy = list(
      sick = rep(0.2, 3),
      die = rep(0.1, 3)
    ),
    sick = list(
      die = rep(0.3, 3)
    )
  )

  m <- generate_m_list(m_spec, tp_list)
  full_trace <- extrapolate_treatseqr(m, m_spec)
  state_names <- c("healthy", "sick", "dead")

  result <- consolidate_treatseqr_trace(full_trace, m_spec, state_names)

  expect_equal(colnames(result$t), state_names)
  expect_equal(colnames(result$d), state_names)
})

test_that("consolidate_treatseqr_trace t matrix rows sum to 1", {
  m_spec <- specify_m(tunnel_lengths = rep(8, 2), pre_tunnel_states = 1)

  tp_list <- list(
    on_tx = list(
      off_tx = rep(0.15, 8),
      prog = rep(0.1, 8),
      die = rep(0.05, 8)
    ),
    off_tx = list(
      prog = rep(0.2, 8),
      die = rep(0.1, 8)
    ),
    prog = list(
      die = rep(0.25, 8)
    )
  )

  m <- generate_m_list(m_spec, tp_list)
  full_trace <- extrapolate_treatseqr(m, m_spec)
  state_names <- c("on_tx", "off_tx", "prog", "dead")

  result <- consolidate_treatseqr_trace(full_trace, m_spec, state_names)

  row_sums <- rowSums(result$t)
  expect_true(all(abs(row_sums - 1) < 1e-10))
})

test_that("consolidate_treatseqr_trace d matrix excludes time 0", {
  m_spec <- specify_m(tunnel_lengths = rep(5, 2), pre_tunnel_states = 1)

  tp_list <- list(
    s1 = list(
      s2 = rep(0.2, 5),
      s3 = rep(0.15, 5),
      die = rep(0.05, 5)
    ),
    s2 = list(
      s3 = rep(0.25, 5),
      die = rep(0.1, 5)
    ),
    s3 = list(
      die = rep(0.2, 5)
    )
  )

  m <- generate_m_list(m_spec, tp_list)
  full_trace <- extrapolate_treatseqr(m, m_spec)
  state_names <- c("s1", "s2", "s3", "dead")

  result <- consolidate_treatseqr_trace(full_trace, m_spec, state_names)

  # d matrix should have one less row than t matrix (excludes time 0)
  expect_equal(nrow(result$d), nrow(result$t) - 1)

  # d matrix is based on within-tunnel time, so check valid probabilities
  expect_true(all(result$d >= 0))
  expect_true(all(result$d <= 1))
})

test_that("consolidate_treatseqr_trace handles different numbers of tunnels", {
  m_spec <- specify_m(tunnel_lengths = rep(4, 3), pre_tunnel_states = 1)

  tp_list <- list(
    s1 = list(
      s2 = rep(0.1, 4),
      s3 = rep(0.1, 4),
      s4 = rep(0.05, 4),
      die = rep(0.05, 4)
    ),
    s2 = list(
      s3 = rep(0.15, 4),
      s4 = rep(0.1, 4),
      die = rep(0.1, 4)
    ),
    s3 = list(
      s4 = rep(0.2, 4),
      die = rep(0.15, 4)
    ),
    s4 = list(
      die = rep(0.2, 4)
    )
  )

  m <- generate_m_list(m_spec, tp_list)
  full_trace <- extrapolate_treatseqr(m, m_spec)
  state_names <- c("s1", "s2", "s3", "s4", "dead")

  result <- consolidate_treatseqr_trace(full_trace, m_spec, state_names)

  # th = ncol(full_trace) - 1 = 4
  th <- ncol(full_trace) - 1L
  expect_equal(ncol(result$t), 5)
  expect_equal(ncol(result$d), 5)
  expect_equal(nrow(result$t), th + 1)
  expect_equal(nrow(result$d), th)
})

test_that("consolidate_treatseqr_trace dead state accumulates over time", {
  m_spec <- specify_m(tunnel_lengths = c(10), pre_tunnel_states = 1)

  tp_list <- list(
    alive = list(
      tun1 = rep(0.05, 10),
      die = rep(0.1, 10)
    ),
    tun1 = list(
      die = rep(0.15, 10)
    )
  )

  m <- generate_m_list(m_spec, tp_list)
  full_trace <- extrapolate_treatseqr(m, m_spec)
  state_names <- c("alive", "tun1", "dead")

  result <- consolidate_treatseqr_trace(full_trace, m_spec, state_names)

  dead_col <- result$t[, "dead"]
  expect_true(all(diff(dead_col) >= 0))
  expect_equal(unname(dead_col[1]), 0)
  expect_true(all(dead_col[-1] > 0))
})

test_that("consolidate_treatseqr_trace d matrix dead column is survival", {
  m_spec <- specify_m(tunnel_lengths = c(6), pre_tunnel_states = 1)

  tp_list <- list(
    alive = list(
      tun1 = rep(0.1, 6),
      die = rep(0.15, 6)
    ),
    tun1 = list(
      die = rep(0.2, 6)
    )
  )

  m <- generate_m_list(m_spec, tp_list)
  full_trace <- extrapolate_treatseqr(m, m_spec)
  state_names <- c("alive", "tun1", "dead")

  result <- consolidate_treatseqr_trace(full_trace, m_spec, state_names)

  # d matrix dead column should be 1 - dead proportion from t matrix
  # (survival, not cumulative death)
  # th = ncol(full_trace) - 1 = 6
  th <- ncol(full_trace) - 1L
  expected_survival <- 1 - result$t[1:th, "dead"]
  expect_equal(unname(result$d[, "dead"]), unname(expected_survival))
})

test_that("consolidate_treatseqr_trace initial state starts at 1", {
  m_spec <- specify_m(tunnel_lengths = c(5), pre_tunnel_states = 1)

  tp_list <- list(
    s1 = list(
      s2 = rep(0.2, 5),
      die = rep(0.1, 5)
    ),
    s2 = list(
      die = rep(0.3, 5)
    )
  )

  m <- generate_m_list(m_spec, tp_list)
  full_trace <- extrapolate_treatseqr(m, m_spec)
  state_names <- c("s1", "s2", "dead")

  result <- consolidate_treatseqr_trace(full_trace, m_spec, state_names)

  expect_equal(unname(result$t[1, "s1"]), 1)
  expect_equal(unname(result$t[1, "s2"]), 0)
  expect_equal(unname(result$t[1, "dead"]), 0)
})

test_that("consolidate_treatseqr_trace zero-pads d for short tunnel in variable-length spec", {
  # Use a hand-crafted full_trace to test zero-padding in d without going
  # through the full pipeline (which is blocked by validate_tp_source until
  # variable-length TP support is added).
  #
  # Spec: tunnel_lengths = c(2, 4), pre_tunnels = 1
  # matrix_size = 1 + 2 + 4 + 1 = 8
  # tunnel_starts = c(2, 4), dead = 8
  # th (from ncol(full_trace) - 1) = 4
  m_spec <- specify_m(tunnel_lengths = c(2, 4), pre_tunnel_states = 1)

  # Build a minimal full_trace (8 rows x 5 cols = t0:t4)
  # All population starts in state 1 and stays there for simplicity
  full_trace <- matrix(0, nrow = m_spec$matrix_size, ncol = 5)
  full_trace[1, ] <- 1   # all in pre-tunnel state throughout
  full_trace[m_spec$matrix_size, ] <- 0  # nobody dead

  state_names <- c("pre", "tun1", "tun2", "dead")

  result <- consolidate_treatseqr_trace(full_trace, m_spec, state_names)

  # th = 4
  th <- ncol(full_trace) - 1L
  expect_equal(nrow(result$d), th)

  # d matrix should contain no NA values (short tunnel was zero-padded)
  expect_false(anyNA(result$d))
  expect_true(all(result$d >= 0))
})
