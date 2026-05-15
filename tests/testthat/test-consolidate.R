test_that("consolidate_treatseqr_trace returns correct structure", {
  # Create a minimal model specification
  m_spec <- specify_m(n_cycles = 5, n_tunnels = 2, pre_tunnel_states = 1)

  # Create transition probabilities
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

  # Generate matrices and run extrapolation
  m <- generate_m_list(m_spec, tp_list)
  full_trace <- extrapolate_treatseqr(m, m_spec)

  # Define state names
  state_names <- c("state1", "state2", "state3", "dead")

  # Run consolidation
  result <- consolidate_treatseqr_trace(full_trace, m_spec, state_names)

  # Check structure
  expect_type(result, "list")
  expect_named(result, c("t", "d"))
  expect_true(is.matrix(result$t))
  expect_true(is.matrix(result$d))
})

test_that("consolidate_treatseqr_trace has correct dimensions", {
  m_spec <- specify_m(n_cycles = 10, n_tunnels = 3, pre_tunnel_states = 1)

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

  # t matrix should have n_cycles + 1 rows (includes time 0)
  expect_equal(nrow(result$t), m_spec$n_cycles + 1)
  # d matrix should have n_cycles rows (no time 0)
  expect_equal(nrow(result$d), m_spec$n_cycles)
  # Both should have correct number of columns (states)
  expect_equal(ncol(result$t), length(state_names))
  expect_equal(ncol(result$d), length(state_names))
})

test_that("consolidate_treatseqr_trace column names are correct", {
  m_spec <- specify_m(n_cycles = 3, n_tunnels = 1, pre_tunnel_states = 1)

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
  m_spec <- specify_m(n_cycles = 8, n_tunnels = 2, pre_tunnel_states = 1)

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

  # All rows in t matrix should sum to 1
  row_sums <- rowSums(result$t)
  expect_true(all(abs(row_sums - 1) < 1e-10))
})

test_that("consolidate_treatseqr_trace d matrix excludes time 0", {
  m_spec <- specify_m(n_cycles = 5, n_tunnels = 2, pre_tunnel_states = 1)

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

  # d matrix is based on within-tunnel time, so it should be calculated differently
  # than the t matrix. Just check it has valid probabilities
  expect_true(all(result$d >= 0))
  expect_true(all(result$d <= 1))
})

test_that("consolidate_treatseqr_trace handles different numbers of tunnels", {
  # Test with 3 tunnels instead of multiple pre-tunnel states
  m_spec <- specify_m(n_cycles = 4, n_tunnels = 3, pre_tunnel_states = 1)

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

  # Check dimensions are correct with 3 tunnels
  expect_equal(ncol(result$t), 5)
  expect_equal(ncol(result$d), 5)
  expect_equal(nrow(result$t), m_spec$n_cycles + 1)
  expect_equal(nrow(result$d), m_spec$n_cycles)
})

test_that("consolidate_treatseqr_trace dead state accumulates over time", {
  m_spec <- specify_m(n_cycles = 10, n_tunnels = 1, pre_tunnel_states = 1)

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

  # Dead state should be monotonically increasing
  dead_col <- result$t[, "dead"]
  expect_true(all(diff(dead_col) >= 0))

  # Dead state at time 0 should be 0
  expect_equal(unname(dead_col[1]), 0)

  # Dead state should be positive after time 0
  expect_true(all(dead_col[-1] > 0))
})

test_that("consolidate_treatseqr_trace d matrix dead column is survival", {
  m_spec <- specify_m(n_cycles = 6, n_tunnels = 1, pre_tunnel_states = 1)

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
  # trace_dead[-(th + 1)] means exclude the LAST element, not the first
  # So it's survival for times 1:n_cycles (excluding time 0)
  expected_survival <- 1 - result$t[1:m_spec$n_cycles, "dead"]
  expect_equal(unname(result$d[, "dead"]), unname(expected_survival))
})

test_that("consolidate_treatseqr_trace initial state starts at 1", {
  m_spec <- specify_m(n_cycles = 5, n_tunnels = 1, pre_tunnel_states = 1)

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

  # At time 0, all cohort should be in first state
  expect_equal(unname(result$t[1, "s1"]), 1)
  expect_equal(unname(result$t[1, "s2"]), 0)
  expect_equal(unname(result$t[1, "dead"]), 0)
})
