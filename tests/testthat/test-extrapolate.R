# Helper: build a minimal valid tp_list for a given spec (uniform tunnel lengths only)
make_tp_list <- function(spec, p_move = 0.1, p_die = 0.05) {
  n_tun <- length(spec$tunnel_lengths)
  n_states <- spec$pre_tunnels + as.integer(n_tun)
  nms <- paste0("s", seq_len(n_states))
  th <- spec$tunnel_lengths[1]

  lapply(seq_len(n_states), function(i) {
    destinations <- if (i < n_states) nms[(i + 1):n_states] else character(0)
    out <- lapply(seq_along(destinations), function(j) rep(p_move, th))
    out <- c(out, list(rep(p_die, th)))
    names(out) <- c(destinations, "die")
    out
  }) |>
    setNames(nms)
}

# ------------------------------------------------------------
# Cohort conservation
# ------------------------------------------------------------

test_that("cohort is conserved at all time points — 1 tunnel", {
  spec <- specify_m(tunnel_lengths = c(10), pre_tunnel_states = 1)
  tp <- list(
    alive = list(
      tun1 = rep(0.05, 10),
      die = rep(0.10, 10)
    ),
    tun1 = list(
      die = rep(0.15, 10)
    )
  )
  # TODO: Needs state_names
  m <- generate_m_list(spec, tp)
  pop <- extrapolate_treatseqr(m, spec)
  sums <- colSums(pop)
  expect_true(
    all(abs(sums - 1) < 1e-10),
    label = paste(
      "colSums deviate from 1:",
      paste(round(sums, 12), collapse = ", ")
    )
  )
})

test_that("cohort is conserved at all time points — 2 tunnels", {
  spec <- specify_m(tunnel_lengths = rep(8, 2), pre_tunnel_states = 1)
  tp <- list(
    s1 = list(s2 = rep(0.15, 8), s3 = rep(0.10, 8), die = rep(0.05, 8)),
    s2 = list(s3 = rep(0.20, 8), die = rep(0.10, 8)),
    s3 = list(die = rep(0.25, 8))
  )
  # TODO: Needs state_names
  m <- generate_m_list(spec, tp)
  pop <- extrapolate_treatseqr(m, spec)
  sums <- colSums(pop)
  expect_true(
    all(abs(sums - 1) < 1e-10),
    label = paste(
      "colSums deviate from 1:",
      paste(round(sums, 12), collapse = ", ")
    )
  )
})

test_that("cohort is conserved — short horizon, high p_stay stresses last tunnel state", {
  # tunnel_lengths = c(4, 4): patients can reach end of last tunnel within horizon.
  # High p_die is deliberately low so p_stay is large, ensuring self-loop carries
  # significant population.
  spec <- specify_m(tunnel_lengths = rep(4, 2), pre_tunnel_states = 1)
  tp <- list(
    s1 = list(s2 = rep(0.30, 4), s3 = rep(0.20, 4), die = rep(0.01, 4)),
    s2 = list(s3 = rep(0.30, 4), die = rep(0.01, 4)),
    s3 = list(die = rep(0.01, 4))
  )
  # TODO: Needs state_names
  m <- generate_m_list(spec, tp)
  pop <- extrapolate_treatseqr(m, spec)
  sums <- colSums(pop)
  expect_true(
    all(abs(sums - 1) < 1e-10),
    label = paste(
      "colSums deviate from 1:",
      paste(round(sums, 12), collapse = ", ")
    )
  )
})

test_that("cohort is conserved — 3 tunnels", {
  spec <- specify_m(tunnel_lengths = rep(6, 3), pre_tunnel_states = 1)
  tp <- make_tp_list(spec, p_move = 0.08, p_die = 0.04)
  # TODO: Needs state_names
  m <- generate_m_list(spec, tp)
  pop <- extrapolate_treatseqr(m, spec)
  sums <- colSums(pop)
  expect_true(
    all(abs(sums - 1) < 1e-10),
    label = paste(
      "colSums deviate from 1:",
      paste(round(sums, 12), collapse = ", ")
    )
  )
})

# ------------------------------------------------------------
# matrix_size correctness (no overflow row)
# ------------------------------------------------------------

test_that("matrix_size equals pre_tunnels + sum(tunnel_lengths) + 1", {
  check <- function(pre, nt, nc) {
    spec <- specify_m(tunnel_lengths = rep(nc, nt), pre_tunnel_states = pre)
    expected <- pre + nt * nc + 1
    expect_equal(
      spec$matrix_size,
      expected,
      label = sprintf("pre=%d, n_tunnels=%d, n_cycles=%d", pre, nt, nc)
    )
  }
  check(1, 1, 5)
  check(1, 2, 5)
  check(1, 3, 10)
  check(2, 2, 8)
  check(1, 1, 1)
})

# ------------------------------------------------------------
# Last tunnel state self-loop
# ------------------------------------------------------------

test_that("last state of last tunnel self-loops in m2", {
  spec <- specify_m(tunnel_lengths = rep(3, 2), pre_tunnel_states = 1)
  # tunnel_lengths = c(3, 3), pre_tunnels = 1
  # tunnel_starts = 1 + cumsum(c(1, 3)) = c(2, 5)
  # last state of last tunnel = tunnel_starts[2] + tunnel_lengths[2] - 1 = 5 + 2 = 7
  # dead = max(5) + 3 = 8
  last_tun_state <- spec$pre_tunnels + sum(spec$tunnel_lengths)

  tp <- list(
    s1 = list(s2 = rep(0.1, 3), s3 = rep(0.1, 3), die = rep(0.05, 3)),
    s2 = list(s3 = rep(0.2, 3), die = rep(0.10, 3)),
    s3 = list(die = rep(0.20, 3))
  )
  # TODO: Needs state_names
  m <- generate_m_list(spec, tp)

  # The self-loop diagonal element must be non-zero (p_stay = 1 - 0.20 = 0.80)
  expect_gt(m$m2[last_tun_state, last_tun_state], 0)

  # No overflow row — matrix is exactly matrix_size x matrix_size
  expect_equal(nrow(m$m2), spec$matrix_size)
  expect_equal(ncol(m$m2), spec$matrix_size)
})

test_that("last tunnel state population persists when forced p_stay is high", {
  # All patients enter tun1 immediately, die rate is very low.
  # With tunnel_lengths = c(3), by cycle 3 most patients are at the last
  # tunnel state. Without self-loop they would be lost; with self-loop they
  # accumulate in the last state.
  spec <- specify_m(tunnel_lengths = c(3), pre_tunnel_states = 1)
  tp <- list(
    alive = list(tun1 = rep(0.90, 3), die = rep(0.005, 3)),
    tun1 = list(die = rep(0.005, 3))
  )
  # TODO: Needs state_names
  m <- generate_m_list(spec, tp)
  pop <- extrapolate_treatseqr(m, spec)

  # last tunnel state = pre_tunnels + tunnel_lengths[1] = 1 + 3 = 4 (not dead)
  # dead = matrix_size = 5
  last_tun_state <- spec$pre_tunnels + spec$tunnel_lengths[1]
  dead_state <- spec$matrix_size

  # By end of horizon, the last tunnel state should hold non-trivial population
  expect_gt(pop[last_tun_state, ncol(pop)], 0)

  # Cohort must still be conserved
  expect_true(all(abs(colSums(pop) - 1) < 1e-10))
})

# ------------------------------------------------------------
# Multiple pre-tunnel states
# ------------------------------------------------------------

test_that("extrapolate works with pre_tunnel_states = 2", {
  spec <- specify_m(tunnel_lengths = rep(5, 2), pre_tunnel_states = 2)

  # With 2 pre-tunnel states the tp_list must have reducing lengths starting
  # from 4: state1 -> 4 destinations, state2 -> 3, tun1 -> 2, tun2 -> 1
  tp <- list(
    state1 = list(
      state2 = rep(0.10, 5),
      tun1 = rep(0.10, 5),
      tun2 = rep(0.05, 5),
      die = rep(0.05, 5)
    ),
    state2 = list(
      tun1 = rep(0.15, 5),
      tun2 = rep(0.05, 5),
      die = rep(0.05, 5)
    ),
    tun1 = list(
      tun2 = rep(0.10, 5),
      die = rep(0.10, 5)
    ),
    tun2 = list(
      die = rep(0.20, 5)
    )
  )

  # TODO: Needs state_names
  m <- generate_m_list(spec, tp)
  pop <- extrapolate_treatseqr(m, spec)

  # dimensions: nrow = matrix_size, ncol = th + 1 (th = tunnel_lengths[1] = 5)
  expect_equal(nrow(pop), spec$matrix_size)
  expect_equal(ncol(pop), spec$tunnel_lengths[1] + 1)

  # cohort conservation
  expect_true(all(abs(colSums(pop) - 1) < 1e-10))

  # initial state: all population in state 1 at t=0
  expect_equal(pop[1, 1], 1)
  expect_true(all(pop[-1, 1] == 0))
})

test_that("cohort is conserved with pre_tunnel_states = 2, short horizon", {
  spec <- specify_m(tunnel_lengths = rep(4, 2), pre_tunnel_states = 2)
  tp <- list(
    s1 = list(
      s2 = rep(0.20, 4),
      t1 = rep(0.20, 4),
      t2 = rep(0.10, 4),
      die = rep(0.02, 4)
    ),
    s2 = list(t1 = rep(0.25, 4), t2 = rep(0.10, 4), die = rep(0.02, 4)),
    t1 = list(t2 = rep(0.25, 4), die = rep(0.02, 4)),
    t2 = list(die = rep(0.02, 4))
  )
  # TODO: Needs state_names
  m <- generate_m_list(spec, tp)
  pop <- extrapolate_treatseqr(m, spec)
  expect_true(all(abs(colSums(pop) - 1) < 1e-10))
})
