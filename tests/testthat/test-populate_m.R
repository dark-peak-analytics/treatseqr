test_that("generate_m_list returns correct structure and dimensions", {
  # Minimal valid model specification
  m_spec <- specify_m(n_cycles = 2, n_tunnels = 2, pre_tunnel_states = 1)
  tp_list <- list(
    state1 = list(
      state2 = rep(0.3, 2),
      state3 = rep(0.2, 2),
      die = rep(0.1, 2)
    ),
    state2 = list(state3 = rep(0.4, 2), die = rep(0.2, 2)),
    state3 = list(die = rep(0.5, 2))
  )

  # generate_m_list expects a validated structure, so names must match m_spec
  names(tp_list) <- c("state1", "state2", "state3")

  result <- generate_m_list(m_spec, tp_list)

  expect_type(result, "list")
  expect_named(result, c("m1", "m2"))
  expect_true(is.list(result$m1))
  expect_true(methods::is(result$m2, "sparseMatrix"))
  expect_length(result$m1, m_spec$n_cycles)
  expect_equal(dim(result$m2), c(m_spec$matrix_size, m_spec$matrix_size))
})

test_that("generate_m_list throws error for invalid transition probabilities", {
  m_spec <- specify_m(n_cycles = 2, n_tunnels = 2, pre_tunnel_states = 1)
  # Invalid: probabilities sum > 1
  tp_list <- list(
    state1 = list(
      state2 = rep(0.8, 2),
      state3 = rep(0.5, 2),
      die = rep(0.1, 2)
    ),
    state2 = list(state3 = rep(0.4, 2), die = rep(0.2, 2)),
    state3 = list(die = rep(0.5, 2))
  )
  names(tp_list) <- c("state1", "state2", "state3")
  expect_error(generate_m_list(m_spec, tp_list), "Probabilities for state")
})

test_that("generate_m_list m2 rows sum to 1", {
  m_spec <- specify_m(n_cycles = 3, n_tunnels = 1, pre_tunnel_states = 1)
  tp_list <- list(
    s1 = list(s2 = rep(0.2, 3), die = rep(0.1, 3)),
    s2 = list(die = rep(0.3, 3))
  )
  names(tp_list) <- c("s1", "s2")
  result <- generate_m_list(m_spec, tp_list)
  m2 <- result$m2
  # All rows except the first should sum to 1
  row_sums <- Matrix::rowSums(m2[2:nrow(m2), ])
  expect_true(all(abs(row_sums - 1) < 1e-10))
})
