test_that("generate_m_list returns correct structure and dimensions", {
  # Minimal valid model specification
  m_spec <- specify_m(tunnel_lengths = c(2, 2), pre_tunnel_states = 1)
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

  # TODO: Needs state_names
  result <- generate_m_list(m_spec, tp_list)

  expect_type(result, "list")
  expect_named(result, c("m1", "m1_dest", "m2", "state_names"))
  expect_true(is.array(result$m1))
  expect_true(methods::is(result$m2, "sparseMatrix"))
  # m1 length = th = max tunnel length = 2
  expect_equal(dim(result$m1), c(2L, 1L, 4L))
  expect_equal(dim(result$m2), c(m_spec$matrix_size, m_spec$matrix_size))
})

test_that("generate_m_list throws error for invalid transition probabilities", {
  m_spec <- specify_m(tunnel_lengths = c(2, 2), pre_tunnel_states = 1)
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
  # TODO: Needs state_names
  expect_error(generate_m_list(m_spec, tp_list), "Probabilities for state")
})

test_that("generate_m_list m2 rows sum to 1", {
  m_spec <- specify_m(tunnel_lengths = c(3), pre_tunnel_states = 1)
  tp_list <- list(
    s1 = list(s2 = rep(0.2, 3), die = rep(0.1, 3)),
    s2 = list(die = rep(0.3, 3))
  )
  names(tp_list) <- c("s1", "s2")
  # TODO: Needs state_names
  result <- generate_m_list(m_spec, tp_list)
  m2 <- result$m2
  # All rows except the first should sum to 1
  row_sums <- Matrix::rowSums(m2[2:nrow(m2), ])
  expect_true(all(abs(row_sums - 1) < 1e-10))
})

# Backwards transition tests
test_that("generate_m_list handles non-zero backwards pre-tunnel transition", {
  # 2 pre-tunnel states, 2 tunnels of length 3.
  # pre2 has a non-zero backwards transition to pre1.
  th <- 3
  m_spec <- specify_m(tunnel_lengths = rep(th, 2), pre_tunnel_states = 2)

  tp <- list(
    pre1 = list(
      pre2 = rep(0.05, th),
      tun1 = rep(0.10, th),
      tun2 = rep(0.05, th),
      die = rep(0.02, th)
    ),
    pre2 = list(
      pre1 = rep(0.03, th), # backwards pre-tunnel transition
      tun1 = rep(0.10, th),
      tun2 = rep(0.05, th),
      die = rep(0.02, th)
    ),
    tun1 = list(tun2 = rep(0.15, th), die = rep(0.10, th)),
    tun2 = list(die = rep(0.20, th))
  )

  # TODO: Needs state_names
  result <- generate_m_list(m_spec, tp, first_state_name = "pre1")

  # All m1 rows sum to 1 across all model cycles
  m1_row_sums <- apply(result$m1, c(1, 2), sum)
  expect_true(all(abs(m1_row_sums - 1) < 1e-10))

  # m2 tunnel + death rows sum to 1
  m2_row_sums <- Matrix::rowSums(result$m2[
    (m_spec$pre_tunnels + 1):nrow(result$m2),
  ])
  expect_true(all(abs(m2_row_sums - 1) < 1e-10))

  # The backwards transition from pre2 (row 2) to pre1 (col 1) is 0.03
  expect_equal(as.numeric(result$m1[1, 2, 1]), 0.03)
})

test_that("generate_m_list handles non-zero backwards tunnel transition", {
  # 1 pre-tunnel, 2 tunnels of length 3.
  # tun2 has a non-zero backwards transition to tun1.
  # tunnel_starts: tun1=col 2, tun2=col 5; death=col 8
  th <- 3
  m_spec <- specify_m(tunnel_lengths = rep(th, 2), pre_tunnel_states = 1)

  tp <- list(
    pre = list(tun1 = rep(0.10, th), tun2 = rep(0.05, th), die = rep(0.02, th)),
    tun1 = list(tun2 = rep(0.15, th), die = rep(0.10, th)),
    tun2 = list(
      tun1 = rep(0.05, th), # backwards tunnel transition
      die = rep(0.20, th)
    )
  )

  # TODO: Needs state_names
  result <- generate_m_list(m_spec, tp, first_state_name = "pre")

  # m2 tunnel + death rows sum to 1
  m2_row_sums <- Matrix::rowSums(result$m2[
    (m_spec$pre_tunnels + 1):nrow(result$m2),
  ])
  expect_true(all(abs(m2_row_sums - 1) < 1e-10))

  # tun2 occupies rows 5:7. Each should have a non-zero entry at col 2 (tun1 start)
  # with value 0.05 (the backwards transition probability)
  expect_equal(as.numeric(result$m2[5, 2]), 0.05)
  expect_equal(as.numeric(result$m2[6, 2]), 0.05)
  expect_equal(as.numeric(result$m2[7, 2]), 0.05)
})

test_that("generate_m_list rejects reordered state_names, places probs positionally", {
  # 1 pre-tunnel, 2 tunnels of length 2.
  # tunnel_starts: tun1=col 2, tun2=col 4; death=col 6
  th <- 2
  m_spec <- specify_m(tunnel_lengths = rep(th, 2), pre_tunnel_states = 1)

  tp <- list(
    pre = list(tun1 = rep(0.11, th), tun2 = rep(0.22, th), die = rep(0.03, th)),
    tun1 = list(tun2 = rep(0.44, th), die = rep(0.04, th)),
    tun2 = list(tun1 = rep(0.55, th), die = rep(0.05, th))
  )

  # Matching order: probabilities land in the columns implied by specify_m()
  result <- generate_m_list(
    m_spec,
    tp,
    first_state_name = "pre",
    state_names = c("pre", "tun1", "tun2", "die")
  )
  # m1's 3rd dimension is in canonical destination order. m1_dest maps those
  # slots onto columns of the full matrix M.
  expect_equal(result$m1_dest, c(1, 2, 4, 6))
  m1_cycle1 <- result$m1[1, 1, ]
  expect_equal(m1_cycle1[match(2, result$m1_dest)], 0.11) # pre -> tun1
  expect_equal(m1_cycle1[match(4, result$m1_dest)], 0.22) # pre -> tun2
  expect_equal(m1_cycle1[match(6, result$m1_dest)], 0.03) # pre -> die

  # Reordered state_names must error, not silently swap the tunnel columns
  expect_error(
    generate_m_list(
      m_spec,
      tp,
      first_state_name = "pre",
      state_names = c("pre", "tun2", "tun1", "die")
    ),
    "'state_names' must match names\\(tp_source\\) in the same order"
  )
})

# This one guards against misordering of the TPs!
test_that("m1_dest matches the sorted j coordinates from specify_m", {
  # 3 pre-tunnel states, 2 tunnels of unequal length.
  # tunnel_starts: tun1 = col 4, tun2 = col 7; death = col 11
  th <- 5
  m_spec <- specify_m(tunnel_lengths = c(3, 4), pre_tunnel_states = 3)

  tp <- list(
    pre1 = list(
      pre2 = rep(0.10, th),
      pre3 = NULL,
      tun1 = rep(0.05, th),
      tun2 = rep(0.02, th),
      die = rep(0.01, th)
    ),
    pre2 = list(
      pre1 = rep(0.03, th), # backwards pre-tunnel
      pre3 = rep(0.20, th),
      tun1 = rep(0.06, th),
      tun2 = rep(0.02, th),
      die = rep(0.02, th)
    ),
    pre3 = list(
      pre1 = NULL,
      pre2 = rep(0.04, th), # backwards pre-tunnel
      tun1 = rep(0.08, th),
      tun2 = rep(0.03, th),
      die = rep(0.02, th)
    ),
    tun1 = list(tun2 = rep(0.15, 3), die = rep(0.10, 3)),
    tun2 = list(tun1 = rep(0.05, 4), die = rep(0.20, 4))
  )

  # TODO: Needs state_names
  result <- generate_m_list(m_spec, tp, first_state_name = "pre1")

  expect_equal(dim(result$m1), c(th, 3, 6))
  expect_equal(result$m1_dest, sort(unique(m_spec$m1_ijx[, "j"])))
  expect_equal(result$m1_dest, c(1, 2, 3, 4, 7, 11))

  # every (cycle, pre-state) slice is a complete probability distribution
  expect_true(all(abs(apply(result$m1, c(1, 2), sum) - 1) < 1e-10))

  # specify_m() hoists each state's self-transition to the front of its j
  # coordinates, so pre3's raw order is (3, 1, 2, 4, 7, 11). The array must be
  # in canonical order instead: slot 2 is pre2 and slot 3 is pre3's p_stay.
  # If the hoisted order leaked through, slot 1 would hold p_stay and these
  # would swap -- and the row would still sum to 1, so only this catches it.
  expect_equal(result$m1[1, 3, 2], 0.04) # pre3 -> pre2
  expect_equal(result$m1[1, 3, 3], 0.83) # pre3 p_stay = 1 - 0.17
  expect_equal(result$m1[1, 2, 1], 0.03) # pre2 -> pre1

  # and it extrapolates to a conserved cohort
  trace <- extrapolate_treatseqr(result, m_spec)
  expect_true(all(abs(colSums(trace) - 1) < 1e-10))
})
