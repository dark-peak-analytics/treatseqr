# Test the specify_m() function for correct output structure and values.
# - Checks that the return type is a list with correct elements.
# - Verifies that the m1_ijx matrix contains the correct indices for the top row(s) of
# the transition matrix.
# - Ensures that the matrix_size is correct for the given tunnel_lengths vector.
# - Confirms that invalid inputs trigger errors.

test_that("specify_m returns correct structure", {
  result <- specify_m(tunnel_lengths = rep(5, 5))

  expect_type(result, "list")
  expect_setequal(
    names(result),
    c("tunnel_lengths", "m1_ijx", "m2_ijx", "matrix_size", "pre_tunnels")
  )
  expect_equal(result$tunnel_lengths, rep(5, 5))
  expect_equal(result$pre_tunnels, 1)
  expect_true(is.matrix(result$m1_ijx))
  expect_true(is.matrix(result$m2_ijx))
  expect_equal(ncol(result$m1_ijx), 3)
  expect_equal(ncol(result$m2_ijx), 3)
  expect_equal(colnames(result$m1_ijx), c("i", "j", "x"))
  expect_equal(colnames(result$m2_ijx), c("i", "j", "x"))
})

test_that("specify_m returns correct values for known case", {
  result <- specify_m(tunnel_lengths = rep(5, 5))

  expect_equal(
    result$m1_ijx,
    structure(
      c(1, 1, 1, 1, 1, 1, 1, 1, 2, 7, 12, 17, 22, 27, 0, 0, 0, 0, 0, 0, 0),
      dim = c(7L, 3L),
      dimnames = list(NULL, c("i", "j", "x"))
    )
  )
  expect_equal(result$matrix_size, 27)
})

test_that("specify_m handles minimal valid input", {
  result <- specify_m(tunnel_lengths = c(1))

  expect_equal(
    result$m1_ijx,
    structure(
      c(1, 1, 1, 1, 2, 3, 0, 0, 0),
      dim = c(3L, 3L),
      dimnames = list(NULL, c("i", "j", "x"))
    )
  )
  expect_equal(
    result$m2_ijx,
    structure(
      c(2, 2, 3, 2, 3, 3, 0, 0, 1),
      dim = c(3L, 3L),
      dimnames = list(NULL, c("i", "j", "x"))
    )
  )
  expect_equal(result$matrix_size, 3)
})

test_that("specify_m variable-length tunnels: matrix_size and tunnel_starts are correct", {
  # tunnel_lengths = c(3, 5, 2): 3 tunnels of different lengths
  # pre_tunnel_states = 1 (default)
  # tunnel_starts: 1 + cumsum(c(1, 3, 5)) = 1 + c(1, 4, 9) = c(2, 5, 10)
  # dead = max(c(2,5,10)) + 2 = 12
  # matrix_size = 12
  result <- specify_m(tunnel_lengths = c(3, 5, 2))

  expect_equal(result$tunnel_lengths, c(3, 5, 2))
  expect_equal(result$matrix_size, 12)

  # All indices within bounds
  expect_true(all(result$m1_ijx[, "i"] <= result$matrix_size))
  expect_true(all(result$m1_ijx[, "j"] <= result$matrix_size))
  expect_true(all(result$m2_ijx[, "i"] <= result$matrix_size))
  expect_true(all(result$m2_ijx[, "j"] <= result$matrix_size))
})

test_that("specify_m errors on non-integer input", {
  expect_error(
    specify_m(tunnel_lengths = c(2.5, 2)),
    "tunnel_lengths must be a vector of positive integers"
  )
})

test_that("specify_m errors on negative input", {
  expect_error(
    specify_m(tunnel_lengths = c(-1, 2)),
    "tunnel_lengths must be a vector of positive integers"
  )
})

test_that("specify_m errors on zero input", {
  expect_error(
    specify_m(tunnel_lengths = c(0, 2)),
    "tunnel_lengths must be a vector of positive integers"
  )
})

test_that("specify_m errors on empty vector", {
  expect_error(
    specify_m(tunnel_lengths = integer(0)),
    "tunnel_lengths must be a vector of positive integers"
  )
})

test_that("specify_m works for multiple pre-tunnel states", {
  result <- specify_m(tunnel_lengths = rep(10, 3), pre_tunnel_states = 10)

  expect_equal(result$pre_tunnels, 10)
  expect_equal(nrow(result$m1_ijx), 95)

  # The last pre-tunnel state (row 10) should not have any pre-tunnel
  # transitions except to itself and to tunnels/death
  last_pre_state_transitions <- result$m1_ijx[result$m1_ijx[, "i"] == 10, "j"]

  expect_true(
    all(last_pre_state_transitions >= 10),
    label = paste0(
      "The last pre-tunnel state should not have any transitions to ",
      "previous pre-tunnel states (only to itself, tunnels, or death)."
    )
  )
})

test_that("specify_m works for large valid input", {
  result <- specify_m(tunnel_lengths = rep(100, 10))

  expect_equal(nrow(result$m1_ijx), 12)
  expect_equal(
    result$m1_ijx,
    structure(
      c(
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
        1, 2, 102, 202, 302, 402, 502, 602, 702, 802, 902, 1002,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
      ),
      dim = c(12L, 3L),
      dimnames = list(NULL, c("i", "j", "x"))
    )
  )
  expect_equal(result$matrix_size, 1002)
})

test_that("specify_m errors if required argument tunnel_lengths is missing", {
  expect_error(
    specify_m(),
    "argument \"tunnel_lengths\" is missing"
  )
})

test_that("specify_m pre_tunnel_states parameter works correctly", {
  # Test default (1 pre-tunnel state)
  result1 <- specify_m(tunnel_lengths = rep(5, 2))
  expect_equal(result1$pre_tunnels, 1)

  # Test explicit pre_tunnel_states
  result2 <- specify_m(tunnel_lengths = rep(5, 2), pre_tunnel_states = 3)
  expect_equal(result2$pre_tunnels, 3)

  # m1_ijx should have more rows with more pre-tunnel states
  expect_true(nrow(result2$m1_ijx) > nrow(result1$m1_ijx))
})

test_that("specify_m matrix dimensions are consistent", {
  result <- specify_m(tunnel_lengths = rep(10, 4), pre_tunnel_states = 2)

  expect_true(all(result$m1_ijx[, "i"] <= result$matrix_size))
  expect_true(all(result$m1_ijx[, "j"] <= result$matrix_size))
  expect_true(all(result$m2_ijx[, "i"] <= result$matrix_size))
  expect_true(all(result$m2_ijx[, "j"] <= result$matrix_size))

  expect_true(all(result$m1_ijx[, "i"] > 0))
  expect_true(all(result$m1_ijx[, "j"] > 0))
  expect_true(all(result$m2_ijx[, "i"] > 0))
  expect_true(all(result$m2_ijx[, "j"] > 0))
})

test_that("specify_m m2_ijx death state is correctly set", {
  result <- specify_m(tunnel_lengths = rep(5, 3))

  last_row <- result$m2_ijx[nrow(result$m2_ijx), ]

  expect_equal(unname(last_row["i"]), result$matrix_size)
  expect_equal(unname(last_row["j"]), result$matrix_size)
  expect_equal(unname(last_row["x"]), 1)
})

test_that("specify_m errors on invalid pre_tunnel_states", {
  expect_error(
    specify_m(tunnel_lengths = rep(5, 2), pre_tunnel_states = 0),
    "'pre_tunnel_states' must be a single positive integer"
  )
  expect_error(
    specify_m(tunnel_lengths = rep(5, 2), pre_tunnel_states = -1),
    "'pre_tunnel_states' must be a single positive integer"
  )
  expect_error(
    specify_m(tunnel_lengths = rep(5, 2), pre_tunnel_states = 2.5),
    "'pre_tunnel_states' must be a single positive integer"
  )
  expect_error(
    specify_m(tunnel_lengths = rep(5, 2), pre_tunnel_states = 11),
    "'pre_tunnel_states' must be a single positive integer. Limit of 10"
  )
})
