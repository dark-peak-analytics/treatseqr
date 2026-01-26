# Test the m_size() function for correct output structure and values.
# - Checks that the return type is a list of length 2.
# - Verifies that the m1 vector contains the correct indices for the top row of
# the transition matrix.
# - Ensures that the matrix_size is correct for the given number of cycles and
# tunnel states.
# - Confirms that invalid input (n_cycles = 0) triggers an error.
test_that("m_size returns correct values", {
  cols <- m_size(n_cycles = 5, n_tunnels = 5)
  expect_type(cols, "list")
  expect_equal(length(cols), 2)
  expect_equal(cols$m1, c(1, 2, 7, 12, 17, 22))
  expect_equal(cols$matrix_size, 22)
  expect_error(m_size(n_cycles = 0))
})

test_that("m_size handles minimal valid input", {
  cols <- m_size(n_cycles = 1, n_tunnels = 1)
  expect_equal(cols$m1, c(1, 2))
  expect_equal(cols$matrix_size, 2)
})

test_that("m_size errors on non-integer input", {
  expect_error(
    m_size(n_cycles = 2.5, n_tunnels = 2),
    "n_cycles must be a single positive integer"
  )
  expect_error(
    m_size(n_cycles = 2, n_tunnels = 2.5),
    "n_tunnels must be a single positive integer"
  )
})

test_that("m_size errors on negative input", {
  expect_error(
    m_size(n_cycles = -1, n_tunnels = 2),
    "n_cycles must be a single positive integer"
  )
  expect_error(
    m_size(n_cycles = 2, n_tunnels = -2),
    "n_tunnels must be a single positive integer"
  )
})

test_that("m_size errors on vector input", {
  expect_error(
    m_size(n_cycles = c(2, 3), n_tunnels = 2),
    "n_cycles must be a single positive integer"
  )
  expect_error(
    m_size(n_cycles = 2, n_tunnels = c(2, 3)),
    "n_tunnels must be a single positive integer"
  )
})

test_that("m_size works for large valid input", {
  cols <- m_size(n_cycles = 100, n_tunnels = 10)
  expect_equal(length(cols$m1), 11)
  expect_equal(cols$m1, c(1, seq(2, (10 * 100) + 1, by = 100)))
  expect_equal(cols$matrix_size, 902)
})

test_that("m_size errors if arguments are missing", {
  expect_error(m_size(n_cycles = 2), "argument \"n_tunnels\" is missing")
  expect_error(m_size(n_tunnels = 2), "argument \"n_cycles\" is missing")
})
