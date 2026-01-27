#' Validate and Sanitize Transition Probability Source List
#'
#' Validates that the transition probability list (`tp_source`) follows the
#' strict structural rules required for sparse matrix extrapolation. It checks
#' for reducing lengths, topological consistency of state names, consistent
#' time horizons, and valid probability values. It also sanitizes the list by
#' replacing `NULL` entries with zero vectors.
#'
#' @param tp_source A named list of lists containing transition probability
#' vectors. The structure must follow specific rules:
#'   \itemize{
#'     \item The list lengths must reduce by 1 sequentially (e.g., 5, 4, 3, 2,
#'     1).
#'     \item The names of the inner lists must be topologically consistent
#'     (subset rule).
#'     \item All probability vectors must have the same length (number of
#'     cycles).
#'     \item Probabilities must be non-negative and sum to <= 1 (row-wise).
#'   }
#'
#' @return The sanitized `tp_source` list where all `NULL` entries have been
#'   replaced with vectors of zeros.
#' @export
validate_tp_source <- function(tp_source) {
  assertthat::assert_that(
    is.list(tp_source),
    !is.null(names(tp_source)),
    msg = "tp_source must be a named list"
  )

  n_states <- length(tp_source)
  assertthat::assert_that(n_states > 0, msg = "tp_source cannot be empty")

  # 1. Length validation (reducing by 1)
  element_lengths <- vapply(tp_source, length, numeric(1))

  # The first element dictates the starting size. The sequence must go down to
  # 1.
  expected_lengths <- seq(
    from = element_lengths[1],
    by = -1,
    length.out = n_states
  )

  assertthat::assert_that(
    element_lengths[n_states] == 1,
    msg = "The last element of tp_source must have length 1"
  )

  assertthat::assert_that(
    all(element_lengths == expected_lengths),
    msg = sprintf(
      "Lengths of tp_source elements must reduce by 1 sequentially. Found: %s. Expected: %s.",
      paste(element_lengths, collapse = ", "),
      paste(expected_lengths, collapse = ", ")
    )
  )

  # 2. Topological Naming Validation ("Backwards" check)
  # The keys of state i+1 must be a subset of the keys of state i.
  if (n_states > 1) {
    for (i in 1:(n_states - 1)) {
      curr_names <- names(tp_source[[i]])
      next_names <- names(tp_source[[i + 1]])

      assertthat::assert_that(
        !is.null(curr_names) && !is.null(next_names),
        msg = "All inner lists in tp_source must be named"
      )

      missing_names <- setdiff(next_names, curr_names)
      assertthat::assert_that(
        length(missing_names) == 0,
        msg = sprintf(
          "Topology error: State '%s' contains transitions to {%s} which are not present in the previous state '%s'.",
          names(tp_source)[i + 1],
          paste(missing_names, collapse = ", "),
          names(tp_source)[i]
        )
      )
    }
  }

  # 3. Cycle Detection & Sanitization
  # Find the number of cycles from the first non-null, non-empty vector
  n_cyc <- NULL

  # Flatten temporarily to find a valid vector
  all_vectors <- unlist(tp_source, recursive = FALSE)
  for (v in all_vectors) {
    if (!is.null(v) && length(v) > 0) {
      n_cyc <- length(v)
      break
    }
  }

  assertthat::assert_that(
    !is.null(n_cyc),
    msg = "Could not determine number of cycles (all vectors are NULL or empty)"
  )

  # Sanitize: Loop through structure, validate lengths/values, replace NULLs
  sanitized_tp <- lapply(names(tp_source), function(state_name) {
    state_list <- tp_source[[state_name]]

    # Process each transition within the state
    new_state_list <- lapply(names(state_list), function(trans_name) {
      vec <- state_list[[trans_name]]

      if (is.null(vec)) {
        return(numeric(n_cyc))
      }

      assertthat::assert_that(
        length(vec) == n_cyc,
        msg = sprintf(
          "Length mismatch in transition %s -> %s. Expected %d, found %d.",
          state_name,
          trans_name,
          n_cyc,
          length(vec)
        )
      )

      assertthat::assert_that(
        is.numeric(vec),
        all(vec >= 0),
        msg = sprintf(
          "Invalid values in transition %s -> %s. Must be numeric and >= 0.",
          state_name,
          trans_name
        )
      )

      return(vec)
    })
    names(new_state_list) <- names(state_list)
    return(new_state_list)
  })
  names(sanitized_tp) <- names(tp_source)

  # 4. Summation Constraint
  # Element-wise sum of all transitions out of a state must be <= 1
  for (state_name in names(sanitized_tp)) {
    trans_list <- sanitized_tp[[state_name]]

    # Efficient element-wise sum of list of vectors
    total_prob <- Reduce(`+`, trans_list)

    # Check with tolerance for floating point arithmetic
    assertthat::assert_that(
      all(total_prob <= 1 + 1e-10),
      msg = sprintf(
        "Probabilities for state '%s' sum to > 1 in some cycles.",
        state_name
      )
    )
  }

  return(sanitized_tp)
}
