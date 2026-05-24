#' Validate and Sanitize Transition Probability List for Sequence Modelling
#'
#' Validates that the transition probability list (`tp_source`) follows the
#' strict structural rules required for sparse matrix extrapolation. It checks
#' for reducing lengths, topological consistency of state names, consistent
#' time horizons, and valid probability values. It also sanitizes the list by
#' replacing `NULL` entries with zero vectors of the correct length.
#'
#' Specifically, this function ensures:
#'   \itemize{
#'     \item The list lengths must reduce by 1 sequentially (e.g., 5, 4, 3, 2,
#'     1).
#'     \item The names of the inner lists must be topologically consistent
#'     (subset rule).
#'     \item All probability vectors must have consistent lengths (see
#'     \code{tunnel_lengths} below).
#'     \item Probabilities must be non-negative and sum to <= 1 (row-wise).
#'     \item Any `NULL` transition vectors are replaced with numeric vectors of
#'     zeros of the correct length.
#'   }
#'
#' @param tp_source A named list of lists containing transition probability
#'   vectors. The structure must follow specific rules as described above.
#' @param tunnel_lengths An optional integer vector giving the length of each
#'   tunnel state block, in order (one element per tunnel). When \code{NULL}
#'   (default), all probability vectors must have the same length and the
#'   function behaves as in the uniform-tunnel case. When supplied, the number
#'   of pre-tunnel states is inferred as
#'   \code{length(tp_source) - length(tunnel_lengths)}, pre-tunnel TP vectors
#'   must all share a common length, and each tunnel \code{k}'s vectors must
#'   have length \code{tunnel_lengths[k]}.
#'
#' @return The sanitized `tp_source` list where all `NULL` entries have been
#'   replaced with vectors of zeros of the correct length.
#' @export
validate_tp_source <- function(tp_source, tunnel_lengths = NULL) {
  assertthat::assert_that(
    is.list(tp_source),
    !is.null(names(tp_source)),
    msg = "tp_source must be a named list"
  )

  n_states <- length(tp_source)
  assertthat::assert_that(
    n_states > 1,
    msg = paste0(
      "'tp_source' must have length greater than 1 otherwise there are no",
      " tunnel states. The death state should not be included as there is",
      " nowhere to transition to from dead"
    )
  )

  # Length validation (reducing n_states:1)
  element_lengths <- vapply(tp_source, length, numeric(1))
  expected_lengths <- seq(from = n_states, by = -1, length.out = n_states)
  assertthat::assert_that(
    all(element_lengths == expected_lengths),
    msg = sprintf(
      "Lengths of tp_source elements must reduce by 1 sequentially. Found: %s. Expected: %s.",
      paste(element_lengths, collapse = ", "),
      paste(expected_lengths, collapse = ", ")
    )
  )

  # Name check: names of each state must be a subset of the previous state's
  # names. Compute all pairwise setdiffs at once, then report the first
  # violation.
  assertthat::assert_that(
    all(vapply(tp_source, function(x) !is.null(names(x)), logical(1L))),
    msg = "All inner lists in tp_source must be named"
  )

  missing_names <- Map(
    function(curr, nxt) setdiff(names(nxt), names(curr)),
    tp_source[-n_states],
    tp_source[-1L]
  )
  violation_idx <- which(lengths(missing_names) > 0L)
  if (length(violation_idx) > 0L) {
    i <- violation_idx[1L]
    stop(sprintf(
      "Topology error: State '%s' contains transitions to {%s} which are not present in the previous state '%s'.",
      names(tp_source)[i + 1L],
      paste(missing_names[[i]], collapse = ", "),
      names(tp_source)[i]
    ))
  }

  # ------------------------------------------------------------------
  # Determine expected vector lengths per state
  # ------------------------------------------------------------------

  if (is.null(tunnel_lengths)) {
    # Uniform case: all vectors must share a single common length. Infer it
    # from the first non-null, non-empty vector.
    all_vectors <- unlist(tp_source, recursive = FALSE)
    vec_lens <- lengths(all_vectors)
    n_cyc <- vec_lens[vec_lens > 0L][1L]
    assertthat::assert_that(
      !is.na(n_cyc),
      msg = "Could not determine number of cycles (all vectors are NULL or empty)"
    )
    # One expected length per state
    expected_vec_lengths <- rep(n_cyc, n_states)

  } else {
    # Variable-length case.
    n_tunnels <- length(tunnel_lengths)
    assertthat::assert_that(
      n_tunnels >= 1,
      n_tunnels <= n_states - 1,
      msg = sprintf(
        paste0(
          "length(tunnel_lengths) must be >= 1 and <= length(tp_source) - 1.",
          " Got %d tunnels for %d states (need at least 1 pre-tunnel state)."
        ),
        n_tunnels,
        n_states
      )
    )

    pre_tunnels <- n_states - n_tunnels

    # Infer th (pre-tunnel horizon) from first non-null, non-empty vector in
    # the pre-tunnel states.
    pre_tun_vectors <- unlist(tp_source[seq_len(pre_tunnels)], recursive = FALSE)
    vec_lens <- lengths(pre_tun_vectors)
    th <- vec_lens[vec_lens > 0L][1L]
    assertthat::assert_that(
      !is.na(th),
      msg = "Could not determine pre-tunnel horizon (all pre-tunnel vectors are NULL or empty)"
    )

    # Pre-tunnel states all share th; tunnel k has length tunnel_lengths[k]
    expected_vec_lengths <- c(rep(th, pre_tunnels), tunnel_lengths)
  }

  # ------------------------------------------------------------------
  # Sanitize: validate lengths/values and replace NULLs
  # ------------------------------------------------------------------

  sanitized_tp <- lapply(seq_along(tp_source), function(state_idx) {
    state_name <- names(tp_source)[state_idx]
    state_list <- tp_source[[state_idx]]
    exp_len    <- expected_vec_lengths[state_idx]

    new_state_list <- lapply(names(state_list), function(trans_name) {
      vec <- state_list[[trans_name]]

      if (is.null(vec)) {
        return(numeric(exp_len))
      }

      assertthat::assert_that(
        length(vec) == exp_len,
        msg = sprintf(
          "Length mismatch in transition %s -> %s. Expected %d, found %d.",
          state_name,
          trans_name,
          exp_len,
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

  # Summation constraint: element-wise sum of all transitions out of a state
  # must be <= 1 for every cycle. Compute all row-sums at once, then report
  # the first violating state.
  total_probs <- lapply(sanitized_tp, function(state_list) Reduce(`+`, state_list))
  violations <- which(vapply(total_probs, function(tp) any(tp > 1 + 1e-10), logical(1L)))
  if (length(violations) > 0L) {
    stop(sprintf(
      "Probabilities for state '%s' sum to > 1 in some cycles.",
      names(sanitized_tp)[violations[1L]]
    ))
  }

  sanitized_tp
}
