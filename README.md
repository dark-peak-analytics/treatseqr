# treatseqr

Health economics sequence modelling using sparse matrices

## Description

This package uses the `Matrix` package as well as standard `base` `R` to generate Markov traces for models of successive tunnel states.

Within a given health state, patients can do one of four things:

- `stay`: The patient can stay in this state (i.e., stay on treatment in the current line)
- `advance`: The patient can move to a (i.e., any) later tunnel-state
- `return`: The patient can move back to the start of a (i.e., any) earlier tunnel-state, for instance where a relapse resolves without accumulating permanent disability
- `die`: The patient can die whilst in this state

The complicating factor is that for each state past the first (e.g., first-line on treatment), the probability `advance` does not depend on time from baseline, `t`, but time within the _current_ state, `d` (for difference between current `t` and `t` at the time of entry to the current state). This means that a time-constant transition probability matrix (`tpm`) cannot be used to accurately extrapolate state residency, unless all event hazards are constant (as is the case for exponential hazards).

One strategy to arrive at an exact solution even in this context is to use "tunnel states". These are states which explicitly track `d` for each state, spreading out the cohort into sub-cohorts by arrival time into the tunnel. In this approach, each treatment line is represented by a series of tunnel states corresponding to a specific situation (e.g., on active treatment in second-line for `d = 0, 1, ..., th`, where `th` is the time horizon).

The simplest tunnel in an oncology setting might be corresponding to Kaplan-Meier (`KM`) data of best supportive care (`BSC`) from the point of arriving there. This cannot be entered into a "typical" Markov modelling framework as the basis for the `KM` data is `d` (time in `BSC`), not `t` (time from baseline in the model). In this case, the only destination for the cohort(s) to `advance` to is `dead`, as there are no subsequent treatment lines.

This complexity is compounded when thinking about _multiple_ tunnels, and the idea that a given tunnel-state can have _multiple_ destinations. For instance, if there are 5 different treatment lines, and patients can "skip" lines (e.g., move from first-line directly to third-line), then the number of possible transitions increases rapidly.

However, constructing such a framework of transition probabilities and using them to produce a Markov trace is possible, using some relatively simple matrix algebra. This package implements this approach, using sparse matrices to keep memory usage manageable. The Matrix approach has the advantage over patient-level simulation because the solution is exact (up to numerical precision), and runs very quickly even for large models in comparison to e.g., a simulated cohort of 5000 patients.

### Specification

There several features to the specification of a model following this framework:

- The top row (or rows — up to 10 pre-tunnel states are supported) of the `tp` changes every model cycle, to reflect the fact that the `advance` probabilities depend on `t` (time from baseline). This reduces the size of the final matrix by `th-1` rows and columns, where `th` is the time horizon. This is referred to as `m1` in the code and documentation.
- The rest of the matrix does **not** change over time, as the `advance` probabilities depend on `d` (time in current state), which is tracked by the tunnel states and is compiled into the matrix itself. This is referred to as `m2` in the code and documentation.
- `matrix_size = pre_tunnel_states + sum(tunnel_lengths) + 1`, where `tunnel_lengths` gives the number of cycles in each tunnel (they need not be equal) and the trailing `1` is the absorbing death state.
- `m1` is a dense 3d array of dimension `c(th, n_pre, n_dest)`, indexed `[cycle, pre-tunnel state, destination]`, and is consequently relatively small even for large `th` — `n_dest` covers only the states a pre-tunnel state can reach, not all `matrix_size` columns.
- `m2` has dimensions `matrix_size` by `matrix_size`, and has an **empty top row**
- Within `m2` the transition probabilities are arranged such that each tunnel-state only has non-zero probabilities to `stay` in the same tunnel-state, move to the start of any other tunnel-state, or `die`, rendering `m2` mostly "empty" (i.e., sparse).
- In `m2`, `stay` probabilities are located in a superdiagonal (1) position (i.e., one to the right of the main diagonal for that tunnel-state), whilst transitions to other tunnels are arranged **vertically** at the column of the destination tunnel's first cell — to the right for `advance`, to the left for `return` — and `die` probabilities are located in the last column. Every tunnel block gets a vertical strip for every other tunnel plus death; unused ones stay at zero and are filtered out before the sparse matrix is built, so they cost essentially nothing. Visually, this looks like one long staircase going from top left to bottom right, with a series of vertical lines:

<img width="653" height="620" alt="image" src="https://github.com/user-attachments/assets/e32b19c8-e425-4a6a-9908-a03b1df4f1c1" />

Note here that the model is for two active treatment lines, best supportive care, and death, with the ability to "skip" treatment lines or go straight to BSC from any position in the pathway. The vertical lines are transitions to subsequent tunnel-states, whilst the diagonal lines (which are actually superdiagonal (1)), are `stay` probabilities.

The functions in this package simply compile `m1` (which is the missing top row of the above image) and `m2` (the image), and then run a Markov engine column-wise using sparse matrix multiplication to generate the Markov trace.

The Markov engine itself is remarkably simple. Here is an example:

```r
pop <- matrix(0, nrow = spec$matrix_size, ncol = th + 1)
pop[1, 1] <- 1
m2t <- Matrix::t(m$m2)
pre_pos <- seq_len(spec$pre_tunnels)

for (cyc in seq_len(th) + 1) {
  pop[, cyc] <- as.numeric(m2t %*% pop[, cyc - 1, drop = FALSE])
  pop[m$m1_dest, cyc] <- pop[m$m1_dest, cyc] +
    as.numeric(pop[pre_pos, cyc - 1] %*% m$m1[cyc - 1, , ])
}
```

Put simply, a population matrix `pop` is made with 100% of the cohort in the first state at baseline. Note that this matrix has time as columns, and states as rows. This is to optimise the performance of `R`, which is "column-major" or in other words substantially faster at accessing and manipulating columns than rows. Each cycle, the whole population vector in `t-1` is multiplied by `m2`, and then the pre-tunnel exits are computed from the `t-1` slice of the `m1` array and added into their destination columns — a "scatter add", which avoids doing a full matrix multiplication for the handful of rows that change every cycle. `m1_dest` records which columns of `M` those destinations are.

This looks very simple, but in reality can produce any number of tunnel states with probabilities of skipping states and within-state tracking to any level of complication. This enables and simplifies model structures which are difficult to implement outside of patient-level simulation.
