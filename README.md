# heseqmat

Health economics sequence modelling using sparse matrices

## Description

This package uses the `Matrix` package as well as standard `base` `R` to generate Markov traces for models of successive tunnel states.

Within a given health state, patients can do one of three things:

- `remain`: The patient can stay in this state (i.e., remain on treatment in the current line)
- `advance`: The patient can move to a (i.e., any) later tunnel-state
- `die`: The patient can die whilst in this state

The complicating factor is that for each state past the first (e.g., first-line on treatment), the probability `advance` does not depend on time from baseline, `t`, but time within the _current_ state, `d` (for difference between current `t` and `t` at the time of entry to the current state). This means that a time-constant transition probability matrix (`tpm`) cannot be used to accurately extrapolate state residency, unless all event hazards are constant (as is the case for exponential hazards).

One strategy to arrive at an exact solution even in this context is to use "tunnel states". These are states which explicitly track `d` for each state, spreading out the cohort into sub-cohorts by arrival time into the tunnel. In this approach, each treatment line is represented by a series of tunnel states corresponding to a specific situation (e.g., on active treatment in second-line for `d = 0, 1, ..., th`, where `th` is the time horizon).

The simplest tunnel in an oncology setting might be corresponding to Kaplan-Meier (`KM`) data of best supportive care (`BSC`) from the point of arriving there. This cannot be entered into a "typical" Markov modelling framework as the basis for the `KM` data is `d` (time in `BSC`), not `t` (time from baseline in the model). In this case, the only destination for the cohort(s) to `advance` to is `dead`, as there are no subsequent treatment lines.

This complexity is compounded when thinking about _multiple_ tunnels, and the idea that a given tunnel-state can have _multiple_ destinations. For instance, if there are 5 different treatment lines, and patients can "skip" lines (e.g., move from first-line directly to third-line), then the number of possible transitions increases rapidly.

However, constructing such a framework of transition probabilities and using them to produce a Markov trace is possible, using some relatively simple matrix algebra. This package implements this approach, using sparse matrices to keep memory usage manageable. The Matrix approach has the advantage over patient-level simulation because the solution is exact (up to numerical precision), and runs very quickly even for large models in comparison to e.g., a simulated cohort of 5000 patients.

### Specification

There several features to the specification of a model following this framework:

- The top row of the `tp` changes every model cycle, to reflect the fact that the `advance` probabilities depend on `t` (time from baseline). This reduces the size of the final matrix by `th-1` rows and columns, where `th` is the time horizon. This is referred to as `m1` in the code and documentation.
- The rest of the matrix does **not** change over time, as the `advance` probabilities depend on `d` (time in current state), which is tracked by the tunnel states and is compiled into the matrix itself. This is referred to as `m2` in the code and documentation.
- `m1` has 1 row and `matrix_size = ( th + 2 ) * n_tunnels` columns, and is consequently relatively small even for large `th`.
- `m2` has dimensions `matrix_size` by `matrix_size`, and has an **empty top row**
- Within `m2` the transition probabilities are arranged such that each tunnel-state only has non-zero probabilities to `remain` in the same tunnel-state, `advance` to later tunnel-states, or `die`, rendering `m2` mostly "empty" (i.e., sparse).
- In `m2`, `remain` probabilities are located in a superdiagonal (1) position (i.e., one to the right of the main diagonal for that tunnel-state), whilst `advance` probabilities are arranged **vertically** in positions further to the right (corresponding to the `+1th` cell of the previous tunnel block), and `die` probabilities are located in the last column. Visually, this looks like one long staircase going from top left to bottom right, with a series of vertical lines:
