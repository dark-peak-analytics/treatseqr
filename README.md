# heseqmat

Health economics sequence modelling using sparse matrices

## Description

This package uses the `Matrix` package as well as standard `base` `R` to generate Markov traces for models of successive tunnel states.

Within a given health state, patients can do one of three things:

- `remain`: The patient can stay in this state (i.e., remain on treatment in the current line)
- `advance`: The patient can move to a (i.e., any) later tunnel-state
- `die`: The patient can die whilst in this state
