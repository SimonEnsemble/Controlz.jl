# Simulation

we wish simulate to simulate the response (output) of a linear, time-invariant system, characterized by a transfer function $g(s)$, to an input. 

we will conduct this simulation in the familiar time domain.

learn by example! in each case, `simulate` returns an array of times `t` and corresponding output values, `y`.

## response of an underdamped second-order system to a unit step input

```
julia> g = 4 / (4 * s ^ 2 + 0.8 * s + 1) # construct transfer function
julia> u = 1 / s # unit step input
julia> Y = g * u # system output
julia> t, y = simulate(Y, (0.0, 50.0)) # simulate from t = 0 to t = 40
```

```@docs
    simulate
```
