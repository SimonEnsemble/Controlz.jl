# Simulation

we wish simulate to simulate the response (output) of a linear, time-invariant system, characterized by a transfer function $g(s)$, to an input. we run the simulation in the familiar time domain by converting the system into state space and using `DifferentialEquations.jl` (see [here](https://github.com/JuliaDiffEq/DifferentialEquations.jl)) to solve the resulting ODE.

learn by example! in each case, `simulate` returns an array of times `t` and corresponding output values, `y`.

## response of an underdamped second-order system to a unit step input

```
g = 4 / (4 * s ^ 2 + 0.8 * s + 1) # construct transfer function
u = 1 / s # unit step input
Y = g * u # system output
t, y = simulate(Y, (0.0, 50.0)) # simulate from t = 0 to t = 50
```

we can then plot the `y` array versus the `t` array:

![](example_response.png)

```@docs
    simulate
```
