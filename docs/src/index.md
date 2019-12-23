# Controlz.jl

`Controlz.jl` is a [Julia](https://julialang.org/) package to explore concepts in the simulation of process dynamics and control of linear, time-invariant (LTI) systems using transfer function representations.

For example, to simulate the unit step response of a second-order, undamped system characterized by the transfer function:

$$g(s) = \dfrac{4}{4s^2 + 0.8s +1}$$

the output $Y(s)$ follows from $g(s)u(s)$.

```julia
julia> using Controlz
julia> g = 4 / (4 * s ^ 2 + 0.8 * s + 1) # construct transfer function
julia> u = 1 / s # unit step input, U(s)
julia> Y = g * u # system output, Y(s)
julia> t, y = simulate(Y, (0.0, 50.0)) # simulate from t = 0 to t = 50
julia> viz_response(t, y, title="SO underdamped step response")
```

![](example_response.png)

# install in Julia

* in the Julia REPL: go into package mode by typing `]`. Then `add Controlz#master`. Then `Backspace` to exit package mode.
* in Jupyter Notebook or Julia: `using Pkg; Pkg.add("Controlz#master")`.

