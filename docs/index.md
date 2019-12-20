# Controlz.jl

`Controlz.jl` is a [Julia](https://julialang.org/) package to explore concepts in the simulation of process dynamics and control of linear, time-invariant (LTI) systems using transfer function representations.

For example, consider an LTI system characterized by a transfer function $g(s) = 4 / (3s+1)$.

```julia
julia> using Controlz

julia> g = 4 / (9 * s ^ 2 + s + 1) # transfer function

julia> z, p, k = zeros_poles_gain(g) # compute zeros, poles, and gain

julia> u = 1 / s # step input

julia> t, y = simulate(g * u, (0.0, 40.0)) # simulate from t = 0 to t = 40

julia> viz_response(t, y, title="SO underdamped step response")
```

![](example_response.png)
