# Controlz.jl

`Controlz.jl` is a pure-[Julia](https://julialang.org/) package to analyze and simulate process dynamics and control systems using transfer function representations.

for example, to simulate the unit step response of a second-order, underdamped system characterized by the transfer function

$$g(s) = \dfrac{4}{4s^2 + 0.8s +1}$$

the output $Y(s)$ follows from $g(s)U(s)$, where $U(s)$ is the input.

```julia
using Controlz

g = 4 / (4 * s ^ 2 + 0.8 * s + 1) # construct transfer function
U = 1 / s                         # unit step input, U(s)
Y = g * U                         # system output, Y(s)

data = simulate(Y, 50.0)          # simulate until t = 50

viz_response(data, plot_title="SO underdamped step response")
```

![](SO_underdamped_step_response.png)

# install the `Controlz.jl` package in Julia

`Controlz.jl` is an officially registered Julia package. install in the Julia REPL by typing `]` to go into package mode, then `add Controlz`.

I recommend interactive [Pluto notebooks](https://github.com/fonsp/Pluto.jl) for coding in Julia, whose automatic package manager installs `Controlz.jl` upon running `using Controlz`.
