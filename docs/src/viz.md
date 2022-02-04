# Visualization

the visualizations rely on [`CairoMakie.jl`](https://github.com/JuliaPlots/Makie.jl).

## poles and zeros of a transfer function

```julia
g = (s + 2) / (s^2 + 1/4)
viz_poles_and_zeros(g)
```

![](example_poles_and_zeros.png)

## response of a system to an input

```julia
g = 4 / (4 * s ^ 2 + 0.8 * s + 1)
U = 1 / s
Y = g * U
data = simulate(Y, (0.0, 50.0))
viz_response(data, title="SO underdamped step response")
```

![](SO_underdamped_step_response.png)

## Nyquist diagram

```julia
g = 1 / (s^2 + s + 1)
nyquist_diagram(g)
```

![](example_nyquist.png)

## Bode plot

```julia
g = 3 / (s + 1)
bode_plot(g, log10_ω_min=-4.0, log10_ω_max=4.0, nb_pts=300)
```

![](example_bode.png)

the range of frequencies presented is determined by `log10_ω_min` and `log10_ω_max`. the resolution of the Bode plot is determined by `nb_pts`.

see [`gain_phase_margins`](@ref) to compute the gain and phase margins and the critical and gain crossover frequencies.

## Root locus plot

```julia
g_ol = 4 / (s + 3) / (s + 2) / (s + 1)
root_locus(g_ol)
```

![](example_root_locus.png)

## modifying the figures

all visualization functions return a `Figure` object from `CairoMakie.jl` that can be further modified. for example:

```julia
g_ol = 4 / (s + 3) / (s + 2) / (s + 1)
fig = root_locus(g_ol)
ax = current_axis(fig)
xlims!(ax, -15, 5)
ax.xlabel = "real numbers"
```

## cool plot theme

the custom plot theme can be invoked in `CairoMakie.jl` for other plots via:
```julia
using Controlz, CairoMakie
set_theme!(cool_theme)
```
more, `CairoMakie.jl` offers other themes [here](https://makie.juliaplots.org/dev/documentation/theming/predefined_themes/).

## detailed docs

```@docs
    viz_response
    viz_poles_and_zeros
    nyquist_diagram
    bode_plot
    root_locus
    mk_gif
```
