# Visualization

## poles and zeros of a transfer function

```
g = (s + 2) / (s^2 + 1/4)
viz_poles_and_zeros(g)
```

![](example_poles_and_zeros.png)

## response of a system to an input

```
g = 4 / (4 * s ^ 2 + 0.8 * s + 1)
U = 1 / s
Y = g * U
data = simulate(Y, (0.0, 50.0))
viz_response(data, plot_title="SO underdamped step response")
```

![](SO_underdamped_step_response.png)

## Nyquist diagram

```
g = 1 / (s^2 + s + 1)
nyquist_diagram(g)
```

![](example_nyquist.png)

## Bode plot

```
g = 3 / (s + 1)
bode_plot(g, log10_ω_min=-4.0, log10_ω_max=4.0, nb_pts=300)
```

![](example_bode.png)

the range of frequencies presented is determined by `log10_ω_min` and `log10_ω_max`. the resolution of the Bode plot is determined by `nb_pts`.

see [`gain_phase_margins`](@ref) to compute the gain and phase margins and the critical and gain crossover frequencies.

## Root locus plot

```
g_ol = 4 / (s + 3) / (s + 2) / (s + 1)
root_locus(g_ol)
```

![](example_root_locus.png)

## hipster plot theme

invoke the hipster plot theme used to make plots for this documentation by:

```
using PyPlot
PyPlot.matplotlib.style.use("https://raw.githubusercontent.com/SimonEnsemble/Controlz.jl/master/src/hipster.mplstyle")
```

## detailed docs

```@docs
    viz_poles_and_zeros
    viz_response
    nyquist_diagram
    bode_plot
    root_locus
    mk_gif
```
