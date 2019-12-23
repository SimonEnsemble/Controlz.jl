# Visualization

## poles and zeros of a transfer function

```
julia> g = (s + 2) / (s^2 + 1/4)
julia> viz_poles_and_zeros(g)
```

![](example_poles_and_zeros.png)

## response of a system to an input

```
julia> g = 4 / (4 * s ^ 2 + 0.8 * s + 1)
julia> u = 1 / s
julia> t, y = simulate(g * u, (0.0, 50.0))
viz_response(t, y, plot_title="SO underdamped step response")
```

![](example_response.png)

## Nyquist diagram

```
julia> g = 1 / (s^2 + s + 1)
julia> nyquist_diagram(g)
```

![](example_nyquist.png)

## Bode plot

```
julia> g = 3 / (s + 1)
julia> bode_plot(g, log10_ω_min=-4.0, log10_ω_max=4.0)
```

![](example_bode.png)

## Root locus plot

```
julia> g_ol = 4 / (s + 3) / (s + 2) / (s + 1)
julia> root_locus(g_ol)
```

![](example_root_locus.png)

## detailed docs

```@docs
    viz_poles_and_zeros
    viz_response
    nyquist_diagram
    bode_plot
    root_locus
```
