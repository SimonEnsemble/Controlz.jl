# Simulation

we wish to simulate the response $y(t)$ (output) of a linear, time-invariant system, characterized by a transfer function $g(s)$, to an input $u(s)$. 

![](assets/simple_input_output.png)

pass the output $Y(s)$ in the frequency domain into the function `simulate` to invert it into the time domain to obtain $y(t)$. 
`simulate` returns a time series data frame (a `DataFrame`, see `DataFrames.jl` [docs](https://dataframes.juliadata.org/latest/)), with a `:t` column for the times and an `:output` column for $y(t)$. we can then plot the time series data, interpolate the data to obtain the value of $y(t)$ at a particular time $t=\Tau$, etc. we provide examples below.

## response of an underdamped second-order system to a unit step input

```julia
g = 4 / (4 * s ^ 2 + 0.8 * s + 1) # second order transfer function, underdamped

U = 1 / s # unit step input
Y = g * U # system output

data = simulate(Y, 50.0) # simulate until t = 50, returns DataFrame
data[:, :t]      # array of times, tᵢ's
data[:, :output] # array of outputs, yᵢ's ≈ y(tᵢ)'s
```

we can then plot the time series via:

```julia
viz_response(data, plot_title="SO underdamped step response")
```

![](SO_underdamped_step_response.png)

## response of a first-order plus time delay system to a unit step input

```julia
K = 2.0 # gain
τ = 4.0 # time constant
θ = 1.5 # time delay
g = K * exp(-θ * s) / (τ * s + 1) # FOPTD transfer function

U = 1 / s # step input
Y = g * U

data = simulate(Y, 15.0) # simulate until t = 15

viz_response(data, plot_title="FOPTD step response")
```

![](FOPTD_step_response.png)

## inverse Laplace transform

to emphasize that our `simulate` function takes a function of the complex frequency `s` and inverts it into the time domain, consider the Laplace transform of $t \cos(at)$:

$$\mathcal{L}[t \cos(at)] = \dfrac{s^2-a^2}{(s^2+a^2)^2}.$$

we can numerically invert $\dfrac{s^2-a^2}{(s^2+a^2)^2}$ as follows:

```julia
a = π
U = (s^2 - a^2) / (s^2 + a^2) ^ 2

data = simulate(U, 8.0, nb_time_points=300) # simulate until t=8, use 300 time points for high resolution

viz_response(data, plot_title=L"inverting an input $U(s)$", plot_ylabel=L"$u(t)$")
```

the `nb_time_points` argument allows us to return a time series with a higher resolution in time. if the plot of the response appears jagged, likely you need to increase `nb_time_points`.

![](tcosat.png)

## $y(t)$ at an arbitrary time $\Tau$

the [`simulate`](@ref) function returns an array of times $t_i$'s and corresponding $y_i=y(t_i)$'s. if we wish to know $y(t)$ at a particular time $t=\Tau$, we can call `interpolate` to linearly interpolate the time series data.

for example, to obtain the output of a first-order system with time constant $\tau$ in response to a unit step input at $t=\tau$:

```julia
τ = 3.45
g = 1 / (τ * s + 1) # FO system
data = simulate(g / s, 10.0)  # unit step response
y_at_τ = interpolate(data, τ) # 0.63 ≈ 1 - 1/e, as we expect
```

## under the hood

under the hood, `simulate` converts the system passed to it into a state space ODE (a system of ODEs) in the time domain and uses `DifferentialEquations.jl` (see [here](https://github.com/JuliaDiffEq/DifferentialEquations.jl)) to numerically solve the resulting ODE.

# detailed docs

```@docs
    simulate
    TransferFunction
    interpolate
```
