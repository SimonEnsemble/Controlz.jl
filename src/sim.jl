"""
    A, B, C = ss_to_tf(tf)

this formulation is based on the notes here: 
http://web.mit.edu/2.14/www/Handouts/StateSpace.pdf

         Y(s)    bₙ sⁿ + ... + b₁s + b₀
 G(s) = ----- = -------------------------
         U(s)    aₙ sⁿ + ... + a₁s + a₀

then we can formulate the system described by transfer
function G(s) as a state space representation known
as the controllable canonical form:

 dx
---- = A * x(t) + B * u(t)
 dt

y(t) = C * x(t) + bₙ/aₙ * u(t)

we compute the matrices A, B, C here.
"""
function tf_to_ss(tf::TransferFunction)
    @assert proper(tf) "tf not proper!"
    
    a_i(i::Int) = tf.denominator[i]
    b_i(i::Int) = tf.numerator[i]

    # get the highest degree in the num/den, which is that in the den if proper.
    n = degree(tf.denominator) # asserted proper above

    A = zeros(n, n)
    for i = 1:n-1
        A[i, i + 1] = 1.0
    end
    for i = 1:n
        A[n, i] = - a_i(i-1)
    end
    A[n, :] /= a_i(n)
    
    B = zeros(n, 1)
    B[n, 1] = 1 / a_i(n)
    
    C = zeros(1, n)
    for i = 1:n
        C[1, i] = b_i(i-1) - b_i(n) * a_i(i-1) / a_i(n)
    end
    return A, B, C
end


@doc raw"""
    data = simulate(Y, final_time, nb_time_points=100) # invert Y(s)

simulate the output $y(t)$ of an LTI system, given the Laplace transform of the output, $Y(s)$, `Y`.

in other words, `simulate` inverts an expression in the frequency domain into the time domain.

# arguments
* `Y::TransferFunction`: the Laplace transform of the output $y(t)$. usually formed by $g(s)U(s)$, where $U(s)$ is the Laplace transform of the input and $g(s)$ is the transfer function governing the dynamics of the system.
* `final_time::Tuple{Float64, Float64}`: the duration over which to simulate the output of the LTI system, starting at time zero.
* `nb_time_points::Int=100`: the number of time points at which to save the solution $y(t)$

two time points preceding $t=0$ are included to illustrate that it is assumed $y(t)=0$ for $t<0$.

# returns
* `data::DataFrame`: data frame containing two columns: `:t` for time $t$ and `:output` for $y(t)$. each row corresponds to a $(t_i, y(t_i))$ pair. i.e., row `i` of the `:t` column is time $i$, $t_i$, and row `i` of the `:output` column is $y_i=y(t_i)$. access the columns by `data[:, :t]` and `data[:, :output]`.

# example
simulate the first order step response to a step, given the Laplace transform of the output, `Y`:

```
g = 4 / (3 * s + 1)      # first-order transfer function g(s)
U = 1 / s                # unit step input U(s)
Y = g / s                # output Y(s)
data = simulate(Y, 12.0) # time series data frame
data[:, :t]              # array of time points tᵢ
data[:, :output]         # array of corresponding outputs y(tᵢ)
```
"""
function simulate(Y::TransferFunction, final_time::Float64; nb_time_points::Int=100)
    if ! proper(Y)
        error("LTI system is not proper...")
    end
    
    # coeffs on denominator, numerator polynomial
    a_i(i::Int) = Y.denominator[i]
    b_i(i::Int) = Y.numerator[i]
    n = degree(Y.denominator)

    # convert tf to state space form
    A, B, C = tf_to_ss(Y)

    x0 = deepcopy(B) # initial condition

    f(x, p, t) = A * x # RHS of ODE (ignore p for params)
    prob = ODEProblem(f, x0, (0.0, final_time))
    sol = solve(prob, d_discontinuities=[0.0])

    t = vcat([-0.05 * final_time, -1e-5], range(1e-5, final_time, length=nb_time_points - 2))
    y = [NaN for i = 1:nb_time_points]
    for (i, t_i) in enumerate(t)
        if t_i < Y.time_delay
            y[i] = 0.0
        else
            y[i] = (C * sol(t_i - Y.time_delay))[1] # put in shifted time
        end
    end
    return DataFrame(t=t, output=y)
end

@doc raw"""
    y_at_Τ = interpolate(data, Τ)

interpolate a data frame containing a time series characterizing $y(t)$, with `:t` and `:output` columns, whose rows are $(t_i, y(t_i))$ pairs.
interpolate the data to approximate the function $y(t)$ at a new time `Τ`, i.e. $y(\Tau)$.

[`simulate`](@ref) returns such a data frame, thus `interpolate` is useful for obtaining the solution at a particular time `Τ` that is not necessarily present in the `:t` column of `data`.

# arguments
* `data::DataFrame`: a data frame of a time series, containing a `:t` column for times and `:output` column for outputs.
* `Τ::Float64`: the new time at which we wish to know $y$. i.e. we wish to know $y(\Tau)$.

# Returns
* `y_at_Τ::Float64`: the value of $y$ when $t$ is equal to `Τ`, $y(\Tau)$, according to linear interpolation.

# example
the unit step response of a first-order process with time constant $\tau$ is $\approx 63\%$ of the final value when $t=\tau$.

```julia
τ = 3.45
g = 1 / (τ * s + 1)           # FO system
data = simulate(g / s, 10.0)  # unit step response
y_at_τ = interpolate(data, τ) # 0.63
```
"""
function interpolate(data::DataFrame, t̃::Float64)
    if t̃ < minimum(data[:, :t]) || t̃ > maximum(data[:, :t])
        error("t̃ = $t̃ ∉ time range of data, [$(minimum(data[:, :t])), $(maximum(data[:, :t]))].\n")
    end
    sort!(data, :t)
    # from Interpolations.jl.
    #   will throw error if t not sorted. 
    return LinearInterpolation(data[:, :t], data[:, :output])(t̃)
end
