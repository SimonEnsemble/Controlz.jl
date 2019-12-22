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
    t, y = simulate(tf, u, tspan, nb_time_points=100) # explicit input function u(t)
    t, y = simulate(Y, tspan, nb_time_points=100) # invert Y(s)

Simulate the output $y(t)$ of an LTI system. `simulate` handles two scenarios:
1. we have the transfer function `tf` that characterizes how the LTI system responds to inputs and an input function `u`, an explicit function of time $u(t)$.
2. we have the Laplace transform of the output, $Y(s)$, `Y`.

# Arguments
* `tf::TransferFunction`: the transfer function describing the relationship between input `u` and output `y`
* `Y::TransferFunction`: the Laplace transform of the output $y(t)$. Usually formed by $g(s)U(s)$, where $U(s)$ is the Laplace transform of the input.
* `u::Function`: the input function u(t)
* `tspan::Tuple{Float64, Float64}`: the time span over which to simulate the LTI system.
* `nb_time_points::Int=100`: the number of time points at which to save the solution $y(t)$

# Returns
* `t::Array{Float64, 1}`: array of times $t$ at which the solution was saved
* `y::Array{Float64, 1}`: array of $y$ values at corresponding times in `t`

# Example

Ex. 1: given transfer function `tf` and input function `u`

One can simulate the first order step response as:
```
julia> tf = 4 / (3 * s + 1)
julia> u(t) = (t < 0.0) ? 0.0 : 1.0
julia> t, y = simulate(tf, u, (0.0, 12.0))
```

Ex. 2: given Laplace transform of the output, `Y`

One can also simulate the first order step response as:
```
julia> tf = 4 / (3 * s + 1)
julia> Y = tf / s
julia> t, y = simulate(Y, (0.0, 12.0))
```
"""
function simulate(tf::TransferFunction, u::Function, tspan::Tuple{Float64, Float64};
        nb_time_points::Int=100)
    if ! proper(tf)
        error("transfer function not proper...")
    end

    a_i(i::Int) = tf.denominator[i]
    b_i(i::Int) = tf.numerator[i]
    n = degree(tf.denominator)

    # convert tf to state space form
    A, B, C = tf_to_ss(tf)

    x0 = zeros(n, 1) # initial condition

    f(x, p, t) = A * x + B * u(t - tf.time_delay) # RHS of ODE (ignore p for params)
    prob = ODEProblem(f, x0, tspan)
    # pass time delay as discontinuity to avoid spurious shift in result
    sol = solve(prob, d_discontinuities=[tf.time_delay])

    t = range(tspan[1], stop=tspan[2], length=nb_time_points)
    y = [NaN for i = 1:nb_time_points]
    for (i, t_i) in enumerate(t)
        if t_i < 0.0
            y[i] = 0.0
        else
            y[i] = (C * sol(t_i))[1] + b_i(n) / a_i(n) * u(t_i - tf.time_delay)
        end
    end
    return collect(t), y
end

function simulate(Y::TransferFunction, tspan::Tuple{Float64, Float64}; nb_time_points::Int=100)
    if ! proper(Y)
        error("transfer function not proper...")
    end

    a_i(i::Int) = Y.denominator[i]
    b_i(i::Int) = Y.numerator[i]
    n = degree(Y.denominator)

    # convert tf to state space form
    A, B, C = tf_to_ss(Y)

    x0 = deepcopy(B) # initial condition

    f(x, p, t) = A * x # RHS of ODE (ignore p for params)
    prob = ODEProblem(f, x0, tspan)
    # pass time delay as discontinuity to avoid spurious shift in result
    sol = solve(prob, d_discontinuities=[Y.time_delay])

    t = range(tspan[1], stop=tspan[2], length=nb_time_points)
    y = [NaN for i = 1:nb_time_points]
    for (i, t_i) in enumerate(t)
        if t_i < 0.0
            y[i] = 0.0
        else
            y[i] = (C * sol(t_i))[1]
        end
    end
    return collect(t), y
end
