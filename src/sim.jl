"""
    A, B, C = ss_to_tf(tf)

this formulation is based on the notes here: 
http://web.mit.edu/2.14/www/Handouts/StateSpace.pdf

         Y(s)    bₙ sⁿ + ... + b₁s + b₀
 G(s) = ----- = -------------------------
         U(s)    aₙ sⁿ + ... + a₁s + a₀

then we can formulate the tf G(s) as a state space ODE:

 dx
---- = A * x(t) + B * u(t)
 dt

y(t) = C * x(t) + bₙ/aₙ * u(t)

we compute A, B, C here.

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


"""
    t, y = simulate(tf, u, tspan, nb_time_points=100)

Simulate LTI system described by transfer function `tf` with input function `u` (a function of time).

# Arguments
* `tf::TransferFunction`: the transfer function describing the relationship between input `u` and output `y`
* `u::Function`: the input function u(t)
* `tspan::Tuple{Float64, Float64}`: the time span over which to simulate the LTI system.
* `nb_time_points::Int=100`: the number of time points at which to save the solution y(t)

# Returns
* `t::Array{Float64, 1}`: array of times at which solution was saved
* `y::Array{Float64, 1}`: array of y values at corresponding times in `t`

# Example

julia> tf = 4 / (3 * s + 1)

one can use anonymous functions as the input:

julia> t, y = simulate(tf, t -> (t < 0.0) ? 0.0 : 1.0, (0.0, 12.0)) # unit step response

or explicit functions:
julia> a_step_input(t) = (t < 0.0) ? 0.0 : 1.0
julia> t, y = simulate(tf, a_step_input, (0.0, 12.0)) # unit step response
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
        y[i] = (C * sol(t_i))[1] + b_i(n) / a_i(n) * u(t_i - tf.time_delay)
    end
    return t, y
end
