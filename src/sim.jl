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

# Returns
* `t::Array{Float64, 1}`
* `y::Array{Float64, 1}`
"""
function simulate(tf::TransferFunction, u::Function, tspan::Tuple{Float64, Float64};
        nb_time_points::Int=100)
    if ! proper(tf)
        error("transfer function not proper...")
    end
    # TODO implement time delays
    if tf.time_delay != 0.0
        error("time delay not implemented yet")
    end

    a_i(i::Int) = tf.denominator[i]
    b_i(i::Int) = tf.numerator[i]
    n = degree(tf.denominator)

    # convert tf to state space form
    A, B, C = tf_to_ss(tf)

    x0 = zeros(n, 1) # initial condition

    f(x, p, t) = A * x + B * u(t) # RHS of ODE (ignore p for params)
    prob = ODEProblem(f, x0, tspan)
    sol = solve(prob)

    t = range(tspan[1], stop=tspan[2], length=nb_time_points)
    y = [NaN for i = 1:nb_time_points]
    for (i, t_i) in enumerate(t)
        y[i] = (C * sol(t_i))[1] + b_i(n) / a_i(n) * u(t_i)
    end
    return t, y
end
