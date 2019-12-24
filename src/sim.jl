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
    t, y = simulate(Y, final_time, nb_time_points=100) # invert Y(s)

Simulate the output $y(t)$ of an LTI system, given the Laplace transform of the output, $Y(s)$, `Y`.

# Arguments
* `Y::TransferFunction`: the Laplace transform of the output $y(t)$. Usually formed by $g(s)U(s)$, where $U(s)$ is the Laplace transform of the input.
* `final_time::Tuple{Float64, Float64}`: the duration over which to simulate the output of the LTI system, starting at time zero.
* `nb_time_points::Int=100`: the number of time points at which to save the solution $y(t)$

Two points before time zero are included to illustrate that it is assumed $y(t)=0$ for $t<0$.

# Returns
* `t::Array{Float64, 1}`: array of times $t$ at which the solution was saved
* `y::Array{Float64, 1}`: array of $y$ values at corresponding times in `t`

# Example

One can simulate the first order step response as, given the Laplace transform of the output, `Y`:

```
julia> g = 4 / (3 * s + 1) # first-order transfer function
julia> u = 1 / s #  unit step input
julia> Y = g / s
julia> t, y = simulate(Y, 12.0)
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

    t = vcat([-0.05 * final_time, 0.0], range(1e-4, final_time, length=nb_time_points - 2))
    y = [NaN for i = 1:nb_time_points]
    for (i, t_i) in enumerate(t)
        if t_i < Y.time_delay
            y[i] = 0.0
        else
            y[i] = (C * sol(t_i - Y.time_delay))[1] # put in shifted time
        end
    end
    return t, y
end
