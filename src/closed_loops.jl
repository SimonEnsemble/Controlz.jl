import Base.*, Base./

"""
a closed-loop transfer function that relates an output `Y` and an input `U` in a feedback loop.

the resulting closed-loop transfer function is:

```
 Y      top
--- = --------
 U    1 + g_ol
```

# example

```jldoctest
g_ol = 4 / (s + 1) * 2 / (s + 2)
top = 5 / (s + 4)
g = ClosedLoopTransferFunction(top, g_ol)
# output
closed-loop transfer function.
      top
    -------
    1 + g_ol

  top =
    5.0
-----------
1.0*s + 4.0

  g_ol =
         8.0
---------------------
1.0*s^2 + 3.0*s + 2.0
```                                             

# attributes
* `top::TransferFunction`: numerator
* `g_ol::TransferFunction`: open-loop transfer function
"""
struct ClosedLoopTransferFunction
	top::TransferFunction
	g_ol::TransferFunction
end

# algebra
*(tf_cl::ClosedLoopTransferFunction, tf::TransferFunction) = ClosedLoopTransferFunction(tf_cl.top * tf, tf_cl.g_ol)
*(tf::TransferFunction, tf_cl::ClosedLoopTransferFunction) = *(tf_cl::ClosedLoopTransferFunction, tf::TransferFunction)
*(tf_cl::ClosedLoopTransferFunction, x::Number) = ClosedLoopTransferFunction(tf_cl.top * x, tf_cl.g_ol)
*(x::Number, tf_cl::ClosedLoopTransferFunction) = *(tf_cl::ClosedLoopTransferFunction, x::Number)
/(tf_cl::ClosedLoopTransferFunction, tf::TransferFunction) = ClosedLoopTransferFunction(tf_cl.top / tf, tf_cl.g_ol)
/(tf_cl::ClosedLoopTransferFunction, x::Number) = ClosedLoopTransferFunction(tf_cl.top / x, tf_cl.g_ol)

# representation of a closed loop system in standard form.
# useful for converting into state space representation.
#        p_a(s) e^ {-θs}
#  ---------------------------
#   p_b(s) + p_c(s) e^ {-ϕs}
struct CLTFStandard
    # numerator
    p_a::Polynomial
    θ::Union{Float64, Int}
    # denominator
    p_b::Polynomial
    p_c::Polynomial
    ϕ::Union{Float64, Int}
end

function CLTFStandard(cl::ClosedLoopTransferFunction)
    # numerator
    p_a = cl.top.numerator * cl.g_ol.denominator
    θ = cl.top.time_delay
    # denominator
    p_b = cl.g_ol.denominator * cl.top.denominator
    p_c = cl.g_ol.numerator   * cl.top.denominator
    ϕ = cl.g_ol.time_delay
    return CLTFStandard(p_a, θ, p_b, p_c, ϕ)
end

order(cl::CLTFStandard) = degree(cl.p_b)
order(cl::ClosedLoopTransferFunction) = degree(CLTFStandard(cl).p_b)

function strictly_proper(cl::CLTFStandard)
    return degree(cl.p_b) > max(degree(cl.p_a), degree(cl.p_c))
end
strictly_proper(cl::ClosedLoopTransferFunction) = strictly_proper(CLTFStandard(cl))

"""
 dx
---- = A * x(t) + B * u(t) + C * x(t - ϕ)
 dt
y(t) = D * x(t)
we compute the matrices A, B, C, D here.
only works if strictly proper.
"""
function tf_to_ss(cl::CLTFStandard)
    @assert strictly_proper(cl) "closed-loop system not strictly proper!"

    a_i(i::Int) = cl.p_a[i]
    b_i(i::Int) = cl.p_b[i]
	c_i(i::Int) = cl.p_c[i]

    n = order(cl)

    A = zeros(n, n)
    for i = 1:n-1
        A[i, i + 1] = 1.0
    end
    for i = 1:n
        A[n, i] = - b_i(i-1) / b_i(n)
    end

    B = zeros(n, 1)
    B[n, 1] = 1 / b_i(n)

    C = zeros(n, n)
    for i = 1:n
        C[n, i] = c_i(i-1) / b_i(n)
    end

	D = zeros(1, n)
	for i = 1:n
		D[i] = a_i(i-1)
	end

    return A, B, C, D
end

# see sim.jl for doc string
function simulate(cl::ClosedLoopTransferFunction, final_time::Float64; nb_time_points::Int=100)
	cls = CLTFStandard(cl)
	
    if ! strictly_proper(cls)
        error("closed loop system is not strictly proper...")
    end

    # convert tf to state space form
    A, B, C, D = tf_to_ss(cls)

    x0 = deepcopy(B) # initial condition
	# h is history
	h(p, t) = zeros(order(cls))
    f(x, h, p, t) = A * x - C * h(p, t - cls.ϕ) # RHS of ODE (ignore p for params)
    prob = DDEProblem(f, x0, h, (0.0, final_time))
    sol = solve(prob, d_discontinuities=[0.0, cls.ϕ, cls.θ, cls.θ + cls.ϕ])

	t = vcat([-0.05 * final_time, -1e-5], 
	range(1e-5, final_time, length=nb_time_points - 2))
	y = [NaN for i = 1:nb_time_points]
	for (i, t_i) in enumerate(t)
		if t_i < cls.θ
			y[i] = 0.0
		else
			y[i] = (D * sol(t_i - cls.θ))[1] # put in shifted time
		end
	end
    return DataFrame(t=t, output=y)
end
