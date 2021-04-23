import Base.*, Base./

"""
closed-loop system with output Y and input D related as:

 Y      top
--- = --------
 D    1 + g_ol

* `top::TransferFunction`
* `g_ol::TransferFunction`
"""
struct ClosedLoopTF
	top::TransferFunction
	g_ol::TransferFunction
end

# algebra
*(tf_cl::ClosedLoopTF, tf::TransferFunction) = ClosedLoopTF(tf_cl.top * tf, tf_cl.g_ol)
*(tf::TransferFunction, tf_cl::ClosedLoopTF) = *(tf_cl::ClosedLoopTF, tf::TransferFunction)
*(tf_cl::ClosedLoopTF, x::Number) = ClosedLoopTF(tf_cl.top * x, tf_cl.g_ol)
*(x::Number, tf_cl::ClosedLoopTF) = *(tf_cl::ClosedLoopTF, x::Number)
/(tf_cl::ClosedLoopTF, tf::TransferFunction) = ClosedLoopTF(tf_cl.top / tf, tf_cl.g_ol)
/(tf_cl::ClosedLoopTF, x::Number) = ClosedLoopTF(tf_cl.top / x, tf_cl.g_ol)

# representation of a closed loop system in standard form.
#        p_a(s) e^ {-θs}
#  ---------------------------
#   p_b(s) + p_c(s) e^ {-ϕs}
struct ClosedLoopTFStandard
    # numerator
    p_a::Poly
    θ::Union{Float64, Int}
    # denominator
    p_b::Poly
    p_c::Poly
    ϕ::Union{Float64, Int}
end

function ClosedLoopTFStandard(cl::ClosedLoopTF)
    # numerator
    p_a = cl.top.numerator * cl.g_ol.denominator
    θ = cl.top.time_delay
    # denominator
    p_b = cl.g_ol.denominator * cl.top.denominator
    p_c = cl.g_ol.numerator   * cl.top.denominator
    ϕ = cl.g_ol.time_delay
    return ClosedLoopTFStandard(p_a, θ, p_b, p_c, ϕ)
end

order(cl::ClosedLoopTFStandard) = degree(cl.p_b)
order(cl::ClosedLoopTF) = degree(ClosedLoopTFStandard(cl).p_b)

function strictly_proper(cl::ClosedLoopTFStandard)
    return degree(cl.p_b) > max(degree(cl.p_a), degree(cl.p_c))
end
strictly_proper(cl::ClosedLoopTF) = strictly_proper(ClosedLoopTFStandard(cl))

"""
 dx
---- = A * x(t) + B * u(t) + C * x(t - ϕ)
 dt
y(t) = D * x(t)
we compute the matrices A, B, C, D here.
only works if strictly proper.
"""
function tf_to_ss(cl::ClosedLoopTFStandard)
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

function simulate(cl::ClosedLoopTF, final_time::Float64; nb_time_points::Int=100)
	cls = ClosedLoopTFStandard(cl)
	
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
