using Test
using Polynomials, DataFrames

push!(LOAD_PATH, joinpath(pwd(), "src")) # run from main dir
using Controlz
    
const s = TransferFunction([1, 0], [1])

@testset "testing TransferFunctions" begin
    ###
    # constructors
    ###
    @test TransferFunction([1], [3, 2], 0.0) == TransferFunction([1], [3, 2])
    tf = TransferFunction([1], [3, 2], 3.2)
    @test tf.time_delay == 3.2
    @test tf.numerator == Polynomial([1], :s)
    @test tf.denominator == Polynomial([2, 3], :s) # reverse order
    # special constructors
    @test first_order_system(5.0, 8.0) == 5 / (8 * s + 1)
    @test second_order_system(5.0, 8.0, 0.1) == 5 / (8^2 * s^2 + 2*8*0.1*s + 1)

    ###
    #  isapprox
    ###
    g1 = 4 * s / (2 * s + 1)
    g2 = 8 * s / (4 * s + 2)
    @test isapprox(g1, g2)
    g2 = 8 * s / (4 * s + 2) * exp(-0.01*s)
    @test ! isapprox(g1, g2)
    @test isapprox((s^2 + 1) / (3 * s + 3), 1/3*(s^2+1)/(s+1))
    @test ! isapprox(s^2 + 1 / (3 * s + 3), s*1/3*(s^2+1)/(s+1))
    @test isapprox((5 * s + 1) / (s ^ 2 + 4 * s + 5), TransferFunction([5, 1], [1, 4, 5]))
    @test isapprox((5 * s + 1) / (s ^ 2 + 4 * s + 5), zeros_poles_k([-1/5], [-2 + im, -2 - im], 5.0, time_delay=0.0))
    
    ###
    # multiply
    ###
    g1 = TransferFunction([1], [3, 1])
    g2 = TransferFunction([2], [5, 1])
    @test isapprox(zero_frequency_gain(g1 * g2), 2.0)
    @test isapprox(g1 * g2, 2/(3*s+1)/(5*s+1))
    @test isapprox(g1 * g2, zeros_poles_k([], [-1/3, -1/5], 2/15))
    g1 = TransferFunction([2], [4, 1], 4.0)
    g2 = TransferFunction([1], [2, 1], 2.0)
    @test isapprox(g1 * g2, TransferFunction([2], [8, 6, 1], 6.0))
    g = TransferFunction([2], [3, 4])
    @test isapprox(5 * g, TransferFunction([10], [3, 4]))

    # gain
    @test isapprox(zero_frequency_gain(5 / (s + 1)), 5.0)
    @test isapprox(zero_frequency_gain(5 / (s + 1) / (4 * s + 2)), 5.0 / 2)
    @test isapprox(zero_frequency_gain(5 * s / (s + 1) / s), 5.0)

    ###
    #  add
    ###
    g1 = 3 / (s + 2)
    g2 = 1 / (s + 4)
    @test isapprox(g1 + g2, (4*s+14) / (s^2 + 6*s + 8))

    ###
    #   subtract
    ###
    @test isapprox(-s, -1.0 * s)
    g = 3 / (s+2) - 2 / (s+2)
    @test isapprox(pole_zero_cancellation(g), 1/(s+2))
    @test isapprox(g, 1/(s+2))
    
    ###
    #  divide
    ###
    g1 = TransferFunction([2], [4, 1], 4.0)
    g2 = TransferFunction([1], [2, 1], 2.0)
    @test isapprox(g1 / g2, TransferFunction([4, 2], [4, 1], 2.0))
    @test isapprox(g1 / 2, TransferFunction([1], [4, 1], 4.0))
    @test isapprox(5 / (s+3), TransferFunction([5], [1, 3]))
    @test isapprox(33 / (2*s+3), TransferFunction([33], [2, 3]))
    @test isapprox((s+3) / ((s - 2) * (2*s-3)), TransferFunction([1, 3], [2, -7, 6]))

    ###
    #  power
    ###
    @test isapprox(s ^ 2, s * s)
    @test isapprox(s ^ 4, s * s * s * s)
    @test isapprox(s ^ 1, s)
    g = TransferFunction([2], [3, 4])
    @test isapprox(g * s ^ 0, g)
    @test_throws MethodError s ^ -1

    ###
    #  tf algebra
    ###
    τ = 1.4
    K = 2.8
    Kc = 3.3
    Kd = 4.5
    τI = 9.0
    τd = 8.2
    gd = Kd / (τd * s + 1)
    gp = K / (τ * s + 1)
    # P-control regulator first-order step response
    gc = Kc
    y_ovr_d = gd / (1 + gc * gp)
    @test isapprox(y_ovr_d, TransferFunction(Kd/(1+Kc*K) * [τ, 1], [τd*τ/(1+Kc*K), (τ+τd+Kc*K*τd)/(1+Kc*K), 1]))
    # PI-control servo first-order step response
    gc = Kc * (τI * s + 1) / (τI * s)
    y_ovr_y_sp = gc * gp / (1 + gc * gp)
    y_ovr_y_sp = pole_zero_cancellation(y_ovr_y_sp)
    @test isapprox(y_ovr_y_sp, TransferFunction([τI, 1], [τI*τ/(Kc*K), τI*(1+Kc*K)/(Kc*K), 1]))

    ###
    #  time delay
    ###
    g = exp(-4.0*s)
    @test g.time_delay == 4.0
    g = 1 / (s + 1) * exp(-6.2 * s)
    @test isapprox(g, TransferFunction([1], [1, 1], 6.2))
    @test isapprox(exp(-3*s), exp(-3.0*s))
    @test isapprox(exp(-s), TransferFunction([1], [1], 1.0))
    
    ###
    # zeros, poles, gain
    ###
    g = TransferFunction([3], [4, 3, 4, 1])
    z, p, k = zeros_poles_gain(g)
    @test k == 3
    g = TransferFunction([3], [4, 3, 4, 2])
    z, p, k = zeros_poles_gain(g)
    @test k == 3/2
    # G(s) = s + 2 / (s^2 - 4)
    tf = TransferFunction([1, 2], [1, 0, -4])
    z, p, k = zeros_poles_gain(tf)
    @test isapprox(k, -1/2, atol=1e-6)
    @test z == [-2.0]
    @test isapprox(p, [-2.0, 2.0])

    ###
    #  zeros, poles, k
    ###
    z = [-3.2, 4.0]
    p = [-3.0, -3.0, 0.0, 2.0]
    k = 23.2
    g = zeros_poles_k(z, p, k)
    _z, _p, _k = zeros_poles_k(g)
    sort!(_z)
    sort!(_p)
    @test isapprox(z, _z)
    @test isapprox(p, _p)
    @test isapprox(k, _k)
    g2 = zeros_poles_k(z, p, k, time_delay=3.0)
    @test ! isapprox(g, g2)
    @test isapprox(g2, g * exp(-3*s))
    g = 45.0 * s * s / (s * (s-23.0) * (s-8.0))
    _z, _p, _k = zeros_poles_k(g)
    sort!(_z)
    sort!(_p)
    @test isapprox(_k, 45.0)
    @test isapprox(_z, [0.0])
    @test isapprox(_p, [8.0, 23.0])
    p = [-2+im, -2-im] # try with complex numbers.
    z = [-1/5]
    @test isapprox(zeros_poles_k(z, p, 5.0), (5*s+1) / (s^2+4*s+5))
    @test ! isapprox(zeros_poles_k(z, p, 5.0, time_delay=1.0), (5*s+1) / (s^2+4*s+5))
    @test isapprox(zeros_poles_k(z, p, 5.0, time_delay=1.0), (5*s+1) / (s^2+4*s+5) * exp(-s))
    g = (6s^2+18*s+12)/(2*s^3+10*s^2+16*s+12)
    z, p, k = zeros_poles_k(g)
    sort!(z)
    @test isapprox(z, [-2.0, -1.0])
    @test isapprox(p, [-3.0, -1.0-im, -1+im])


    ###
    #  evaluate
    ###
    g = TransferFunction([1], [3, 1]) # 1 / (3s+1)
    z, p, k = zeros_poles_gain(g)
    @test evaluate(g, 0.0) == k
    @test evaluate(g, 1.0) == 1/4
    @test evaluate(g, 0.0) == 1.0
    # https://www.mathworks.com/help/control/ref/evalfr.html
    g = (s - 1) / (s ^ 2 + s + 1)
    z, p, k = zeros_poles_gain(g)
    @test isapprox(k, -1, atol=1e-6)
    @test evaluate(g, 0.0) == k
    @test isapprox(evaluate(g, 1+im), 0.2308 + 0.1538im, atol=0.0001)
    g = 1.0 / (s^2 + 2*s + 1)
    @test isapprox(evaluate(g, 0.1*im), 0.9705 - 0.1961im, atol=0.0001)
    

    ###
    #  properness
    ###
    tf = 1 / (2 * s + 1)
    @test proper(tf)
    @test strictly_proper(tf)
    tf = 3 * s / (2 * s + 1)
    @test proper(tf)
    @test ! strictly_proper(tf)
    tf = (3 * s * s) / (2 * s + 1)
    @test ! proper(tf)
    @test ! strictly_proper(tf)

    ###
    #  zeros, poles cancellation
    ###
    @test isapprox(pole_zero_cancellation((s-1)^2 / (s-1)), s-1)
    @test isapprox(pole_zero_cancellation((s-1) / (s-1)^2), 1/(s-1))
    @test isapprox(pole_zero_cancellation((2*s-2.0) / (s-1)^2), 2/(s-1))
    @test isapprox(pole_zero_cancellation(s * (9*s+3) / (s * (3*s-3))), (9*s+3)/(3*s-3))
    @test isapprox(pole_zero_cancellation(s^5/s), s^4)
    @test isapprox(pole_zero_cancellation(s^5*(s-1)/s/(s-1)), s^4)
    @test isapprox(pole_zero_cancellation(s^5*(s-1)/s/(s-1)*(s+2)), (s+2)*s^4)
    g = s * (s+1) / ((s+3) * s * (s+1) ^ 2)
    @test isapprox(pole_zero_cancellation(g), 1 / ((s+3) * (s+1)))
    g = (s+1) / ((s^2-2*s+2) * (s+1) * s) # some imaginary poles
    @test isapprox(pole_zero_cancellation(g), 1/(s^2-2*s+2)/s)

    ###
    #  characteristic eqn.
    ###
    g_ol = 4 / (s + 3) / (s + 2) / (s + 1)
    @test isapprox(characteristic_polynomial(g_ol), Polynomial([10.0, 11, 6, 1], :s))

    ###
    #   order
    ###
    g = (s + 1) / (s + 1) ^ 2 # make sure it cancels the pole and zero.
    @test system_order(g) == (0, 1)
    @test system_order(first_order_system(1.0, 2.0)) == (0, 1)
    @test system_order(second_order_system(1.0, 2.0, 1.9)) == (0, 2)

    ###
    #  time constant, damping coefficient
    ###
    g = 4 / (6 * s + 2)
    @test isapprox(time_constant(g), 3.0)
    @test_throws ErrorException damping_coefficient(g)
    g = 1.0 / (8 * s^2 + 0.8 * s + 2)
    @test isapprox(time_constant(g), 2.0)
    @test_throws ErrorException time_constant(s / (s+4)) # only makes sense with first, second order...
    g = 1.0 / (8 * s^2 + 0.8 * s + 2)
    @test isapprox(damping_coefficient(g), 0.1) # 0.1
    τ = 4.3
    ξ = 1.2
    g = second_order_system(20.0, τ, ξ)
    @test isapprox(time_constant(g), τ)
    @test isapprox(damping_coefficient(g), ξ)
end

@testset "testing Simulation" begin
    # from http://web.mit.edu/2.14/www/Handouts/StateSpace.pdf
    tf = (13*s+26) / (1*s*s*s+7*s*s+19*s+13)
    A, B, C = Controlz.tf_to_ss(tf)
    @test isapprox(A, [0.0 1.0 0.0; 0  0 1; -13 -19 -7])
    @test isapprox(B, [0.0, 0.0, 1.0])
    @test isapprox(C, [26.0 13.0 0.0])
    
    # from https://www.engr.mun.ca/~millan/Eng6825/canonicals.pdf
    tf = TransferFunction([1, 3], [1, 3, 2])
    A, B, C = Controlz.tf_to_ss(tf)
    @test isapprox(A, [0.0 1.0; -2.0 -3.0])
    @test isapprox(B, [0.0, 1.0])
    @test isapprox(C, [3.0 1.0])

    # from https://people.kth.se/~demirel/State_Space_Representation_of_Transfer_Function_Systems.pdf
    tf = TransferFunction([1], [1, 6, 11, 6])
    A, B, C = Controlz.tf_to_ss(tf)
    @test isapprox(A, [0.0 1.0 0.0; 0.0 0.0 1.0; -6.0 -11.0 -6.0])
    @test isapprox(B, [0.0, 0.0, 1.0])
    @test isapprox(C, [1.0 0.0 0.0])
    
    tf = TransferFunction([1, 3, 3], [1, 2, 1])
    A, B, C = Controlz.tf_to_ss(tf)
    @test isapprox(A, [0.0 1.0; -1.0 -2.0])
    @test isapprox(B, [0.0, 1.0])
    @test isapprox(C, [2.0 1.0])

    # some inputs
    # L[cos(at)] = s/(s^2+a^2)
    a = 2.3
    U = s / (s^2+a^2)
    u_data = simulate(U, 12.0)
    _cosat(t::Float64) = t < 0.0 ? 0.0 : cos(a*t)
    @test isapprox(u_data[:, :output], _cosat.(u_data[:, :t]), rtol=0.001)
    # L[shifted step] = e^{-a*s}/s
    U = exp(-a*s) / s
    u_data = simulate(U, 12.0)
    _shifted_step(t::Float64) = t < a ? 0.0 : 1.0
    @test isapprox(u_data[:, :output], _shifted_step.(u_data[:, :t]), rtol=0.001)
    # L[t sin(at)] = 2as/(s^2+a^2)^2
    a = 2.1
    U = 2 * a * s / (s^2 + a^2)^2
    u_data = simulate(U, 12.0)
    u_truth = [t_i > 0.0 ? t_i * sin(a * t_i) : 0.0 for t_i in u_data[:, :t]]
    isapprox(u_data[:, :output], u_truth, rtol=0.01)
    # L{exp[-3(t-2)]S(t-2)} = e^{-2s} / (s+3)
    U = exp(-2*s) / (s + 3)
    u_data = simulate(U, 12.0)
    u_truth = [t_i > 2.0 ? exp(-3*(t_i - 2)) : 0.0 for t_i in u_data[:, :t]]
    isapprox(u_data[:, :output], u_truth, rtol=0.001)

    # first order step response
    K = 4.3
    τ = 2.8
    g = K / (τ * s + 1)
    u = 1 / s
    data = simulate(g * u, 12.0)
    y_truth = K * (1.0 .- exp.(-data[:, :t] ./ τ))
    y_truth[data[:, :t] .< 0.0] .= 0.0
    @test isapprox(y_truth, data[:, :output], rtol=0.0001)

    # first order impulse response
    data = simulate(g, 12.0)
    y_truth = K / τ * exp.(-data[:, :t]/τ)
    y_truth[data[:, :t] .< 0.0] .= 0.0
    @test isapprox(y_truth, data[:, :output], rtol=0.0001)

    # first order ramp input
    a = 2.0 # slope of ramp
    Y = g * a / s^2
    data = simulate(Y, 10.0)
    y_truth = K * a * data[:, :t] .+ K * a * τ * (exp.(- data[:, :t] ./ τ) .- 1.0)
    y_truth[data[:, :t] .< 0.0] .= 0.0
    @test isapprox(y_truth, data[:, :output], rtol=0.0001)

    # first order sinusoidal input
    ω = 2.0
    A = 4.5
    Y = g * A * ω / (s^2 + ω^2)
    data = simulate(Y, 10.0)
    ϕ = atan(-τ*ω)
    y_truth = τ * ω / (1 + (τ*ω)^2) * exp.(-data[:, :t] ./ τ) .+ 1 / sqrt(1+(τ*ω)^2) * sin.(ω*data[:, :t] .+ ϕ)
    y_truth *= K * A
    y_truth[data[:, :t] .< 0.0] .= 0.0
    @test isapprox(y_truth, data[:, :output], rtol=0.001)

    # FOPTD
    M = 3.3
    θ = 2.3
    τ = 0.7
    K = 9.0
    g = K / (τ * s + 1) * exp(-θ * s)
    u = M / s
    data = simulate(g * u, 10.0)
    y_truth = K * M * (1.0 .- exp.(-(data[:, :t] .- θ) ./ τ))
    y_truth[data[:, :t] .< θ] .= 0.0
    @test isapprox(y_truth, data[:, :output], atol=0.001)

    ##
    # second order, overdamped
    K = 4.3
    τ = 2.8
    ξ = 1.2
    g = K / (τ ^ 2 * s ^ 2 + 2 * τ * ξ * s + 1)
    # response to step
    M = 3.3
    Y = g * M / s 
    data = simulate(Y, 20.0)
    y_truth = 1.0 .- exp.(- ξ / τ * data[:, :t]) .* (
        ξ / sqrt(ξ^2 - 1) * sinh.(sqrt.(ξ^2 - 1) / τ * data[:, :t]) .+ cosh.(sqrt.(ξ^2 - 1) / τ * data[:, :t]))
    y_truth *= K * M
    y_truth[data[:, :t] .< 0.0] .= 0.0
    @test isapprox(y_truth, data[:, :output], rtol=0.001)
    # response to impulse
    data = simulate(g, 20.0)
    y_truth = K / τ / sqrt(ξ^2-1) * exp.(-ξ/τ*data[:, :t]) .* sinh.(sqrt(ξ^2-1)/τ*data[:, :t])
    y_truth[data[:, :t] .< 0.0] .= 0.0
    @test isapprox(y_truth, data[:, :output], rtol=0.001)

    ##
    # second order, underdamped
    K = 4.3
    τ = 2.8
    ξ = 0.2
    g = K / (τ ^ 2 * s ^ 2 + 2 * τ * ξ * s + 1)
    # impulse response
    data = simulate(g, 20.0)
    y_truth = K / τ / sqrt(1-ξ^2) * exp.(-ξ/τ*data[:, :t]) .* sin.(sqrt(1-ξ^2)/τ*data[:, :t])
    y_truth[data[:, :t] .< 0.0] .= 0.0
    @test isapprox(y_truth, data[:, :output], rtol=0.001)

    ##
    # second order with zeros
    K = 32.0
    τₐ = 2.3
    τ₁ = 1.1
    τ₂ = 2.5
    g = K * (τₐ * s + 1) / (τ₁ * s + 1) / (τ₂ * s + 1)
    # step response
    data = simulate(g / s, 10.0) # direct inversion
    y_truth = 1 .- (τ₁ - τₐ) / (τ₁ - τ₂) * exp.(-data[:, :t]/τ₁) .- (τ₂ - τₐ) / (τ₂ - τ₁) * exp.(-data[:, :t]/τ₂)
    y_truth *= K
    y_truth[data[:, :t] .< 0.0] .= 0.0
    @test isapprox(y_truth, data[:, :output], rtol=0.001)

    # inerpolate
    t = [0.0, 1.0, 2.0]
    y = [0.0, 5.0, 10.0]
    data = DataFrame(t=t, output=y)
    @test isapprox(interpolate(data, 1.0), 5.0)
    @test isapprox(interpolate(data, 1.2), 1.2*5)
    τ = 3.45
    g = 1 / (τ * s + 1)
    data = simulate(g / s, 10.0)
    @test isapprox(interpolate(data, τ), 0.632, atol=0.002)
end

@testset "testing Controls" begin
    @test isapprox(TransferFunction(PIController(1.3, 2.0)), 1.3 * (1+1/(2*s)))
    @test isapprox(TransferFunction(PIDController(1.3, 2.0, 0.0)), 1.3 * (1+1/(2*s)))
    @test isapprox(TransferFunction(PIDController(1.3, 2.0, 0.1)), 1.3 * (1+1/(2*s) + 0.1*s))
end

@testset "testing margins" begin
    g_ol = 2 * exp(-s) / (5 * s + 1)
    m = gain_phase_margins(g_ol)
    @test isapprox(m.ω_c, 1.6887, atol=0.0001)
    @test isapprox(m.ω_g, 0.3464, atol=0.0001)
    @test isapprox(m.gain_margin, 4.2512, atol=0.0001)
    @test isapprox(m.phase_margin, 100.1535 / 180 * π, atol=0.0001)

    g_ol = TransferFunction([0.25], [1.0, 2, 1, 1])
    m = gain_phase_margins(g_ol)
    @test isapprox(m.ω_c, 1.0, atol=0.0001)
    @test isnan(m.ω_g)
    @test isapprox(m.gain_margin, 4.0, atol=0.0001)
    @test isnan(m.phase_margin)
    
    g_ol = 20*(1+1/(3*s)) * 4 /(s+4)/(s+6)*exp(-0.001*s)/(s+2)
    m = gain_phase_margins(g_ol)
    @test isapprox(m.ω_c, 6.32, atol=0.01)
    @test isapprox(m.gain_margin, 5.397, atol=0.01)
end

@testset "closed loop stuff" begin
    g_d =  3 / (s + 1) * exp(-5 * s)
    g_ol = 18 / (100 * s + 1) * exp(-s)
    cl = ClosedLoopTransferFunction(g_d, g_ol)
    
    # test conversion into standard form
    cls = Controlz.CLTFStandard(cl)
    @test cls.p_a ≈ Polynomial([3, 300], :s)
    @test cls.p_b ≈ Polynomial([1, 101, 100], :s)
    @test cls.p_c ≈ Polynomial([18, 18], :s)
    @test cls.ϕ ≈ 1.0
	@test cls.θ ≈ 5.0

    @test Controlz.order(cl) == 2
    @test Controlz.strictly_proper(cl)
    @test ! strictly_proper(ClosedLoopTransferFunction(g_d, g_ol * (4 * s + 1)))

    # algebra
    cl2 = cl * 4 / s
    @test cl2.g_ol == cl.g_ol
    @test cl2.top == cl.top * 4 / s

    # run without a time delay to compare to other simulate function.
    gp = 3 / (s + 2)
    gc = TransferFunction(PIController(1.0, 3.0))
    data_old = simulate(gp * gc / (1 + gp * gc), 10.0)
    data_new = simulate(ClosedLoopTransferFunction(gp * gc, gp * gc), 10.0)
    @test isapprox(data_old[:, :output], data_new[:, :output], atol=0.0001)
end

@test "example notebook (also to generate images for docs)" begin
    include(joinpath("..", "examples", "examples.jl"))
    @test true
end
