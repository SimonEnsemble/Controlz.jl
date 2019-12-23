using Test
using Polynomials

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
    @test tf.numerator == Poly([1], :s)
    @test tf.denominator == Poly([2, 3], :s) # reverse order

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

    ###
    #   subtract
    ###
    @test isapprox(-s, -1.0 * s)
    g = 3 / (s+2) - 2 / (s+2)
    @test isapprox(pole_zero_cancellation(g), 1/(s+2))
    
    ###
    #  divide
    ###
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
    @test k == -1/2
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
    @test isapprox(_z, [0.0, 0.0])
    @test isapprox(_p, [0.0, 8.0, 23.0])
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
    @test k == -1
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
    @test characteristic_polynomial(g_ol) == Poly([10.0, 11, 6, 1], :s)
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

    # first order step response
    K = 4.3
    τ = 2.8
    g = K / (τ * s + 1)
    t, y = simulate(g, unit_step, (0.0, 12.0))
    y_truth = K * (1.0 .- exp.(-t ./ τ))
    @test isapprox(y_truth, y, rtol=0.0001)
    Y = g / s # invert Y(s) way
    t, y = simulate(Y, (0.0, 12.0))
    @test isapprox(y_truth, y, rtol=0.0001)

    # first order impulse response
    t, y = simulate(g, (0.0, 12.0))
    y_truth = K / τ * exp.(-t/τ)
    @test isapprox(y_truth, y, rtol=0.0001)

    # first order ramp input
    a = 2.0 # slope of ramp
    t, y = simulate(g, t -> (t < 0.0) ? 0.0 : a * t, (0.0, 10.0))
    y_truth = K * a * t .+ K * a * τ * (exp.(- t ./ τ) .- 1.0)
    @test isapprox(y_truth, y, rtol=0.0001)
    Y = g * a / s^2
    t, y = simulate(Y, (0.0, 10.0))
    @test isapprox(y_truth, y, rtol=0.0001)

    # first order sinusoidal input
    ω = 2.0
    A = 4.5
    t, y = simulate(g, t -> (t < 0.0) ? 0.0 : A * sin(ω * t), (0.0, 10.0))
    ϕ=atan(-τ*ω)
    y_truth = τ * ω / (1 + (τ*ω)^2) * exp.(-t ./ τ) .+ 1 / sqrt(1+(τ*ω)^2) * sin.(ω*t .+ ϕ)
    y_truth *= K * A
    @test isapprox(y_truth, y, rtol=0.001)
    Y = g * A * ω / (s^2 + ω^2)
    t, y = simulate(Y, (0.0, 10.0))
    @test isapprox(y_truth, y, rtol=0.001)

    # FOPTD
    M = 3.3
    θ = 2.3
    τ = 0.7
    g = K / (τ * s + 1) * exp(-θ * s)
    t, y = simulate(g, t -> (t < 0.0) ? 0.0 : M, (0.0, 10.0))
    y_truth = K * M * (1.0 .- exp.(-(t .- θ) ./ τ))
    y_truth[t .< θ] .= 0.0
    @test isapprox(y_truth, y, rtol=0.001)

    ##
    # second order, overdamped
    K = 4.3
    τ = 2.8
    ξ = 1.2
    g = K / (τ ^ 2 * s ^ 2 + 2 * τ * ξ * s + 1)
    # response to step
    M = 3.3
    t, y = simulate(g, t -> (t < 0.0) ? 0.0 : M, (0.0, 20.0))
    y_truth = 1.0 .- exp.(- ξ / τ * t) .* (
    ξ / sqrt(ξ^2 - 1) * sinh.(sqrt.(ξ^2 - 1) / τ * t) .+ cosh.(sqrt.(ξ^2 - 1) / τ * t))
    y_truth *= K * M
    @test isapprox(y_truth, y, rtol=0.001)
    Y = g * M / s # invert Y(s) way
    t, y = simulate(Y, (0.0, 20.0))
    @test isapprox(y_truth, y, rtol=0.001)
    # response to impulse
    t, y = simulate(g, (0.0, 20.0))
    y_truth = K / τ / sqrt(ξ^2-1) * exp.(-ξ/τ*t) .* sinh.(sqrt(ξ^2-1)/τ*t)
    @test isapprox(y_truth, y, rtol=0.001)

    ##
    # second order, underdamped
    K = 4.3
    τ = 2.8
    ξ = 0.2
    g = K / (τ ^ 2 * s ^ 2 + 2 * τ * ξ * s + 1)
    # impulse response
    t, y = simulate(g, (0.0, 20.0))
    y_truth = K / τ / sqrt(1-ξ^2) * exp.(-ξ/τ*t) .* sin.(sqrt(1-ξ^2)/τ*t)
    @test isapprox(y_truth, y, rtol=0.001)

    ##
    # second order with zeros
    K = 32.0
    τₐ = 2.3
    τ₁ = 1.1
    τ₂ = 2.5
    g = K * (τₐ * s + 1) / (τ₁ * s + 1) / (τ₂ * s + 1)
    # step response
    t, y = simulate(g, unit_step, (0.0, 10.0)) # with explicit input as a function of time
    y_truth = 1 .- (τ₁ - τₐ) / (τ₁ - τ₂) * exp.(-t/τ₁) .- (τ₂ - τₐ) / (τ₂ - τ₁) * exp.(-t/τ₂)
    y_truth *= K
    @test isapprox(y_truth, y, rtol=0.001)
    t, y = simulate(g / s, (0.0, 10.0)) # direct inversion
    @test isapprox(y_truth, y, rtol=0.001)
end
