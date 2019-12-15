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
    # multiply
    ###
    g1 = TransferFunction([1], [3, 1])
    g2 = TransferFunction([2], [5, 1])
    @test isapprox(g1 * g2, zeros_poles_gain([], [-1/3, -1/5], 2))
    g1 = TransferFunction([2], [4, 1], 4.0)
    g2 = TransferFunction([1], [2, 1], 2.0)
    @test isapprox(g1 * g2, TransferFunction([2], [8, 6, 1], 6.0))
    g = TransferFunction([2], [3, 4])
    @test isapprox(5 * g, TransferFunction([10], [3, 4]))
    
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
    #  time delay
    ###
    g = exp(-4.0*s)
    @test g.time_delay == 4.0
    g = 1 / (s + 1) * exp(-6.2 * s)
    @test isapprox(g, TransferFunction([1], [1, 1], 6.2))
    
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
    # reconstruct, make sure the same
    tf2 = zeros_poles_gain(z, p, k)
    @test isapprox(tf, tf2)
    tf2 = zeros_poles_gain(z, p, k, time_delay=3.0)
    @test ! isapprox(tf, tf2)

    ###
    #  evaluate
    ###
    g = TransferFunction([1], [3, 1])
    z, p, k = zeros_poles_gain(g)
    evaluate(g, 0.0) == k
    evaluate(g, 1.0) == 1/3

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
