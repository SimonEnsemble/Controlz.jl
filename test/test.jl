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
