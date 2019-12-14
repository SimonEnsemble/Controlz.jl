module Controlz

using Polynomials
using DifferentialEquations

include("tf.jl")
include("show.jl")
include("sim.jl")
include("input_zoo.jl")

export
    # tf.jl
    TransferFunction, zeros_poles_gain, evaluate, s, proper, strictly_proper,
    # show.jl
    # sim.jl
    simulate,
    # input_zoo.jl
    unit_step
end
