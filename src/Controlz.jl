module Controlz

using Polynomials
using DifferentialEquations

include("tf.jl")
include("systems.jl")
include("show.jl")
include("sim.jl")
include("input_zoo.jl")
include("viz.jl")

export
    # tf.jl
    TransferFunction, zeros_poles_gain, zeros_poles_k, zero_frequency_gain, evaluate, s, proper, strictly_proper, pole_zero_cancellation,
    # systems.jl
    characteristic_polynomial,
    # show.jl
    # sim.jl
    simulate,
    # input_zoo.jl
    unit_step,
    # viz.jl
    viz_response, viz_poles_and_zeros, viz_nyquist_diagram, viz_root_locus
end
