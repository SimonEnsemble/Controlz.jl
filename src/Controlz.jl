module Controlz

using Polynomials
using DifferentialEquations
using Interpolations
using Printf
using Roots
using DataFrames
using CairoMakie

# TODO set plot theme
include("beavs_theme.jl")

include("tf.jl")
include("special_tfs.jl")
include("systems.jl")
include("sim.jl")
include("viz.jl")
include("controls.jl")
include("margins.jl")
include("closed_loops.jl")
include("show.jl")

export
    # tf.jl
    TransferFunction, zeros_poles_gain, zeros_poles_k, zero_frequency_gain, evaluate, s, proper, strictly_proper, pole_zero_cancellation, zpk_form, system_order,
    # special_tfs.jl
    first_order_system, second_order_system, time_constant, damping_coefficient,
    # systems.jl
    characteristic_polynomial,
    # show.jl
    # sim.jl
    simulate, interpolate,
    # viz.jl
    viz_response, viz_poles_and_zeros, nyquist_diagram, root_locus, bode_plot, mk_gif,
    # controlz.jl
    PController, PIController, PIDController,
    # margins.jl
    gain_phase_margins,
    # closedloops.jl
    ClosedLoopTransferFunction
end
