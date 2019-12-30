module Controlz

using Polynomials
using DifferentialEquations
using PyPlot
using Interpolations
using Printf

# hipster plot theme
PyPlot.matplotlib.style.use(normpath(joinpath(pathof(Controlz), "..", "hipster.mplstyle")))

include("tf.jl")
include("systems.jl")
include("sim.jl")
include("viz.jl")
include("controls.jl")
include("show.jl")

export
    # tf.jl
    TransferFunction, zeros_poles_gain, zeros_poles_k, zero_frequency_gain, evaluate, s, proper, strictly_proper, pole_zero_cancellation, zpk_form,
    # systems.jl
    characteristic_polynomial,
    # show.jl
    # sim.jl
    simulate, interpolate,
    # viz.jl
    viz_response, viz_poles_and_zeros, nyquist_diagram, root_locus, bode_plot,
    # controlz.jl
    PController, PIController, PIDController
end
