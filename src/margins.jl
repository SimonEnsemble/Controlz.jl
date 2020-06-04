struct Margins
    # critical frequency
    ω_c::Float64
    # gain cross-over frequency
    ω_g::Float64
    # gain margin
    gain_margin::Float64
    # phase margin (radians)
    phase_margin::Float64
end

"""
    margins = gain_phase_margins(g_ol, ω_c_guess=0.001, ω_g_guess=0.001)

compute critical frequency (radians / time), gain crossover frequency (radians / time), 
gain margin, and phase margin (radians) of a closed loop, given its
closed loop transfer function `g_ol::TransferFunction`.

if ω_c or ω_g is not found (i.e. if either are `NaN`), but the `bode_plot` clearly shows 
a critical/gain crossover frequency, adjust `ω_c_guess` or `ω_g_guess` to find the root.

# Example
```
g_ol = 2 * exp(-s) / (5 * s + 1)
margins = gain_phase_margins(g_ol)
margins.ω_c # critical freq. (radians / time)
margins.ω_g # gain crossover freq. (radians / time)
margins.gain_margin # gain margin
margins.phase_margin # phase margin (radians)
```
"""
function gain_phase_margins(g_ol::TransferFunction; ω_c_guess::Float64=0.001, ω_g_guess::Float64=0.001)
    # ∠ G(i ω_c) = -π
    ω_c = NaN
    try
        ω_c = fzero(ω -> angle(evaluate(g_ol, im * ω)) + π, ω_c_guess, atol=2*sqrt(eps()))
    catch da_error
        if isa(da_error, Roots.ConvergenceFailed)
            ω_c = NaN
        else
            println(da_error)
            error("something went wrong when computing ω_c")
        end
    end
    # | G(i ω_g) | = 1
    ω_g = NaN
    try
        ω_g = fzero(ω -> abs(evaluate(g_ol, im * ω)) - 1.0, ω_g_guess, atol=2*sqrt(eps()))
    catch da_error
        if isa(da_error, Roots.ConvergenceFailed)
            ω_g = NaN
        else
            println(da_error)
            error("something went wrong when computing ω_g")
        end
    end
    # gain margin = 1 / | G(i ω_c) |
    gm = 1 / abs(evaluate(g_ol, im * ω_c))
    # phase margin = ∠ G(i ω_g) + π
    pm = π + angle(evaluate(g_ol, im * ω_g))
    return Margins(ω_c, ω_g, gm, pm)
end
