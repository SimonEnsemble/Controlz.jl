@doc raw"""
    gc = PID(Kc, τI, τD, α=0.0)

Construct Proportional-Integral-Derivative (PID) controller transfer function, $g_c(s)$.

$$g_c(s)=K_c \bigl[1+\frac{\tau_I}{s}+\tau_D s \frac{1}{\alpha \tau_D s + 1}\bigr]$$

# Arguments
* `Kc::Float64`: controller gain
* `τI::Float64`: integral time constant
* `τD::Float64`: derivative time constant
* `α::Float64`: derivative filter

# Returns
PID controller transfer function `gc::TransferFunction` that takes error signal as input and outputs the manipulated variable.
"""
function PID(Kc::Float64, τI::Float64, τD::Float64; α::Float64=0.0)
    return Kc * (1 + τI / s + τD * s / (α * τD * s + 1))
end

@doc raw"""
    gc = PI(Kc, τI)

Construct Proportional-Integral (PI) controller transfer function, $g_c(s)$.

$$g_c(s)=K_c \bigl[1+\frac{\tau_I}{s}\bigr]$$

# Arguments
* `Kc::Float64`: controller gain
* `τI::Float64`: integral time constant

# Returns
PI controller transfer function `gc::TransferFunction` that takes error signal as input and outputs the manipulated variable.
"""
PI(Kc::Float64, τI::Float64) = Kc * (1 + τI / s)
