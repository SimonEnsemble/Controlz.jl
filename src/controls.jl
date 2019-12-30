@doc raw"""
    pc = PController(Kc)

Construct a Proportional (P) controller by specifying the controller gain defined under the following transfer function representation:

$$g_c(s)=K_c$$

# Arguments
* `Kc::Float64`: controller gain

# Example
```julia
pc = PController(1.0) # specify P controller gain
gc = TransferFunction(pc) # construct transfer function with this P-controller gain
```
"""
mutable struct PController
    Kc::Float64
end

TransferFunction(pc::PController) = TransferFunction([pc.Kc], [1.0])

@doc raw"""
    pic = PIController(Kc, τI)

Construct a Proportional-Integral (PI) controller by specifying the controller gain and integral time constant defined under the following transfer function representation:

$$g_c(s)=K_c \left[1+\frac{1}{\tau_I s}\right]$$

# Arguments
* `Kc::Float64`: controller gain
* `τI::Float64`: integral time constant

# Example
```julia
pic = PIController(1.0, 3.0) # specify PI controller params
gc = TransferFunction(pic) # construct transfer function with these PI-controller params
```
"""
mutable struct PIController
    Kc::Float64
    τI::Float64
end

TransferFunction(pic::PIController) = pic.Kc * (1 + 1 / (pic.τI * s))

mutable struct PIDController
    Kc::Float64
    τI::Float64
    τD::Float64
    α::Float64
end

@doc raw"""
    pidc = PIDController(Kc, τI, τD, α=0.0)

Construct a Proportional-Integral-Derivative (PID) controller by specifying the controller gain, integral time constant, derivative time constant, and derivative filter defined under the following transfer function representation:

$$g_c(s)=K_c \left[1+\frac{1}{\tau_I s}+\tau_D s \frac{1}{\alpha \tau_D s + 1}\right]$$

# Arguments
* `Kc::Float64`: controller gain
* `τI::Float64`: integral time constant
* `τD::Float64`: derivative time constant
* `α::Float64`: derivative filter

# Example
```julia
pidc = PIDController(1.0, 3.0, 0.1) # specify PID controller params
gc = TransferFunction(pidc) # construct transfer function with these PID-controller params
```
"""
PIDController(Kc::Float64, τI::Float64, τD::Float64; α::Float64=0.0) = PIDController(Kc, τI, τD, α)

TransferFunction(pidc::PIDController) = pidc.Kc * (1 + 1.0 / (pidc.τI * s) + pidc.τD * s / (pidc.α * pidc.τD * s + 1))
