@doc raw"""
    g = first_order_system(K, τ)

construct a first-order transfer function with gain `K` and time constant `τ`:

$$g(s)=\frac{K}{\tau s+1}$$

# example
```julia
K = 1.0
τ = 3.0
g = first_order_system(K, τ) # 1 / (3 * s + 1)
```

# returns
* `g::TransferFunction`: the first order transfer function. well, (0, 1) order.
"""
first_order_system(K::Float64, τ::Float64) = K / (τ * s + 1)

@doc raw"""
    g = second_order_system(K, τ, ξ)

construct a second-order transfer function with gain `K`, time constant `τ`, and damping coefficient `ξ`:

$$g(s)=\frac{K}{\tau^2 s^2 + 2\tau \xi s +1}$$

# example
```julia
K = 1.0
τ = 2.0
ξ = 0.1
g = second_order_system(K, τ, ξ) # 1 / (4 * s^2 + 0.4 * s + 1)
```

# returns
* `g::TransferFunction`: the second order transfer function. well, (0, 2) order.
"""
second_order_system(K::Float64, τ::Float64, ξ::Float64) = K / (τ^2 * s^2 + 2 * τ * ξ * s + 1)

@doc raw"""
    τ = time_constant(g)

compute the time constant τ of an order (0, 1) or order (0, 2) transfer function.

order (0, 1) representation:

$$g(s)=\frac{K}{\tau s+1}$$

order (0, 2) representation:

$$g(s)=\frac{K}{\tau^2 s^2 + 2\tau \xi s +1}$$

# returns
`τ::Float64`: the time constant.

# examples
```julia
g = 4 / (6 * s + 2)
time_constant(g) # 3.0

g = 1.0 / (8 * s^2 + 0.8 * s + 2)
time_constant(g) # 2.0
```
"""
function time_constant(g::TransferFunction)
    if system_order(g) == (0, 1)
        # xx / (a s + b) ==> τ = a / b
        return g.denominator[1] / g.denominator[0]
    elseif system_order(g) == (0, 2)
        # xx / (a s^2 + b s + c) ==> τ = sqrt(a/c)
        return sqrt(g.denominator[2] / g.denominator[0])
    else
        error("`time_constant` only supports order (0, 1) and (0, 2) transfer functions.")
    end
end

@doc raw"""
    ξ = damping_coefficient(g)

compute the damping coefficient ξ of an order (0, 2) transfer function.

order (0, 2) representation:

$$g(s)=\frac{K}{\tau^2 s^2 + 2\tau \xi s +1}$$

# returns
`ξ::Float64`: the damping coefficient

# examples
```julia
g = 1.0 / (8 * s^2 + 0.8 * s + 2)
damping_coefficient(g) # 0.1
```
"""
function damping_coefficient(g::TransferFunction)
    if system_order(g) != (0, 2)
        error("`damping_coefficient` only pertains to (0, 2) transfer functions.")
    end
    τ = time_constant(g)
    # xx / (a s² + b s + c) then 2τξ = b / c
    τξ2 = g.denominator[1] / g.denominator[0]
    return τξ2 / (2 * τ)
end
