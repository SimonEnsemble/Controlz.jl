@doc raw"""
    tf = TransferFunction([1, 2], [3, 5, 8])
    tf = TransferFunction([1, 2], [3, 5, 8], 3.0)

Construct a transfer function representing a linear, time-invariant system.

# Example
construct the transfer function:

```math
$G(s) = \frac{4e^{-2.2s}}{2s+1}$
```

```
julia> tf = TransferFunction([4], [2, 1], 2.2)
```

# Attributes
* `numerator::Poly`: the polynomial in the numerator of the transfer function
* `denominator::Poly`: the polynomial in the denominator of the transfer function
* `time_delay::Float64`: the associated time delay
"""
struct TransferFunction
    numerator::Poly
    denominator::Poly
    time_delay::Float64
end

ArrayOfReals = Union{Array{Int64, 1}, Array{Float64, 1}}

TransferFunction(num::ArrayOfReals, den::ArrayOfReals) = 
    TransferFunction(Poly(reverse(num), :s), Poly(reverse(den), :s), 0.0)

TransferFunction(num::ArrayOfReals, den::ArrayOfReals, td::Float64) = 
    TransferFunction(Poly(reverse(num), :s), Poly(reverse(den), :s), td)

"""
    tf_zeros(tf)

Compute the zeros of a transfer function.
"""
tf_zeros(tf::TransferFunction) = roots(tf.numerator)

"""
    tf_poles(tf)

Compute the poles of a transfer function.

# Example
```
julia> tf = Tf([1], [2, 4])
```
"""
tf_poles(tf::TransferFunction) = roots(tf.denominator)

"""
    K = gain(tf)

Compute the gain of a transfer function. The gain is computed by evaluating the transfer function G(s) at s = 0.

# Example
julia> gain(TransferFunction([2], [4, 2])) # 1
julia> gain(TransferFunction([1], [2, 1])) # 1
"""
gain(tf::TransferFunction) = polyval(tf.numerator, 0.0) / polyval(tf.denominator, 0.0)

"""
    standard_K_τ_form(tf)

Convert the transfer function to standard, gain (K)/ time constant (τ) form.

# Example
julia> tf = TransferFunction([2], [4, 2]) # 2 / (4s + 2)
julia> tf = standard_K_τ_form(tf) # 1 / (2s + 1)
"""
function standard_K_τ_form(tf::TransferFunction)
    # multiply by one in a fancy way.
    constant_in_denominator = tf.denominator.a[1]
    # divide by numerator and denominator in this

    return TransferFunction(
        tf.numerator / constant_in_denominator,
        tf.denominator / constant_in_denominator,
        tf.time_delay)
end
