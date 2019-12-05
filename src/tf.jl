"""
    tf = TransferFunction([1, 2], [3, 5, 8])
    tf = TransferFunction([1, 2], [3, 5, 8], 3.0)

A transfer function representing a linear, time-invariant system.

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
Compute the zeros of a transfer function.
"""
zeros(tf::TransferFunction) = roots(tf.numerator)

"""
Compute the poles of a transfer function.
"""
poles(tf::TransferFunction) = roots(tf.denominator)
