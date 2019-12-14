import Base.*, Base./, Base.+, Base.-, Base.==, Base.^, Base.isapprox

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

ArrayOfReals = Union{Array{Float64, 1}, Array{Int64, 1}}

TransferFunction(num::ArrayOfReals, den::ArrayOfReals) = 
    TransferFunction(Poly(reverse(num), :s), Poly(reverse(den), :s), 0.0)

TransferFunction(num::ArrayOfReals, den::ArrayOfReals, td::Float64) = 
    TransferFunction(Poly(reverse(num), :s), Poly(reverse(den), :s), td)

*(tf1::TransferFunction, tf2::TransferFunction) =
    TransferFunction(tf1.numerator * tf2.numerator,
                     tf1.denominator * tf2.denominator,
                     tf1.time_delay + tf2.time_delay)

*(k::Number, tf::TransferFunction) = TransferFunction(k * tf.numerator, tf.denominator, tf.time_delay)
*(tf::TransferFunction, k::Number) = *(k, tf)
function ^(tf::TransferFunction, x::Int)
    # x< 0 case handled by MethodError
    if x == 0
        return TransferFunction([1], [1])
    else
        tfp = tf
        for i = 2:x
            tfp = tfp * tf
        end
        return tfp
    end
end

function +(tf1::TransferFunction, tf2::TransferFunction)
    if (tf1.time_delay != tf2.time_delay)
        error("cannot add two transfer functions that have different time delays")
    end
    return TransferFunction(tf1.numerator * tf2.denominator + tf1.denominator * tf2.numerator,
                            tf1.denominator * tf2.denominator,
                            tf1.time_delay)
end

+(tf::TransferFunction, x::Number) = tf + TransferFunction([x], [1.0])
+(x::Number, tf::TransferFunction) = +(tf, x)
-(tf::TransferFunction, x::Number) = +(tf, -x)
-(x::Number, tf::TransferFunction) = +(-1 * tf, x)

/(tf1::TransferFunction, tf2::TransferFunction) =
   TransferFunction(tf1.numerator * tf2.denominator,
                    tf1.denominator * tf2.numerator,
                    tf1.time_delay - tf2.time_delay)
/(tf::TransferFunction, x::Number) = TransferFunction(tf.numerator / x, tf.denominator, tf.time_delay)
/(x::Number, tf::TransferFunction) = TransferFunction([x], [1.0], 0.0) / tf

# note: must have polynomails written in same form
==(tf1::TransferFunction, tf2::TransferFunction) = (tf1.numerator == tf2.numerator) && (tf1.denominator == tf2.denominator) && (tf1.time_delay == tf2.time_delay)

# make sure we have a +1 in the denominator
function _standardize(tf::TransferFunction)
    # multiply by one in a fancy way.
    constant_in_denominator = tf.denominator.a[1]
    # divide by numerator and denominator in this

    return TransferFunction(
        tf.numerator / constant_in_denominator,
        tf.denominator / constant_in_denominator,
        tf.time_delay)
end

function isapprox(tf1::TransferFunction, tf2::TransferFunction)
    tf1s = _standardize(tf1)
    tf2s = _standardize(tf2)
    return isapprox(tf1s.time_delay, tf2s.time_delay) && isapprox(tf1s.numerator, tf2s.numerator) && isapprox(tf1s.denominator, tf2s.denominator)
end

_zeros(tf::TransferFunction) = roots(tf.numerator)

_poles(tf::TransferFunction) = roots(tf.denominator)

_gain(tf::TransferFunction) = polyval(tf.numerator, 0.0) / polyval(tf.denominator, 0.0)

"""
    z, p, k = zeros_poles_gain(tf) # compute zeros, poles, gain for a `TransferFunction`

Compute the zeros, poles, and gain of a transfer function. 

* the gain is computed by evaluating the transfer function G(s) at s = 0.
* the zeros are computed as the zeros of the numerator of the transfer function.
* the poles are computed as the zeros of the denominator of the transfer function.

# Example
julia> tf = TransferFunction([1], [4, 1])
julia> z, p, k = zeros_poles_gain(tf) # ([], [-0.25], 1)

---

    tf = zeros_poles_gain(z, p, k, time_delay=0.0) # construct a `TransferFunction` with given zeros, poles, and gain.

Construct a `TransferFunction` by passing an array of the zeros, array of the poles, and a gain.

# Example
julia> tf = zeros_poles_gain([], [-0.25], 1) # 1 / (s + 0.25)
"""
zeros_poles_gain(tf::TransferFunction) = _zeros(tf), _poles(tf), _gain(tf)

# docstring above. this overloaded function is for *constructing* tfs
function zeros_poles_gain(zeros::Array, poles::Array, gain::Union{Float64, Int};
                          time_delay::Union{Int64, Float64}=0.0)
    s = Poly([0, 1], :s)

    # construct numerator polynomial 
    num = Poly([1], :s)
    for z in zeros
        num *= (s - z)
    end
    
    # construct denominator polynomial
    den = Poly([1], :s)
    for p in poles
        den *= (s - p)
    end

    # account for gain
    #  ... gain is not necessary 1.0 at this point...
    current_gain = _gain(TransferFunction(num, den, 0.0))
    num *= gain / current_gain

    return TransferFunction(num, den, time_delay)
end

"""
    evaluate(tf, z)

Evaluate a `TransferFunction`, `tf`, at a particular number `z`.

# Examples
julia> tf = TransferFunction([1], [3, 1])
julia> evaluate(tf, 1.0) # 0.25
julia> evaluate(tf, 2.0 + 3.0im) # also takes imaginary numbers as input
"""
function evaluate(tf::TransferFunction, z::Number)
    return polyval(tf.numerator, z) / polyval(tf.denominator, z) * exp(-tf.time_delay * z)
end

"""
    proper(tf)

Return `true` if transfer function `tf` is proper and `false` otherwise.
"""
proper(tf::TransferFunction) = degree(tf.numerator) <= degree(tf.denominator)

"""
    strictly_proper(tf)

Return `true` if transfer function `tf` is strictly proper and `false` otherwise.
"""
strictly_proper(tf::TransferFunction) = degree(tf.numerator) < degree(tf.denominator)

const s = TransferFunction([1, 0], [1])
