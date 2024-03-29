import Base.*, Base./, Base.+, Base.-, Base.==, Base.^, Base.exp, Base.isapprox

@doc raw"""
    tf = TransferFunction([1, 2], [3, 5, 8])
    tf = TransferFunction([1, 2], [3, 5, 8], 3.0)

construct a transfer function representing a linear, time-invariant system.

# example
to construct the transfer function

```math
G(s) = \frac{4e^{-2.2s}}{2s+1}
```

in Julia:

```jldoctest
tf = TransferFunction([4], [2, 1], 2.2)
# output
    4.0
----------- e^(-2.2*s)
2.0*s + 1.0
```

# attributes
* `numerator::Polynomial{Float64, :s}`: the polynomial in the numerator of the transfer function
* `denominator::Polynomial{Float64, :s}`: the polynomial in the denominator of the transfer function
* `time_delay::Float64`: the associated time delay
"""
struct TransferFunction
    numerator::Polynomial{Float64, :s}
    denominator::Polynomial{Float64, :s}
    time_delay::Float64
end

TransferFunction(num::Union{Vector{Int}, Vector{Float64}}, 
                 den::Union{Vector{Int}, Vector{Float64}}, 
                 td::Union{Int, Float64}) = 
    TransferFunction(Polynomial(reverse(convert(Vector{Float64}, num)), :s), 
                     Polynomial(reverse(convert(Vector{Float64}, den)), :s), 
                     convert(Float64, td))

TransferFunction(num::Union{Vector{Int}, Vector{Float64}}, 
                 den::Union{Vector{Int}, Vector{Float64}}) =
    TransferFunction(num, den, 0.0)

const s = TransferFunction([1, 0], [1])

function *(tf1::TransferFunction, tf2::TransferFunction)
    g = TransferFunction(tf1.numerator * tf2.numerator,
                         tf1.denominator * tf2.denominator,
                         tf1.time_delay + tf2.time_delay)
    return pole_zero_cancellation(g)
end

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

"""
    tf_time_delay = exp(- 3 * s)

Conveniently create a time delay by exp(- θ * s).

# Example
g = 1 / (s + 1) * exp(-2.0 * s) # introduce time delay of 2.0
"""
function exp(tf::TransferFunction)
    # if numerator is zero, then the value is 1.0 and no time delay is introduced
    if degree(tf.numerator) == -1
        return TransferFunction([1], [1])
    end
    if ! ((degree(tf.numerator) == 1) && (degree(tf.denominator) == 0) && (tf.time_delay == 0.0))
        error("exp only functions with exp(-3.0 * s) for example, to introduce time delays")
    end
    return TransferFunction([1.0], [1.0], -tf.numerator[1])
end

function +(tf1::TransferFunction, tf2::TransferFunction)
    if (tf1.time_delay != tf2.time_delay)
        error("""we cannot add two transfer functions with two different time delays.

        if you are trying to construct a closed-loop transfer function with a time delay,
        instead use `ClosedLoopTransferFunction`.""")
    end
    g = TransferFunction(tf1.numerator * tf2.denominator + tf1.denominator * tf2.numerator,
                            tf1.denominator * tf2.denominator,
                            tf1.time_delay)
    return pole_zero_cancellation(g)
end

+(tf::TransferFunction, x::Number) = tf + TransferFunction([x], [1.0])
+(x::Number, tf::TransferFunction) = +(tf, x)
-(tf::TransferFunction, x::Number) = +(tf, -x)
-(x::Number, tf::TransferFunction) = +(-1 * tf, x)
-(tf::TransferFunction) = -1.0 * tf
-(tf1::TransferFunction, tf2::TransferFunction) = +(tf1, -1*tf2)

function /(tf1::TransferFunction, tf2::TransferFunction)
    g = TransferFunction(tf1.numerator * tf2.denominator,
                        tf1.denominator * tf2.numerator,
                        tf1.time_delay - tf2.time_delay)
    return pole_zero_cancellation(g)
end
/(tf::TransferFunction, x::Number) = TransferFunction(tf.numerator / x, tf.denominator, tf.time_delay)
/(x::Number, tf::TransferFunction) = TransferFunction([x], [1.0], 0.0) / tf

# note: must have polynomails written in same form
==(tf1::TransferFunction, tf2::TransferFunction) = (tf1.numerator == tf2.numerator) && (tf1.denominator == tf2.denominator) && (tf1.time_delay == tf2.time_delay)

@doc raw"""
    tf = zpk_form(tf)

write transfer function `tf` in zeros, poles, k-factor form:

$$g(s)=k\dfrac{\Pi_j (s-z_j)}{\Pi_j (s-p_j)}$$
where $z_j$ is zero $j$, $p_j$ is pole $j$, and $k$ is a constant factor (not equal to the zero-frequency gain) that uniquely specifies the transfer function.

this is achieved by multiplying by 1.0 in a fancy way such that the highest power of $s$ in the denominator is associated with a coefficient of $1$.

# Example

```jldoctest
g = 8.0 / (2 * s^2 + 3 * s + 4)
g_zpk = zpk_form(g)
# output
         4.0
---------------------
1.0*s^2 + 1.5*s + 2.0
```
"""
function zpk_form(tf::TransferFunction)
    coeff_highest_power_denom = tf.denominator[end]

    return TransferFunction(
        tf.numerator / coeff_highest_power_denom,
        tf.denominator / coeff_highest_power_denom,
        tf.time_delay)
end

function isapprox(tf1::TransferFunction, tf2::TransferFunction; atol=0.0)
    # to directly compare numerators and denominators, put in zpk form
    tf1s = zpk_form(tf1)
    tf2s = zpk_form(tf2)
    return isapprox(tf1s.time_delay, tf2s.time_delay, atol=atol) && isapprox(
	    tf1s.numerator.coeffs, tf2s.numerator.coeffs, atol=atol) && isapprox(tf1s.denominator.coeffs, tf2s.denominator.coeffs, atol=atol)
end

_zeros(tf::TransferFunction) = roots(tf.numerator)

_poles(tf::TransferFunction) = roots(tf.denominator)

_k(tf::TransferFunction) = tf.numerator[end] / tf.denominator[end]

@doc raw"""
    K = zero_frequency_gain(tf)

compute the (signed) zero frequency gain of a transfer function $g(s)$, which is:

```math
K := \lim_{s\rightarrow 0} G(s)
```
    
the zero-frequency gain "represents the ratio of the steady state value of the output with respect to a step input" [source](http://www.cds.caltech.edu/~murray/books/AM05/pdf/am06-xferfcns_16Sep06.pdf)

# example

```jldoctest
g = 5 / (3 * s + 1)
K = zero_frequency_gain(g)
# output
5.0
```

# arguments
* `tf::TransferFunction`: the transfer function

# returns
* `K::Float64`: the zero-frequency gain of the transfer function
"""
function zero_frequency_gain(tf::TransferFunction)
    tf_simplified = pole_zero_cancellation(tf)
    return evaluate(tf_simplified, 0.0)
end

@doc raw"""
    # compute the zeros, poles, and k-factor of a transfer function
    z, p, k = zeros_poles_k(tf)
    # construct a transfer function from its zeros, poles, and k-factor
    tf = zeros_poles_k(z, p, k, time_delay=0.0)

the representation of a transfer function in this context is:

$$g(s)=k\dfrac{\Pi_j (s-z_j)}{\Pi_j (s-p_j)}$$
where $z_j$ is zero $j$, $p_j$ is pole $j$, and $k$ is a constant factor (not equal to the zero-frequency gain) that uniquely specifies the transfer function.

* the zeros are the zeros of the numerator of the transfer function.
* the poles are the zeros of the denominator of the transfer function.
"""
zeros_poles_k(tf::TransferFunction) = _zeros(tf), _poles(tf), _k(tf)

# docstring above. this overloaded function is for *constructing* tfs
function zeros_poles_k(zeros::Array, poles::Array, k::Union{Float64, Int};
                       time_delay::Union{Int64, Float64}=0)
    # construct numerator polynomial 
    top = length(zeros) == 0 ? Polynomial(1.0, :s) : fromroots(zeros, var=:s)
    
    # construct denominator polynomial
    bottom = length(poles) == 0 ? Polynomial(1.0, :s) : fromroots(poles, var=:s)

    # owing to numerical errors, the coeff's could be imaginary with tiny imaginary parts
    #  let's check that the imaginary parts are tiny, then convert the coefficients to real.
    @assert maximum(abs.(imag.(top.coeffs))) < 1e-6
    @assert maximum(abs.(imag.(bottom.coeffs))) < 1e-6
    top = Polynomial(real.(top.coeffs), :s) # reconstruct a polynomial with real coeffs
    bottom = Polynomial(real.(bottom.coeffs), :s)
    
    return TransferFunction(k * top, bottom, time_delay)
end

@doc raw"""
    z, p, gain = zeros_poles_gain(tf)

Compute the zeros, poles, and zero-frequency gain of a transfer function.

* the zeros are the zeros of the numerator of the transfer function.
* the poles are the zeros of the denominator of the transfer function.
* the zero-frequency gain is the transfer function evaluated at $s=0$
"""
zeros_poles_gain(tf::TransferFunction) = _zeros(tf), _poles(tf), zero_frequency_gain(tf)

"""
    evaluate(tf, z)

evaluate a `TransferFunction`, `tf`, at a particular number `z` (could be complex).

# example
```jldoctest eval
tf = TransferFunction([1], [3, 1])
evaluate(tf, 1.0)
# output
0.25
```
"""
function evaluate(tf::TransferFunction, z::Number)
    return tf.numerator(z) / tf.denominator(z) * exp(-tf.time_delay * z)
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


"""
    tf = pole_zero_cancellation(tf, verbose=false, digits=8)

find (pole, zero) pairs such that pole = zero and return a new transfer function with those pairs cancelled.
this is achieved by comparing the poles and zeros with `isapprox`, with poles and zeros rounded to 
`digits` digits (also applies to reconstruction).

# arguments
* `tf::TransferFunction`: the transfer function
* `verbose::Bool=false`: print off which poles, zeros are cancelled.
* `digits::Int`: number of digits to round poles and zeros to, for (i) cancelling and (ii) reconstruction.

# example
```jldoctest
pole_zero_cancellation(s * (s - 1) / (s * (s + 1)))
# output
1.0*s - 1.0
-----------
1.0*s + 1.0
```
"""
function pole_zero_cancellation(tf::TransferFunction; verbose::Bool=false, digits::Int=6)
    # compute poles, zeros, and k-factor of the transfer function
    zs, ps, k = zeros_poles_k(tf)
    zs = round.(zs, digits=digits)
    ps = round.(ps, digits=digits)

    # store in these boolean arrays whether we will cancel them or not
    canceled_zeros = [false for i = 1:length(zs)]
    canceled_poles = [false for i = 1:length(ps)]

    # loop through poles and zeros. if they are equal, cancel!
    for (i_z, z) in enumerate(zs)
        for (i_p, p) in enumerate(ps)
            # the pole and zero are equal...
            # *and* if this pole has not already been canceled...
	    if isapprox(p, z) && (! canceled_poles[i_p])
                canceled_zeros[i_z] = true
                canceled_poles[i_p] = true
                break
            end
        end
    end
    @assert sum(canceled_zeros) == sum(canceled_poles)
    if verbose && sum(canceled_poles) != 0
        println("canceling the following poles and zeros: ", ps[canceled_poles])
    end
    
    if sum(canceled_zeros) == 0
        return tf
    else    
        # reconstruct
        return zeros_poles_k(zs[.! canceled_zeros], ps[.! canceled_poles], k, time_delay=tf.time_delay)
    end
end

"""
    o = system_order(tf::TransferFunction)

return the order of the numerator and denominator of the transfer function `tf`.

use [`pole_zero_cancellation`](@ref) first if you wish to cancel poles and zeros that are equal before determining the order.

# returns
`o::Tuple{Int, Int}`: (order of numerator, order of denominator)

# examples
```jldoctest
g = 1 / (s + 1)
system_order(g)
# output 
(0, 1)

g = (s + 1) / ((s + 2) * (s + 3))
system_order(g)
# output
(1, 2)
```
"""
function system_order(tf::TransferFunction)
    return (degree(tf.numerator), degree(tf.denominator))
end
