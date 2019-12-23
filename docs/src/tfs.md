# transfer functions

the response [output $Y(s)$] of a linear, time-invariant system to any input [$U(s)$] is characterized by a transfer function $g(s)=Y(s)/U(s)$.

## constructing a transfer function

a data structure, `TransferFunction`, represents a transfer function.

for example, to construct the transfer function
$$g(s)=\dfrac{5s+1}{s^2 + 4s+5}$$
in an intuitive way that resembles the algebraic expression:

```julia
julia> g = (5 * s + 1) / (s ^ 2 + 4 * s + 5)
```

alternatively, construct a `TransferFunction` using the powers of the coefficients of the $s$ variables of $g(s)$ in the numerator and denominator polynomials, respectively. The coefficients of the highest powers go first.
```julia
julia> g = TransferFunction([5, 1], [1, 4, 5]) # way 2
```

note that `s == TransferFunction([1, 0], [1])`.

the `TransferFunction` data structure has a `numerator` (a polynomial), `denominator` (a polynomial), and `time_delay` (a number) attribute that we can access.

```julia
julia> g.numerator # 5s + 1, a `Poly`
julia> g.denominator # s² + 4s + 5, a `Poly`
julia> g.time_delay # 0.0, a `Float64`
```

The numerator and denominator (polynomials) are `Poly` types from the package [Polynomials.jl](https://github.com/JuliaMath/Polynomials.jl).

## time delays

add a time delay to a transfer function as follows. 

```julia
julia> θ = 2.0 # time delay
julia> g = 3 / (2 * s + 1) * exp(-θ * s) # way 1
julia> g = TransferFunction([3], [2, 1], θ) # way 2
```

the transfer function `g` represents:
$$g(s)=\dfrac{3}{2s+1}e^{-2s}$$.

## zeros, poles, k-factor representation

we can write any transfer function $g(s)$ in terms of its poles ($p_j$), zeros ($z_j$), and k-factor ($k$):

$$g(s)=k\dfrac{\Pi_j (s-z_j)}{\Pi_j(s-p_j)}e^{-\theta s}$$

the scalar factor $k$ allows us to uniquely specify a transfer function in terms of its poles, zeros, and time delay.

for example:

$$g(s)=\dfrac{5s+1}{s^2 + 4s+5}=5\dfrac{(s+1/5)}{(s-2+i)(s-2-i)}$$

#### construting a transfer function from its zeros, poles and k-factor

```julia
julia> g = zeros_poles_k([-1/5], [-2 + im, -2 - im], 5.0, time_delay=0.0)  # way 3
```

#### computing the poles, zeros, and k-factor from a transfer function

```
julia> g = (5 * s + 5) / (s ^ 2 + 4 * s + 5)
julia> zeros_poles_k(g) # [-1/5], [-2+im, 2+im], 5
```

## transfer function algebra

we can add `+`, subject `-`, multiply `*`, and divide `/` transfer functions.

```julia
julia> g1 = 3 / (s + 2)
julia> g2 = 1 / (s + 4)
julia> g = g1 * g2 # 3 / (s^2 + 6s + 8)
```

## evaluate a transfer function at a complex number

for example, to evaluate $g(s)=\dfrac{4}{s+2}$ at $s=1-i$:
```julia
julia> g = 4 / (s + 2)
julia> evaluate(g, 2 * im) # 1 - im
```
see the Julia documentation on imaginary numbers [here](https://docs.julialang.org/en/v1/manual/complex-and-rational-numbers/).

## zero-frequency gain of a transfer function

compute the zero-frequency gain of a transfer function $g(s)$, which is $g(s)$ evaluated at $s=0$, as follows:

```julia
julia> g = (5 * s + 1) / (s ^ 2 + 4 * s + 5)
julia> zero_frequency_gain(g) # 0.2 = 1/5
```

## poles, zeros, and zero-frequency gain of a transfer function

compute the poles, zeros, and zero-frequency gain of a transfer function as follows:

```julia
julia> g = (5 * s + 5) / (s ^ 2 + 4 * s + 5)
julia> z, p, gain = zeros_poles_gain(g)
# z = [-1.0]
# p = [-2-im, -2+im]
# gain  = 1.0
```

## cancel poles and zeros

we can cancel pairs of identical poles and zeros in a transfer function.

```
julia> g = s * (s+1) / ((s+3) * s * (s+1) ^ 2)
julia> pole_zero_cancellation(g) # 1 / ((s+3) * (s+1))
```
## detailed docs

```@docs
    TransferFunction
    zero_frequency_gain
    zeros_poles_gain
    zeros_poles_k
    pole_zero_cancellation
    evaluate
    proper
    strictly_proper
```
