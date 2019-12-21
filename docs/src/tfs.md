# transfer functions

the response [output $Y(s)$] of a linear, time-invariant system to any input [$U(s)$] is characterized by a transfer function $g(s)=Y(s)/U(s)$.

## constructing a transfer function

a data structure, `TransferFunction`, is used to represent transfer functions.

for example, we can construct the transfer function:
$$g(s)=\dfrac{5s+1}{s^2 + 4s+1}$$
in an intuitive way that resembles the algebraic expression:

```julia
g = (5 * s + 5) / (s ^ 2 + 4 * s + 5)
```

alternatively, we can pass the powers of the coefficients of the $s$ variables of $g(s)$ in the numerator and denominator, respectively. The coefficients of the highest powers go first.
```julia
g = TransferFunction([5, 5], [1, 4, 5]) # way 2
```

note that `s == TransferFunction([1, 0], [1])`.

we can also represent $g(s)$ by specifying its zeros, poles, and gain:

```julia
g = zeros_poles_gain([-1.0], [-2 + im, -2 - im], 1.0)  # way 3
```

the `TransferFunction` data structure has a `numerator` (a polynomial), `denominator` (a polynomial), and `time_delay` (a number) attribute that we can access.

```julia
g = (5 * s + 5) / (s ^ 2 + 4 * s + 5)
g.numerator # 5s + 1, a `Poly`
g.denominator # s² + 4s + 5, a `Poly`
g.time_delay # 0.0, a `Float64`
```

The numerator and denominator are `Poly` types from the package [Polynomials.jl](https://github.com/JuliaMath/Polynomials.jl).

## time delays

add a time delay to a transfer function as follows. 

```julia
θ = 2.0 # time delay
g = 3 / (2 * s + 1) * exp(-θ * s) # way 1
g = TransferFunction([3], [2, 1], θ) # way 2
```

the transfer function `g` represents:
$$g(s)=\dfrac{3}{2s+1}e^{-2s}$$.

## transfer function algebra

we can add `+`, subject `-`, mutiply `*`, and divide `/` transfer functions.

```julia
g1 = 3 / (s + 2)
g2 = 1 / (s + 4)
g = g1 * g2
```

## evaluate a transfer function at a complex number

for example, to evaluate $g(s)=\dfrac{4}{s+2}$ at $s=1-i$:
```julia
g = 4 / (s + 2)
evaluate(g, 2 * im) # 1 - im
```
see the Julia documentation on imaginary numbers [here](https://docs.julialang.org/en/v1/manual/complex-and-rational-numbers/).

## poles, zeros and gain of a transfer function

compute the poles, zeros, and gain of a transfer function as follows:

```julia
g = (5 * s + 5) / (s ^ 2 + 4 * s + 5)
z, p, k = zeros_poles_gain(g)
# z = [-1.0]
# p = [-2-im, -2+im]
# k  = 1.0
```

```@docs
    TransferFunction
    zeros_poles_gain
    evaluate
    proper
    strictly_proper
```
