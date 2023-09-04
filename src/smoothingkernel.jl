abstract type AbstractSmoothingKernel end

"""
    SmoothingKernel

- `h` is the bandwidth parameter.
- `f` is a kernel function.

Exemple of kernels are (but you can define your own)

```julia
uniform(u::Real) = abs(u) > 1 ? 0 : 1 / 2
epanechnikov(u::Real) = abs(u) > 1 ? 0 : 3 / 4 * (1 - u^2)
gaussian(u::Real) = abs(u) > 1 ? 0 : (1 / sqrt(2 * pi)) * exp(-1 / 2 * u^2)
triangular(u::Real) = abs(u) > 1 ? 0 : (1 - abs(u))
biweight(u::Real) = abs(u) > 1 ? 0 : 15 / 16 * (1 - u^2)^2
triweight(u::Real) = abs(u) > 1 ? 0 : 35 / 32 * (1 - u^2)^3
tricube(u::Real) = abs(u) > 1 ? 0 : 70 / 81 * (1 - abs(u)^3)^3
cosine(u::Real) = abs(u) > 1 ? 0 : (pi / 4) * cos((pi / 2) * u)
logistic(u::Real) = abs(u) > 1 ? 0 : 1 / (exp(u) + 2 + exp(-u))
```

```math
Y(x) = \\sum_{x_i \\in X} \\left( y_i f((x-x_i)/h) \\right) / \\sum_{x_i \in X} \\left( f((x-x_i)/h) \\right)
```

"""
struct SmoothingKernel
    f::Function
    h::Real
end

uniform(u::Real) = abs(u) > 1 ? 0 : 1 / 2

triangular(u::Real) = abs(u) > 1 ? 0 : (1 - abs(u))

gaussian(u::Real) = abs(u) > 1 ? 0 : (1 / sqrt(2 * pi)) * exp(-1 / 2 * u^2)

epanechnikov(u::Real) = abs(u) > 1 ? 0 : 3 / 4 * (1 - u^2)

biweight(u::Real) = abs(u) > 1 ? 0 : 15 / 16 * (1 - u^2)^2

triweight(u::Real) = abs(u) > 1 ? 0 : 35 / 32 * (1 - u^2)^3

tricube(u::Real) = abs(u) > 1 ? 0 : 70 / 81 * (1 - abs(u)^3)^3

cosine(u::Real) = abs(u) > 1 ? 0 : (pi / 4) * cos((pi / 2) * u)

logistic(u::Real) = abs(u) > 1 ? 0 : 1 / (exp(u) + 2 + exp(-u))