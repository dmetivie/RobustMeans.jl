abstract type AbstractSmoothingKernel end

"""
    SmoothingKernel

- `h` is the bandwidth parameter.
- `f` is a kernel function.

Exemple of kernels are (but you can define your own)

```julia
uniform(u::T) where T<:Real = abs(u) > one(T) ? zero(T) : 1 / 2
epanechnikov(u::T) where T<:Real = abs(u) > one(T) ? zero(T) : 3 / 4 * (1 - u^2)
gaussian(u::T) where T<:Real = abs(u) > one(T) ? zero(T) : (1 / sqrt(2 * pi)) * exp(-1 / 2 * u^2)
triangular(u::T) where T<:Real = abs(u) > one(T) ? zero(T) : (1 - abs(u))
biweight(u::T) where T<:Real = abs(u) > one(T) ? zero(T) : 15 / 16 * (1 - u^2)^2
triweight(u::T) where T<:Real = abs(u) > one(T) ? zero(T) : 35 / 32 * (1 - u^2)^3
tricube(u::T) where T<:Real = abs(u) > one(T) ? zero(T) : 70 / 81 * (1 - abs(u)^3)^3
cosine(u::T) where T<:Real = abs(u) > one(T) ? zero(T) : (pi / 4) * cos((pi / 2) * u)
logistic(u::T) where T<:Real = abs(u) > one(T) ? zero(T) : 1 / (exp(u) + 2 + exp(-u))
```

"""
struct SmoothingKernel{T<:Real}
    f::Function
    h::T
end

#NOTE: previously type unstable version lead to horrible performances! Hope this one is better

uniform(u::T) where {T<:Real} = abs(u) > one(T) ? zero(T) : 1 / 2

triangular(u::T) where {T<:Real} = abs(u) > one(T) ? zero(T) : (1 - abs(u))

gaussian(u::T) where {T<:Real} = (1 / sqrt(2 * pi)) * exp(-1 / 2 * u^2)

epanechnikov(u::T) where {T<:Real} = abs(u) > one(T) ? zero(T) : 3 / 4 * (1 - u^2)

biweight(u::T) where {T<:Real} = abs(u) > one(T) ? zero(T) : 15 / 16 * (1 - u^2)^2

triweight(u::T) where {T<:Real} = abs(u) > one(T) ? zero(T) : 35 / 32 * (1 - u^2)^3

tricube(u::T) where {T<:Real} = abs(u) > one(T) ? zero(T) : 70 / 81 * (1 - abs(u)^3)^3

cosine(u::T) where {T<:Real} = abs(u) > one(T) ? zero(T) : (pi / 4) * cos((pi / 2) * u)

logistic(u::T) where {T<:Real} = abs(u) > one(T) ? zero(T) : 1 / (exp(u) + 2 + exp(-u))