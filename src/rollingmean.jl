import StatsBase: mean

mean(y::AbstractVector, kernel::SmoothingKernel, args...) = mean(y, 1:length(y), kernel, args...)

function mean(y::AbstractVector, x::AbstractVector, kernel::SmoothingKernel, args...)
    N = length(y)
    N == length(x) || throw(DimensionMismatch("x and y must have same length"))
    w = similar(y)
    m = similar(y)

    for (i, xx) in enumerate(x)
        @views w[:] .= kernel.f.((x .- xx) / kernel.h)
        m[i] = mean(y, weights(w), args...)
    end
    return m
end