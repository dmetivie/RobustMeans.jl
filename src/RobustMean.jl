"""
    RobustMean

Simple package implementing some Robust mean estimator.
"""
module RobustMean

using Statistics

import Statistics: mean
# Write your package code here.


abstract type MeanEstimator end

"""
    mean(A::AbstractArray, MeanEstimator::EmpiricalMean)
    
The usual empirical mean estimator. Nothing new.
"""
struct EmpiricalMean <: MeanEstimator end
function mean(A::AbstractArray, MeanEstimator::EmpiricalMean)
    return mean(A)
end

"""
mean(A::AbstractArray, MeanEstimator::MoM)

The Median of Mean estimator.
"""
struct MoM{T<:Integer} <: MeanEstimator
    k::T # number of blocks
end
function mean(A::AbstractArray, MeanEstimator::MoM)
    return MedianOfMean(A, MeanEstimator.k)
end

"""
mean(A::AbstractArray, MeanEstimator::Z_Estimator)

A Z estimator given an influence function `x->ψ(x)` and a scaling parameter `α`.
"""
struct Z_Estimator{F,T} <: MeanEstimator
    α::T # scaling parameter
    ψ::F # Influence function
end
function mean(A::AbstractArray, MeanEstimator::Z_Estimator; kwargs...)
    return Z_estimator(A, MeanEstimator.α, MeanEstimator.ψ; kwargs...)
end
include("bounds.jl")
end

