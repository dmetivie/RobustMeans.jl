"""
    RobustMeans

Simple package implementing some Robust mean estimator.
"""
module RobustMeans

using Statistics
using StatsBase
using NonlinearSolve

import Statistics: mean

abstract type MeanEstimator end
abstract type RobustMean <: MeanEstimator end

include("means.jl")
include("mean_delta.jl")
# include("delta_mean.jl")
include("bounds.jl")

export MeanEstimator, RobustMean
export EmpiricalMean
export MedianOfMean, TrimmedMean
export Z_Estimator, Catoni, Huber, ψ_Huber, ψ_Catoni
export LeeValiant, MinskerNdaoud
export bound
export mean

include("smoothingkernel.jl")
include("rollingmean.jl")
include("means_weighted.jl")
include("robustnonparamregression.jl")
export SmoothingKernel,
    uniform,
    triangular,
    gaussian,
    epanechnikov,
    biweight,
    triweight,
    tricube,
    cosine,
    logistic,
    r̂ₙMOM, r̂ₙ

end

