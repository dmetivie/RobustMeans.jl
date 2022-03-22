"""
    RobustMean

Simple package implementing some Robust mean estimator.
"""
module RobustMean

using Statistics
using StatsBase: harmmean
using Roots

import Statistics: mean

# Write your package code here.


abstract type MeanEstimator end

include("means.jl")
include("delta_mean.jl")
include("bounds.jl")

export EmpiricalMean, MoM, TrimM, Z_Estimator
export δMoM, δTrimM, δCatoni, δHuber, δLeeValiant, δMinskerNdaoud
export bound
end

