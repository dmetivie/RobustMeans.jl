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
include("mean_delta.jl")
# include("delta_mean.jl")
include("bounds.jl")

export EmpiricalMean
export MedianOfMean, TrimmedMean
export Z_Estimator, Catoni, Huber
export LeeValiant, MinskerNdaoud
export bound
end

