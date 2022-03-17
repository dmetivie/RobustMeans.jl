"""
bound(estimator, δ, n::Int)
Give the theoritical bound B such that for the estimator 
    #TODO! correction math formula
```math
\\mathbb{P}( |X- \\mathbb{E}(X)|/\\sigma \\leq B) \\leq \\delta
```

"""
bound(δ, n::Integer, M::EmpiricalMean) = sqrt(1 / (n * δ))
bound(δ, n::Integer, M::δMoM) = sqrt(32log(1 / δ) / n)
bound(δ, n::Integer, M::δTrimmedMean) = 9sqrt(2log(8 / δ) / n)
bound(δ, n::Integer, M::δCatoni) = sqrt(2log(2 / δ) / (n - 2log(2 / δ)))
bound(δ, n::Integer, M::δHuber) = 8sqrt(2log(4 / δ) / n)
""""
sqrt(2log(1 / δ) / n) # Multiply by a factor (1 + o(1)) where 
o(1) = (1+O(sqrt(log(1/δ)/n)))*(1+log(log(1/δ))/log(1/δ)) goes to zeros with (log(1/δ)/n, δ)→ (0, 0)
see https://www.youtube.com/watch?v=Kr0Kl_sXsJM Q&A
"""
bound(δ, n::Integer, M::δLeeValiant) = 8sqrt(2log(4 / δ) / n)


function bound(n::Integer,)
    if estimator == :EM
        return sqrt(1 / (n * δ))
    elseif estimator == :MM
        return sqrt(32log(1 / δ) / n)
    elseif estimator == :TM
        return 9sqrt(2log(8 / δ) / n)
    elseif estimator == :MN
        return sqrt(2log(1 / δ) / n) # not true bound since the cst is unknown
    elseif estimator == :LV
        return sqrt(2log(1 / δ) / n) # Multiply by a factor (1 + o(1)) where o(1) goes to zeros with (log(1/δ)/n, δ)→ (0, 0)
    elseif estimator == :CA
        return sqrt(2log(2 / δ) / (n - 2log(2 / δ)))
    elseif estimator == :HU
        return 8sqrt(2log(4 / δ) / n)#8sqrt(2 * (log(2 / δ) - log(1 - exp(-n / 8) / δ)) / n)
    else
        println("$(estimator) Not done")
    end
end
