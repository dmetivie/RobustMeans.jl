"""
bound(estimator, δ, n::Int)
Give the theoritical bound B such that for the estimator 
    #TODO! correction math formula
```math
\\mathbb{P}( |X- \\mathbb{E}(X)|/\\sigma \\leq B) \\leq \\delta
```

"""
bound(n, δ, M::EmpiricalMean) = sqrt(1 / (n * δ))
bound(n, M::EmpiricalMean) = sqrt(1 / (n * M.δ))

function bound(n, δ, M::δ_MeanEstimator)
    M.δ = δ
    return bound(n, M)
end
bound(n, M::δMoM) = sqrt(32log(1 / M.δ) / n)
bound(n, M::δTrimM) = 19sqrt(2log(8 / M.δ) / n) # Typo in the review by Lugosie. It seems that the ouragous 19 does th job. It is probably improvable.
bound(n, M::δCatoni) = sqrt(2log(2 / M.δ) / (n - 2log(2 / M.δ)))
bound(n, M::δHuber) = 8sqrt(2log(4 / M.δ) / n)
""""
sqrt(2log(1 / δ) / n) # Multiply by a factor (1 + o(1)) where 
o(1) = (1+O(sqrt(log(1/δ)/n)))*(1+log(log(1/δ))/log(1/δ)) goes to zeros with (log(1/δ)/n, δ)→ (0, 0)
see https://www.youtube.com/watch?v=Kr0Kl_sXsJM Q&A
"""
bound(n, M::δLeeValiant) = 8sqrt(2log(4 / M.δ) / n)
