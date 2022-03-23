"""
bound(estimator, δ, n::Int)
Give the theoritical bound B such that for the estimator 
    #TODO! correction math formula
```math
\\mathbb{P}( |X- \\mathbb{E}(X)|/\\sigma \\leq B) \\leq \\delta
```

"""
bound(δ::Real, n, M::EmpiricalMean) = sqrt(1 / (n * δ))
bound(δ::Real, n, M::MedianOfMean) = sqrt(32log(1 / δ) / n)
bound(δ::Real, n, M::TrimmedMean) = 19sqrt(2log(8 / δ) / n) # Typo in the review by Lugosie. It seems that the ouragous 19 does th job. It is probably improvable.
bound(δ::Real, n, M::Catoni) = sqrt(2log(2 / δ) / (n - 2log(2 / δ)))
bound(δ::Real, n, M::Huber) = 8sqrt(2log(4 / δ) / n)
""""
sqrt(2log(1 / δ) / n) # Multiply by a factor (1 + o(1)) where 
o(1) = (1+O(sqrt(log(1/δ)/n)))*(1+log(log(1/δ))/log(1/δ)) goes to zeros with (log(1/δ)/n, δ)→ (0, 0)
see https://www.youtube.com/watch?v=Kr0Kl_sXsJM Q&A
"""
bound(δ::Real, n, M::LeeValiant) = 8sqrt(2log(4 / δ) / n)
