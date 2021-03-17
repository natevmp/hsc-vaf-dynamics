module Theory

export extendParams!, getλFromTotalDivisions

expGrowthRateFromNT(Nf, t) = log(Nf)/t

function extendParams!(params::Dict)
    params["ρ"] = params["λ"]*(1-params["p"])
    params["ϕ"] = params["λ"]*params["p"]
    params["N"] = params["N final"]
    γ = expGrowthRateFromNT(params["N final"], params["mature time"])
    params["growth rate"] = γ
    return params
end

function extendParamsFromData!(params::Dict)
    params["λ"] = getλFromTotalDivisions(params)
    extendParams!(params)
end

function getEffectiveDivisions(params)
    t = params["evolve time"]
    tM = params["mature time"]
    γ = params["growth rate"]
    ρ = params["ρ"]
    ϕ = params["ϕ"]
    Nf = params["N final"]
    Ni = params["N initial"]
    λEff = ϕ*t + 2γ*tM + 2ρ*t - (ρ/γ)*(1/Ni - 1/Nf) - ρ*(t-tM)/Nf
    return λEff
end

function getCPRateParamsFromBurdenStats(mean, var; μKnown=nothing)
    if isnothing(μKnown)
        μ = var / mean - 1; r = mean^2/(var-mean)
        return r, μ
    else 
        μ = μKnown; r = mean/μ
        return r, μ
    end
end

function getλFromTotalDivisions(params::Dict)
    r = params["divisions"]
    t = params["evolve time"]
    tM = params["mature time"]
    Nf = params["N final"]
    Ni = params["N initial"]
    p = params["p"]
    γ = expGrowthRateFromNT(params["N final"], params["mature time"])

    λ = (r - 2log(Nf)) / 
        ( p*t + (1-p)*( (2-1/Nf)*t + tM/Nf - (Nf/Ni-1)/(Nf*γ) ) )
    return λ
end




end