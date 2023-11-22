module Theory

export extendParams!, getλFromTotalDivisions

expGrowthRateFromNT(Nf, t) = log(Nf)/t
expGrowthRateFromNT(Ni, Nf, t) = log(Nf/Ni)/t

function extendParams!(params::Dict)
    if !haskey(params, "ρ")
        params["ρ"] = params["λ"]*(1-params["p"])
        params["ϕ"] = params["λ"]*params["p"]
    end
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
    Np = params["pure births"]
    tP = Np > 0 ? log(Np)/γ : 0

    λ = ( r + log(4) - 2*log(Nf+1) ) / 
    ( (2-p)*(t-tP) - (1-p)*(t-tM)/Nf - (1-p)*( exp(-γ*tP)-1/Nf )/γ )

    return λ
end

function getSeλFromTotalMutations(params::Dict)
    # seM = params["ste mutations"]
    σM = params["σ mutations"]
    t = params["evolve time"]
    tM = params["mature time"]
    Nf = params["N final"]
    Ni = params["N initial"]
    p = params["p"]
    μ = params["μ"]
    γ = expGrowthRateFromNT(params["N final"], params["mature time"])
    Np = params["pure births"]
    S = params["sample size"]
    tP = Np > 0 ? log(Np)/γ : 0

    seλ = σM * 1/sqrt(S) * 1/μ /
    ( (2-p)*(t-tP) - (1-p)*(t-tM)/Nf - (1-p)*( exp(-γ*tP)-1/Nf )/γ )

    return seλ
end


end