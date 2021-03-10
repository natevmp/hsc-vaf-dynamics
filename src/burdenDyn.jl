module BurdenDyn

using Random, Distributions, DifferentialEquations
export evolveBurden, evolveBurdenGrowth, logisticGrowth, exponentialCappedGrowth, linearGrowth


""" ODE evolve the mutational Burden
params          parameter values
evolveTime      time to evolve system
nMax            maximum number of mutations to allow
ϵ               threshold for minimal probability to account for number of new mutations per division
"""
function evolveBurden(params::Dict, evolveTime, mMax::Integer, ϵ::Real)

    μ = params["μ"]
    λ = params["λ"]
    p = params["p"]
    N = params["N"]
    ρ = λ*(1-p)
    ϕ = λ*p

    probMuts_m, mutTh = mutationProbs(μ, ϵ)

    function step!(dn_m, n_m, p, t)
        dn_m .= -(2ρ+ϕ)*n_m
        for m in 1:(mMax+1-mutTh)
            for k in 0:mutTh
                dn_m[m+k] += (2ρ+ϕ) * n_m[m] * probMuts_m[1+k]
            end
        end
    end

    n0_m = zeros(Float64, mMax+1)
    n0_m[1] = N
    prob = ODEProblem(step!, n0_m, (0.0, evolveTime))
    alg = TRBDF2() #stablest for stiff PDE
    sol = solve(prob, alg, save_everystep=false)

    nCells_m = sol.u[2]

    return nCells_m
end

""" ODE evolve the mutational Burden in a growing population
params          parameter values
evolveTime      time to evolve system
nMax            maximum number of mutations to allow per cell
ϵ               threshold for minimal probability to account for number of new mutations per division
"""
function evolveBurdenGrowth(params::Dict, evolveTime; 
    mMax::Union{Integer,Nothing}=nothing, 
    ϵ::Real=0.0001,
    growthType::String="cappedExponential")

    λ = params["λ"]
    p = params["p"]
    μ = params["μ"]
    Ni = params["N initial"]
    K = params["N final"]
    r = params["growth rate"]
    ρ = λ*(1-p)
    ϕ = λ*p

    # γ(n) = r * ( 1 - n/K )   # logistic growth rate
    # Nt(t) = K / ( 1 + (K-Ni)/Ni * exp(-r*t) )   # logistic pop size
    # γ(n) = r   # exponential growth rate
    # Nt(t) = Ni*exp(r*t)   # exponential pop size
    γ(n) = n<K ? r : 0   # exponential capped growth rate
    Nt(t) = Ni*exp(r*t)<K ? Ni*exp(r*t) : K     # exponential capped pop size
    # γ(n) = Ni*r/n   # linear growth rate
    # Nt(t) = Ni * (1 + r*t)   # linear pop size

    # γ(n) = n<K ? params["exp growth rate"] : K*params["lin growth rate"]/n
    # Nt(t) = exponentialToLinearGrowth(Ni, K, params["exp growth rate"], params["lin growth rate"], t)

    # function setGrowthType()
    #     if growthType=="cappedExponential"
    #         println("growth type: capped exp")
    #         γ(n) = n<K ? r : 0   # exponential capped growth rate
    #         Nt(t) = Ni*exp(r*t)<K ? Ni*exp(r*t) : K     # exponential capped pop size
    #         return γ, Nt
    #     elseif growthType=="exponentialToLinear"
    #         println("growth type: exp to lin")
    #         γ(n) = n<K ? params["exp growth rate"] : K*params["lin growth rate"]/n
    #         Nt(t) = exponentialToLinearGrowth(Ni, K, params["exp growth rate"], params["lin growth rate"], t)
    #         return γ, Nt
    #     else
    #         println("error: growth type not set")
    #     end
    # end
    # γ, Nt = setGrowthType()

    if isnothing(mMax)
        mMax = Int(round(2*λ*evolveTime*μ*3))
    end
    
    probMuts_m, mutTh = mutationProbs(μ, ϵ)
    
    function step!(dn_m, n_m, p, t)
        dn_m .= -( ρ*(2-1/Nt(t))+ϕ+γ(Nt(t)) ) * n_m
        for m in 1:(mMax+1-mutTh)
            for k in 0:mutTh
                dn_m[m+k] += ( ρ*(2-1/Nt(t))+ϕ+2*γ(Nt(t)) ) * n_m[m] * probMuts_m[1+k]
            end
        end
    end
    
    n0_m = zeros(Float64, mMax+1)
    n0_m[1] = Ni
    prob = ODEProblem(step!, n0_m, (0.0, evolveTime))
    alg = TRBDF2() #stablest for stiff PDE
    sol = solve(prob, alg, save_everystep=false)

    nCells_m = sol.u[2]

    return nCells_m
end


""" Return vector of probabilies of number of mutations occurring in a daughter cell
μ       Mutations rate
ϵ       minimal probability threshold
"""
function mutationProbs(μ, ϵ)
    Pμ = Poisson(μ)
    mutTh = Int(round(μ))   #start threshold at expected value
    probMutTh = pdf(Pμ, mutTh)  # probability of current threshold
    while probMutTh > ϵ
        mutTh += 1
        probMutTh = pdf(Pμ, mutTh)  # update probability of current threshold
    end
    # vector of probabilities up until threshold
    probMuts_m = pdf.(Pμ, 0:mutTh)
    # println(probMutTh)
    # println(mutTh)
    # println(probMuts_m)
    return probMuts_m, mutTh
end

function logisticGrowth(Ni, K, r, t)
    return K / ( 1 + (K-Ni)/Ni * exp(-r*t) )
end

exponentialGrowth(Ni, r, t) = Ni*exp(r*t)
linearGrowth(Ni, r, t) = Ni + Ni*r*t

function exponentialCappedGrowth(Ni, K, r, t)
    return Ni*exp(r*t)<K ? Ni*exp(r*t) : K
end

function exponentialToLinearGrowth(Ni, K, rExp, rLin, t)
    tK = log(K/Ni)/rExp
    t<tK ? n=exponentialGrowth(Ni, rExp, t) : n=linearGrowth(K, rLin, t-tK)
    return n
end


end