module BurdenDyn

using Random, Distributions, DifferentialEquations
export evolveBurden, evolveBurdenGrowth, logisticGrowth, exponentialCappedGrowth


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
nMax            maximum number of mutations to allow
ϵ               threshold for minimal probability to account for number of new mutations per division
"""
function evolveBurdenGrowth(params::Dict, evolveTime, mMax::Integer, ϵ::Real)

    λ = params["λ"]
    p = params["p"]
    μ = params["μ"]
    Ni = params["N init"]
    K = params["N max"]
    r = params["growth rate"]
    ρ = λ*(1-p)
    ϕ = λ*p

    γ(n, t) = r * ( 1 - n/K )   # rate of Poisson growth divisions
    Nt(t) = K / ( 1 + (K-Ni)/Ni * exp(-r*t) )   # size of total population
    # γ(n, t) = r   # exponential growth rate
    # Nt(t) = Ni*exp(r*t)   # exponential pop size
    # γ(n, t) = n<K ? r : 0   # exponential capped growth rate
    # Nt(t) = Ni*exp(r*t)<K ? Ni*exp(r*t) : K     # exponential capped pop size

    probMuts_m, mutTh = mutationProbs(μ, ϵ)
    
    function step!(dn_m, n_m, p, t)
        dn_m .= -( 2ρ+ϕ+γ(Nt(t),t) ) * n_m
        for m in 1:(mMax+1-mutTh)
            for k in 0:mutTh
                dn_m[m+k] += ( 2ρ+ϕ+2*γ(Nt(t),t) ) * n_m[m] * probMuts_m[1+k]
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
    println(probMutTh)
    println(mutTh)
    println(probMuts_m)
    
    return probMuts_m, mutTh
end

function logisticGrowth(Ni, K, r, t)
    return K / ( 1 + (K-Ni)/Ni * exp(-r*t) )
end

function exponentialCappedGrowth(Ni, K, r, t)
    return Ni*exp(r*t)<K ? Ni*exp(r*t) : K
end



end