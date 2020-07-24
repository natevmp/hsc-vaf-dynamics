module BurdenDyn

using Random, Distributions, DifferentialEquations
using Plots
gr()


"""
Draw samples from the compound Poisson distribution with Poisson distributed elements.
λ           Poisson parameter for distribution of sum length
μ           Poisson parameter for (Poisson) distribution of sum elements
nSamples    number of samples to generate
returns array of generated samples
"""

function evolveBurden(params::Dict, t, mMax::Integer, ϵ::Real)

    μ = params["μ"]
    λ = params["λ"]
    p = params["p"]
    N = params["N"]
    ρ = λ*p
    ϕ = λ*(1-p)

    Pμ = Poisson(μ)
    #start threshold at expected value
    mutTh = Int(round(μ))
    probMutTh = pdf(Pμ, mutTh)
    while probMutTh > ϵ
        mutTh += 1
        probMutTh = pdf(Pμ, mutTh)
    end
    # number of cells per mutational burden m
    n0_m = zeros(Float64, mMax+1)
    n0_m[1] = N
    probMuts_m = pdf.(Pμ, 0:mutTh)
    println(probMutTh)
    println(mutTh)
    println(probMuts_m)

    function step!(dn_m, n_m, p, t)
        dn_m .= -(2ρ+ϕ)*n_m
        for m in 1:(mMax+1-mutTh)
            # dn_m[m:m+mutTh] += λ/N * 2*n_m[m] * probMuts_m
            for k in 0:mutTh
                dn_m[m+k] += (2ρ+ϕ) * n_m[m] * probMuts_m[1+k]
            end
        end
        # for i in 1:mutTh
        #     dn_m[mMax+1-mutTh + i] = -λ/N * 2*n_m[mMax+1-mutTh+i]/N
        #     dn_m[mMax+1-mutTh+i:end] += λ/N * (2*n_m[mMax+1-mutTh+i:end] .* probMuts_m[1:end-i])
        # end
    end

    prob = ODEProblem(step!, n0_m, (0.0, t))
    alg = TRBDF2() #stablest for stiff PDE
    sol = solve(prob, alg, save_everystep=false)

    nCells_m = sol.u[2]

    return nCells_m
end

end
