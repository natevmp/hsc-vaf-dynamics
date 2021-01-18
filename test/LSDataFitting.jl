include("../src/vafdyn.jl")

using DelimitedFiles
using ApproxBayes
using Distributions, Distances
using .VAFDyn

using Plots
gr()

## ============================== Load Data ==============================

agePatient = 59

untypedM = readdlm("data/Shearwater_calls_FDR0.95_all_muts.txt", '\t', Any; skipstart=1)
untypedM = untypedM[:, 5: end-1]
replace!(untypedM, "NA"=>0)

variants_var_col = Array{Int}(untypedM)

## ===== Order data =====
HSCMask_col = fill(false, size(variants_var_col, 2))
HSCMask_col[1:73] .= true
# HSCMask_col[125:end] .= true
HPCMask_col = .!HSCMask_col
sampleSize = sum(HSCMask_col)

# mutational burden
mutBurden_col = sum(variants_var_col, dims=1)

mutHSCBurdenAv = mean(mutBurden_col[HSCMask_col])
mutHSCBurdenVar = var(mutBurden_col[HSCMask_col])


println(mutHSCBurdenAv)
println(mutHSCBurdenVar)
println(mutHSCBurdenVar/mutHSCBurdenAv - 1)
# println(length(mutBurden_col[HSCMask_col]))

## ===== variant allele frequencies =====
vaf_var = sum(variants_var_col, dims=2) / 140
prevHSC_var = sum(variants_var_col[:, HSCMask_col], dims=2)

function makeHist(varPrev_var, N::Integer)
    n_f = zeros(Int64, N+1)
    for m in 1:N
        n_f[1+m] = length(findall(x -> x==m, varPrev_var))
    end
    return n_f
end

vafHSCData_f = makeHist(prevHSC_var, sampleSize)
freqs_f = range(0, 1; length=sampleSize+1)

p0 = plot(freqs_f, vafHSCData_f, yaxis=:log10,
ylims=(10^-0.3, 10^5),
label="", linewidth=3)
xlabel!("variant frequency")
ylabel!("number of variants")
title!("sampled data VAF")
display(p0)


## ============================== ABC pipeline ==============================

# === normal fitting ===
# dataParams = Dict(
#     "evolveTime" => agePatient,
#     "sampleSize" => sum(HSCMask_col),
#     "mutMean" => mutHSCBurdenAv,
#     "mutVar" => mutHSCBurdenVar,
#     "μ" => mutHSCBurdenVar/mutHSCBurdenAv - 1,
#     "λPar" => mutMean^2 / ((mutVar - mutMean))
# )

# === use known mu=1.2 ===
dataParams = Dict(
    "evolveTime" => agePatient,
    "sampleSize" => sum(HSCMask_col),
    "mutMean" => mutHSCBurdenAv,
    "mutVar" => mutHSCBurdenVar,
    "μ" => 1.2,
    "λPar" => mutHSCBurdenAv / 1.2
)

simParams = Dict(
    "length"=>1001,
    "sMax"=>1
)

function simLV(paramsEst, constants, targetvaf_f)

    # ABC parameters to estimate
    N = Integer(round(paramsEst[1]))
    p = paramsEst[2]

    # parameters from data
    mutMean = dataParams["mutMean"]
    mutVar = dataParams["mutVar"]
    t = dataParams["evolveTime"]
    S = dataParams["sampleSize"]
    λ = dataParams["λPar"] / ((2-p)*t)
    l = simParams["length"]
    sMax = simParams["sMax"]

    # parameters for simulation
    partParams = Dict(
        "N" => N,
        "ρ" => λ*(1-p),
        "ϕ" => λ*p,
        "μ" => dataParams["μ"]
    )

    # PDE simulation
    vfs = VAFDyn.VFreqspace(N, l)
    VAFDyn.evolveVAFfd(vfs, partParams, t)
    dfs = VAFDyn.makeDFSfromVFS(vfs, N)
    sDfs = VAFDyn.sampler(dfs, S, sMax=sMax)
    # compare with data VAF
    # d = sum( (log10.(abs.(sDfs.n_f[2:sMax])) .- log10.(targetvaf_f[2:sMax])).^2 )
    # d = sum( (abs.(sDfs.n_f[2:sMax]) .- targetvaf_f[2:sMax]).^2 )
    d = abs(sDfs.n_f[2] - targetvaf_f[2])
    # println(d)
    return d, 1
end

#define ABC setup type
setup = ABCRejection(
    simLV,
    2,
    500.0,
    Prior([DiscreteUniform(10^4, 4*10^5), Uniform(0.01, 0.99)]),
    maxiterations = 1*10^4,
    nparticles = 40
    )


# ===== Run ABC =====
println("running ABC")
#run ABC SMC algorithm
@time simResult = runabc(setup, vafHSCData_f, verbose=true, progress=true, parallel=true)
println("finished ABC")
#show results
show(simResult)


# ================ Basic plotting ================

param1Accepted_part = []
param2Accepted_part = []

for particle in simResult.particles
    push!(param1Accepted_part, particle.params[1])
    push!(param2Accepted_part, particle.params[2])
end

h1 = histogram(param1Accepted_part, label="", color=1, bins=20)
xlabel!("N")
ylabel!("accepted values")
h2 = histogram(param2Accepted_part, label="", color=2, bins=20)
xlabel!("p")
# ylabel!("accepted values")
p1 = plot(h1, h2, layout=2)

p2 = scatter(param2Accepted_part, param1Accepted_part, label="accepted particle")
xlabel!("p")
ylabel!("N")

p3 = plot(p1, p2, layout=(2,1), size=(600, 800))
display(p3)

savefig(p3, "ABC_N-p_Sampled.pdf")

## ========== Saving data ==========

using JLD2
# saveName = "./test/ABC_VAF1fit_simResult.jld2"
# @save saveName simResult

saveName = "./test/ABC_VAF1fit_mu"*string(round(dataParams["μ"],digits=2))*".jld2"
@save saveName param1Accepted_part param2Accepted_part dataParams

# using JLD

# save("./test/ABC_VAF1fit.jld", "param1Accepted_part", param1Accepted_part, "param2Accepted_part", param2Accepted_part, "dataParams", dataParams)


# save("./ABC_VAF1fit_simResult.jld", "simResultJLD", simResult)

