#!/usr/bin/env julia

module ABCVaf

include("../src/vafdyn.jl")
using .VAFDyn
using GpABC, OrdinaryDiffEq, Distances, Distributions, JLD
using Plots
plotly()

function main()

# measurement:
paramsTrue = Dict(
    "N"=>450,
    "ρ"=>20.0,
    "ϕ"=>6.0,
    "μ"=>1.2
)
sampleSize = 80
evolveTime = 59
dfsTrue = VAFDyn.DFreqspace(paramsTrue["N"])
VAFDyn.evolveVAF(dfsTrue, paramsTrue, evolveTime)
sampFsTrue = VAFDyn.sampler(dfsTrue, sampleSize)

referenceData = zeros(Float64, 1, length(sampFsTrue.n_f)-2)
referenceData[1, :] = sampFsTrue.n_f[2:end-1]
# referenceData[1, :] = log.(sampFsTrue.n_f[2:end-1])

# h = plot(referenceData[1,:])
# display(h)

# ABC settings
# est params: N, ρ, ϕ
# priors = [Uniform(100, 1000), Uniform(1.0, 50.0), Uniform(1.0, 50.0)]
priors = [Uniform(100, 1000), Uniform(1.0, 50.0)]
# priors = [Uniform(100, 1000)]

# A function that simulates the model
function simulatorFunction(paramsEst)
    params = Dict(
        "N"=> Integer(round(paramsEst[1])),
        # "ρ"=>paramsTrue["ρ"],
        "ρ"=>paramsEst[2],
        "ϕ"=>paramsTrue["ϕ"],
        # "ϕ"=>paramsEst[3],
        "μ"=>paramsTrue["μ"]
    )
    dfs = VAFDyn.DFreqspace(params["N"])
    VAFDyn.evolveVAF(dfs, params, evolveTime)
    sampFs = VAFDyn.sampler(dfs, sampleSize)

    # ABC algorithm really wants a 2d matrix for some reason...
    simData = zeros(Float64, 1, length(sampFs.n_f)-2)
    simData[1, :] = sampFs.n_f[2:end-1]
    # simData[1, :] = log.(sampFs.n_f[2:end-1])
    return simData
end

# Simulation
nParticles = 200
threshold = 7.0
maxIter = 1e4

println("running ABC")

@time simResult = SimulatedABCRejection(referenceData, simulatorFunction, priors, threshold, nParticles; max_iter=convert(Int, maxIter), write_progress=false)

# nDesignPoints = 200
# @time emuResult = EmulatedABCRejection(referenceData, simulatorFunction, priors, threshold, nParticles, nDesignPoints; max_iter=convert(Int, 5e3), write_progress=false)
save("ABCrun3d_"*string(Integer(maxIter))*"iter.jld", "paramsTrue", paramsTrue, "sampleSize", sampleSize, "evolveTime", evolveTime, "simResult", simResult)
# save("ABCrun3d.jld", "paramsTrue", paramsTrue, "sampleSize", sampleSize, "evolveTime", evolveTime, "emuResult", emuResult)
p1 = plot(simResult)
# p1 = plot(emuResult)
display(p1)

end




main()

end
