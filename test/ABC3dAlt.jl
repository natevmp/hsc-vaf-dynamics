include("../src/vafdyn.jl")
using .VAFDyn
using OrdinaryDiffEq, Distances, Distributions, JLD
using ApproxBayes


# create reference measurement:
paramsTrue = Dict(
    "N"=>450,
    "ρ"=>2.0,
    "ϕ"=>6.0,
    "μ"=>1.2
)
sampleSize = 80
evolveTime = 59
dfsTrue = VAFDyn.DFreqspace(paramsTrue["N"])
VAFDyn.evolveVAF(dfsTrue, paramsTrue, evolveTime)
sampFsTrue = VAFDyn.sampler(dfsTrue, sampleSize)

reference_f = sampFsTrue.n_f[2:end-1]

#simulations function for ABC. return distance (sum of squared distances) and solution
function simLV(paramsEst, constants, targetdata)
    params = Dict(
        "N"=> Integer(round(paramsEst[1])),
        "ρ"=>paramsEst[2],
        # "ϕ"=>paramsTrue["ϕ"],
        "ϕ"=>paramsEst[3],
        "μ"=>paramsTrue["μ"]
    )
    dfs = VAFDyn.DFreqspace(params["N"])
    VAFDyn.evolveVAF(dfs, params, evolveTime)
    sampFs = VAFDyn.sampler(dfs, sampleSize)

    d = sum(((@view sampFs.n_f[2:end-1]) .- reference_f).^2 )
    return d, 1
end

#define ABC setup type
setup = ABCRejection(
    simLV,
    3,
    500.0,
    Prior([Uniform(100.0, 1000.0), Uniform(0.0, 10.0), Uniform(0.0, 10.0)]),
    # Prior([Uniform(100.0, 1000.0), Uniform(0.5, 10.0)]),
    maxiterations = 3*10^4,
    nparticles = 50
    )

# setup = ABCSMC(
#     simLV,
#     2,
#     1.0,
#     Prior([Uniform(100.0, 1000.0), Uniform(1.0, 10.0)]),
#     maxiterations = 10^2,
#     convergence = 0.001,
#     nparticles = 400
#     )
println("running ABC")
#run ABC SMC algorithm
@time resSim = runabc(setup, reference_f, verbose = true, progress = true, parallel=false)

println("finished ABC")
#show results
show(resSim)

# save("ABCrun3d_new.jld", "paramsTrue", paramsTrue, "sampleSize", sampleSize, "evolveTime", evolveTime, "simResult", resSim)
save("ABCrun3d_new.jld", "paramsTrue", paramsTrue, "sampleSize", sampleSize, "evolveTime", evolveTime)
