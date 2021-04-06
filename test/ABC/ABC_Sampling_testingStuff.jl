# %% codecell

include("src/vafdyn.jl")
using .VAFDyn
using OrdinaryDiffEq, Distances, Distributions
using JLD2
using ApproxBayes

# %% codecell
# create reference measurement:
paramsTrue = Dict(
    "N"=>11000,
    "ρ"=>1.0,
    "ϕ"=>2.0,
    "μ"=>5.8
)
evolveTime = 59
sampleSize = 79

function mutBurdenStats(params, evolveTime)
    ρ = params["ρ"]
    ϕ = params["ϕ"]
    λ = ρ + ϕ
    p = 1/(ρ/ϕ+1)
    t = evolveTime
    μ = params["μ"]
    mean = λ*(2-p)*t*μ
    var = λ*(2-p)*t*(μ+μ^2)
    return mean, var
end

display(paramsTrue)
println(mutBurdenStats(paramsTrue, evolveTime))

# %% codecell
mutMean, mutVar = mutBurdenStats(paramsTrue, evolveTime)


lSize = 80
# dfsTrue = VAFDyn.DFreqspace(paramsTrue["N"])
# @time VAFDyn.evolveVAF(dfsTrue, paramsTrue, evolveTime)
# sampFsTrue = VAFDyn.sampler(dfsTrue, sampleSize)

vfsTrue = VAFDyn.VFreqspace(paramsTrue["N"], lSize)
VAFDyn.evolveVAFfd(vfsTrue, paramsTrue, evolveTime)

reference_f = log10.(vfsTrue.n_f[2:end-1])

p1 = plot(dfsTrue.freqs_f[2:end-1], log10.(dfsTrue.n_f[2:end-1]*paramsTrue["N"]))
plot!(vfsTrue.freqs_f[2:end-1], reference_f, linestyle=:dash)
display(p1)

println(1/(paramsTrue["ρ"]/paramsTrue["ϕ"]+1))

# %% codecell
#simulations function for ABC. return distance (sum of squared distances) and solution
function simLV(paramsEst, constants, targetdata)
    p = paramsEst[2]
    t = evolveTime
    λ = mutMean^2 / ((mutVar - mutMean)*(2-p)*t)
    params = Dict(
        # "N"=> paramsTrue["N"],
        "N" => paramsEst[1],
        "ρ" => λ*(1-p),
        "ϕ" => λ*p,
        "μ" => (mutVar-mutMean)/mutMean
    )
    vfs = VAFDyn.VFreqspace(Integer(round(params["N"])), lSize)
    VAFDyn.evolveVAFfd(vfs, params, evolveTime)
    # sampFs = VAFDyn.sampler(dfs, sampleSize)

    # d = sum(((@view vfs.n_f[2:end-1]) .- reference_f).^2 )
    d = sum(((log10.(vfs.n_f[2:end-1])) .- reference_f).^2 )
    return d, 1
end

#define ABC setup type
setup = ABCRejection(
    simLV,
    2,
    1.0,
    Prior([DiscreteUniform(100, 1000), Uniform(0.001, 0.999)]),
    maxiterations = 1*10^5,
    nparticles = 500
    )

println("running ABC")
#run ABC SMC algorithm
@time simResult = runabc(setup, reference_f, verbose=true, progress=true, parallel=false)

println("finished ABC")
#show results
show(simResult)

# %% codecell

# save("ABC-MJWrun_test1.jld", "paramsTrue", paramsTrue, "evolveTime", evolveTime, "simResult", simResult)
@save "ABC-MJWrun_infer1.jld2" paramsTrue evolveTime simResult

# %% codecell
using Plots
gr()

# %% codecell

param1Accepted_part = []
param2Accepted_part = []

for particle in simResult.particles
    push!(param1Accepted_part, particle.params[1])
    push!(param2Accepted_part, particle.params[2])
end

h1 = histogram(param1Accepted_part, label="", color=1)
xlabel!("N")
ylabel!("accepted values")

h2 = histogram(param2Accepted_part, label="", color=2)
xlabel!("p")
# ylabel!("accepted values")

p2 = plot(h1, h2, layout=2)
savefig(p2, "ABC_N-p-lambdaInferred_noSample.pdf")

display(p2)

# %% codecell
# Inferred division rate
