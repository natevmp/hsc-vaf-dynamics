include("../vafdyn.jl")
include("../theory.jl")
using Distributions
using .VAFDyn, .Theory

##
params = Dict(
    "N initial" => 400,
    "N final" => 400,
    "μ" => 1.2 * 500 * 1e-9,
    "ρ" => 5,
    "ϕ" => 5,
    "mature time" => 5,
    "evolve time" => 60,
)
params["growth rate"] = Theory.expGrowthRateFromNT(params["N initial"], params["N final"], params["mature time"])

## ===== Discrete evolve =====
dfs = VAFDyn.DFreqspace(params["N final"])
dfs.n_f[2] = 1
@time vec = VAFDyn.evolveCloneGrowingPopAlt(dfs, params, params["evolve time"])
display(vec)
println("sum: ")
println(sum(vec[2:end-2]))
println("")

##
dfs = VAFDyn.DFreqspace(params["N final"])
dfs.n_f[2] = 1
@time VAFDyn.evolveCloneGrowingPop(dfs, params, params["evolve time"])
display(dfs.n_f)
println("sum: ")
println(sum(dfs.n_f[2:end-1]))
println("")

## ===== Compare with Direct Poisson prediction of occurrences
driverProbRate(t) = params["N final"] * (2params["ρ"] + params["ϕ"]) * params["μ"] * t
driverProbAge(t) = 1 - pdf(Poisson(driverProbRate(t)), 0)

println("Poisson prob: ", driverProbAge(60))
println("MC evolved prob: ", 1-vec[1])

## ===== PDE evolve =====

lVfs = 401
vfs = VAFDyn.VFreqspace(params["N final"], lVfs)
vfs.n_f[2] = 1
@time VAFDyn.evolveCloneGrowingPop(vfs, params, params["evolve time"])
# sum(vfs.n_f[2:end] .* (vfs.freqs_f[2:end] .- vfs.freqs_f[1:end-1]) * params["N final"])
display(vfs.n_f)
sum(vfs.n_f[2:end-1])
#
# dfs2 = VAFDyn.makeDFSfromVFS(vfs, params["N final"], normalize=false)
# display(dfs2.n_f)
# sum(dfs2.n_f[2:end-1])
# ##
# # println(vfs.n_f)


##
using Plots

plot(dfs.n_f[2:end-1], yscale=:log10)
plot!(vec[2:end-2])
plot!(vfs.n_f[2:end-1])
plot!(dfs2.n_f[2:end-1])