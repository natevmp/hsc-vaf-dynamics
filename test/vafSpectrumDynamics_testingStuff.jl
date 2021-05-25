
include("../src/vafdyn.jl")
using .VAFDyn
include("../src/theory.jl")
using .Theory
using OrdinaryDiffEq
using JLD2, Glob
using Plots
pyplot()
using Statistics

##
fileNames_ = glob("singlePatientFullSim_Ni10000_Nf10000_*.jld2", "./data/Nf10000")
nFiles = length(fileNames_)

vaf1Sim_sim = Int[]
for fname in fileNames_
    @load fname nVarSim_f paramsTrue
    global params = paramsTrue
    push!(vaf1Sim_sim, nVarSim_f[2])
end
vaf1SimAv = mean(vaf1Sim_sim)
vaf1SimVar = var(vaf1Sim_sim)
√vaf1SimVar

println(params)
Theory.extendParams!(params)

##

_lVfs = 100:50:500
vafS1_lVfs = Float64[]
vfs1_lVfs = Float64[]
for lVfs in _lVfs
    vfs = VAFDyn.VFreqspace(params["N final"], lVfs)
    @time VAFDyn.evolveGrowingVAF(vfs, params, params["evolve time"])
    dfs = VAFDyn.makeDFSfromVFS(vfs, params["N final"])
    # dfsS = VAFDyn.sampler(dfs, params["sample size"])
    push!(vfs1_lVfs, vfs.n_f[2])
    push!(vafS1_lVfs, dfs.n_f[2])
end


##

plot(_lVfs, vfs1_lVfs/ params["N final"])
xlabel!("lVFS")
ylabel!("vaf1")


