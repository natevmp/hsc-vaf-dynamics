include("../src/vafdyn.jl")
using .VAFDyn
include("../src/theory.jl")
using .Theory
##
using Plots
pyplot()
##

params1 = Dict(
    "N initial" => 500,
    "N final" => 500,
    "λ" => 6,
    "p" => 0.5,
    "μ" => 2,
    "evolve time" => 50,
    "mature time" => 5,
    "sample size" => 25
)

params2 = deepcopy(params1)
params2["mature time"] = 10

params2 = deepcopy(params1)
params2["mature time"] = 15

params3 = deepcopy(params1)
params3["mature time"] = 20

extendParams!(params1)
extendParams!(params2)
extendParams!(params3)

##

dfs1 = VAFDyn.DFreqspace(params1["N final"])
@time VAFDyn.evolveGrowingVAF(dfs1, params1, params1["evolve time"])
# dfsS1 = VAFDyn.sampler(dfs1, params1["sample size"])

##
dfs2 = VAFDyn.DFreqspace(params2["N final"])
@time VAFDyn.evolveGrowingVAF(dfs2, params2, params2["evolve time"])
# dfsS2 = VAFDyn.sampler(dfs2, params2["sample size"])

dfs3 = VAFDyn.DFreqspace(params3["N final"])
@time VAFDyn.evolveGrowingVAF(dfs3, params3, params3["evolve time"])
# dfsS3 = VAFDyn.sampler(dfs3, params3["sample size"])

##
vfs1 = VAFDyn.VFreqspace(params1["N final"], 50)
@time VAFDyn.evolveGrowingVAF(vfs1, params1, params1["evolve time"])
vDfs1 = VAFDyn.makeDFSfromVFS(vfs1, params1["N final"])
# vDfsS1 = VAFDyn.sampler(vDfs1, params1["sample size"])

##
fig1a = plot(dfs1.freqs_f[2:end-1], dfs1.n_f[2:end-1],
    yscale=:log10,
    linewidth=2,
    label="discrete evolve"
)
plot!(vDfs1.freqs_f[2:end-1], vDfs1.n_f[2:end-1],
linewidth=1.5,
    linestyle=:dash,
    label="continuous evolve"
)
# ylims!(1E-1, 5E3)
ylabel!("number of variants")

fig1b = plot(dfs1.freqs_f[2:end], (dfs1.n_f-vDfs1.n_f)[2:end])
ylabel!("continuous evolve error")

fig1 = plot(fig1a, fig1b, layout=(2,1), size=(600,800))
display(fig1)

##
fig2a = plot(dfs1.freqs_f[2:end-1], dfs1.n_f[2:end-1],
    yscale=:log10,
    label="t_m = "*string(params1["mature time"])
)

plot!(dfs2.freqs_f[2:end-1], dfs2.n_f[2:end-1],
    label="t_m = "*string(params2["mature time"])
)

plot!(dfs3.freqs_f[2:end-1], dfs3.n_f[2:end-1],
    label="t_m = "*string(params3["mature time"])
)

fig2b = plot(dfs1.freqs_f[2:end-1], (dfs1.n_f-dfs2.n_f)[2:end-1],
    color=2,
    label="t_m="*string(params1["mature time"])*" - t_m="*string(params2["mature time"])
)
plot!(dfs1.freqs_f[2:end-1], (dfs1.n_f-dfs3.n_f)[2:end-1],
    color=3,
    label="t_m="*string(params1["mature time"])*" - t_m="*string(params3["mature time"])
)

# ylims!(1E-5, 3.2E3)
fig2 = plot(fig2a, fig2b, layout=(2,1), size=(600,800))
display(fig2)