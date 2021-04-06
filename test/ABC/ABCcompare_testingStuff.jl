include("../src/vafdyn.jl")
using .VAFDyn
using Plots
plotly()

# ===== True DFS =====
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

# ===== Fitted DFS =====
paramsFit = Dict(
    "N"=> 1000,
    "ρ"=>6.2,
    "ϕ"=>1.8,
    "μ"=>paramsTrue["μ"]
)
dfsFit = VAFDyn.DFreqspace(paramsFit["N"])
VAFDyn.evolveVAF(dfsFit, paramsFit, evolveTime)
sampFsFit = VAFDyn.sampler(dfsFit, sampleSize)

p1 = plot(dfsTrue.freqs_f, dfsTrue.n_f, yaxis=:log10, label="true")
plot!(dfsFit.freqs_f, dfsFit.n_f, label="fit")
plot!(sampFsTrue.freqs_f[2:end], sampFsTrue.n_f[2:end], label="true sampled")
plot!(sampFsFit.freqs_f[2:end], sampFsFit.n_f[2:end], label="fit sampled")
display(p1)