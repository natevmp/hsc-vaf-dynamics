include("../src/vafdyn.jl")
using .VAFDyn
include("../src/vafSim.jl")
using .VAFSim
include("../src/theory.jl")
using .Theory

using Plots
pyplot()

##

params = Dict(
    "N initial" => 1,
    "N final" => 500,
    "μ" => 2.,
    "λ" => 5.,
    "p" => 0.5,
    "sample size" => 100,
    "mature time" => 10,
    "evolve time" => 10
)
Theory.extendParams!(params)

println(params)

##
dfs = DFreqspace(params["N final"])
@time VAFDyn.evolveGrowingVAF(dfs, params, params["evolve time"])

## sims

_, _, nVarSim_f, nVarSimS_f, _, _, _, _ = VAFSim.birthDeathFixedGrowth(params, params["evolve time"], 1)



## 1/f^2 model

# nDetmodelF2(f) = 2*params["μ"]*( 1+(params["ρ"]+params["ϕ"]/2)/params["growth rate"] )*( 1/f^2 - 1 )

function nDetGenModelF2(f)
    μ = params["μ"]
    ρ = params["ρ"]
    ϕ = params["ϕ"]
    γ = params["growth rate"]
    2μ*(1+(ρ+ϕ/2)/γ) * (1/f^2 - 1)
end

nDetMoranModel(f) = 2/(params["N final"]*f) - 2/(params["N final"]*f+1)

nDetGenModel_f = nDetGenModelF2.(dfs.freqs_f)
nDetMoranModel_f = nDetMoranModel.(dfs.freqs_f)

##

fig1 = plot(dfs.freqs_f[2:end-1], dfs.n_f[2:end-1]/dfs.n_f[2], yscale=:log10)
plot!(dfs.freqs_f[2:end-1], nVarSim_f[2:end-1]/nVarSim_f[2], seriestype=:sticks)
plot!(dfs.freqs_f[2:end-1], nDetGenModel_f[2:end-1]/nDetGenModel_f[2])
plot!(dfs.freqs_f[2:end-1], nDetMoranModel_f[2:end-1]/nDetMoranModel_f[2])
xlims!(0,1)
ylims!(1E-5, 1E0)
display(fig1)
##
fig1 = plot(dfs.freqs_f[2:end-1], dfs.n_f[2:end-1]/dfs.n_f[2], yscale=:log10)
plot!(dfs.freqs_f[2:end-1], nVarSim_f[2:end-1]/nVarSim_f[2], seriestype=:sticks)
plot!(dfs.freqs_f[2:end-1], nDetGenModel_f[2:end-1]/nDetGenModel_f[2])
plot!(dfs.freqs_f[2:end-1], nDetMoranModel_f[2:end-1]/nDetMoranModel_f[2])
xlims!(0,0.1)
ylims!(1E-5, 1E0)
display(fig1)
