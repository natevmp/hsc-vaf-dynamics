
include("../src/vafdyn.jl")
include("../src/vafSim.jl")

using Revise, ProgressMeter, JLD2
using .VAFDyn
using .VAFSim
using Plots
gr()

params = Dict(
    "N"=>450,
    "ρ"=>1.,
    "ϕ"=>3.,
    "μ"=>1.2
)

evolveTime = 100

# regular evolve
dfs = VAFDyn.DFreqspace(params["N"])
@time VAFDyn.evolveVAF(dfs, params, evolveTime)
#
# println("running regulars")
# @time for i in 1:10
#     dfs = VAFDyn.DFreqspace(params["N"])
#     VAFDyn.evolveVAF(dfs, params, evolveTime)
# end

dfs2 = VAFDyn.DFreqspace(params["N"])
@time VAFDyn.evolveVAFev(dfs2, params, evolveTime)

dfsVar = VAFDyn.DFreqspace(params["N"])
@time VAFDyn.evolveVAFvar(dfsVar, params, evolveTime)

println("running ev")
@time for i in 1:10
    dfs = VAFDyn.DFreqspace(params["N"])
    VAFDyn.evolveVAFev(dfs, params, evolveTime)
end

nSims = 500
# vafSimAv_n = zeros(Float64, params["N"]+1)
vaf_n_sim = zeros(Float64, params["N"]+1, nSims)
@showprogress for i in 1:nSims
    vafSim_n =
        VAFSim.birthDeathAlt(params["N"], params["μ"], 0.75, evolveTime, evolveTime, 4.)[1]
    vaf_n_sim[:, i] .= vafSim_n
end
vafSimAv_n = sum(vaf_n_sim, dims=2) / nSims
vafSimVar_n = sum(vaf_n_sim.^2, dims=2)/nSims - vafSimAv_n.^2

@save "vafSimsMeanVar.jld2" params vafSimAv_n vafSimVar_n dfs dfsVar

# @load "vafSimsMeanVar.jld2" vafSimAv_n vafSimVar_n

# plotting
p1 = scatter(dfs.freqs_f[2:end], vafSimAv_n[2:end], yaxis=:log10, label="simulation average", markersize=2.3, color=3, markerstrokewidth=0.2)

scatter!(dfs.freqs_f[2:end], dfs.n_f[2:end] .+ sqrt.(vafSimVar_n[2:end]), label="simulation std", linewidth=1, color=4, markersize=2.3, markerstrokewidth=0.2)
scatter!(dfs.freqs_f[2:end], dfs.n_f[2:end] .- sqrt.(vafSimVar_n[2:end]), label="", linewidth=1, color=4, markersize=2.3, markerstrokewidth=0.2)

plot!(dfs.freqs_f[2:end], dfs.n_f[2:end], label="Markov Chain evolve", linewidth=3, color=1)

plot!(dfs2.freqs_f[2:end], dfs2.n_f[2:end], label="expected value evolve", linewidth=2, linestyle=:dash, color=2)

plot!(dfs.freqs_f[2:end], dfs.n_f[2:end] .+ sqrt.(dfsVar.n_f[2:end]), label="standard deviation evolve", linewidth=2, linestyle=:dashdot, color=2)
plot!(dfs.freqs_f[2:end], dfs.n_f[2:end] .- sqrt.(dfsVar.n_f[2:end]), label="", linewidth=2, linestyle=:dashdot, color=2)

xlims!((0, 0.05))
ylims!((10^2, 10^4))
xlabel!("variant frequency")
ylabel!("number of variants")
display(p1)

p3 = scatter(dfs.freqs_f[2:end], dfsVar.n_f[2:end], yaxis=:log10, linewidth=2, label="prediction")
scatter!(dfs.freqs_f[2:end], vafSimVar_n[2:end], label="simulations")
xlims!((0, 0.2))
ylims!((10^1, 2*10^5))
xlabel!("variant frequency")
ylabel!("σ^2 on #variants")
display(p3)

# savefig(p3, "figures/varianceCompare.pdf")
#
# p4 = plot(dfs.freqs_f[2:end], sqrt.(dfsVar.n_f[2:end]), yaxis=:log10, linewidth=2)
# scatter!(dfs.freqs_f[2:end], sqrt.(vafSimAv_n[2:end]))
# xlims!((0, 0.2))
# ylims!((10^0, 2*10^2))
# xlabel!("variant frequency")
# ylabel!("σ on number of variants")
# display(p4)

# savefig("figures/varianceComparison.pdf")

# p4 = plot(dfs2.freqs_f[2:end], dfs2.n_f[2:end], yaxis=:log10, linewidth=2, label="ev")
# plot!(dfs.freqs_f[2:end], dfsVar.n_f[2:end], label="var", linewidth=2)
# xlims!((0, 0.2))
# ylims!((10^1, 2*10^5))
# xlabel!("variant frequency")
# ylabel!("σ^2 on #variants")
# display(p4)

# p2 = scatter(dfs.freqs_f[2:end], vafSimAv_n[2:end], label="simulation average", markersize=2.3, color=3, markerstrokewidth=0.2)
#
# scatter!(dfs.freqs_f[2:end], vafSimAv_n[2:end] .+ sqrt.(vafSimVar_n[2:end]), label="simulation std", linewidth=1, color=4, markersize=2.3, markerstrokewidth=0.2)
# scatter!(dfs.freqs_f[2:end], vafSimAv_n[2:end] .- sqrt.(vafSimVar_n[2:end]), label="", linewidth=1, color=4, markersize=2.3, markerstrokewidth=0.2)
#
# plot!(dfs.freqs_f[2:end], dfs.n_f[2:end], label="Markov Chain evolve", linewidth=3, color=1)
#
# plot!(dfs2.freqs_f[2:end], dfs2.n_f[2:end], label="expected value evolve", linewidth=2, linestyle=:dash, color=2)
#
# plot!(dfs.freqs_f[2:end], dfs.n_f[2:end] .+ sqrt.(dfsVar.n_f[2:end]), label="standard deviation evolve", linewidth=2, linestyle=:dashdot, color=2)
# plot!(dfs.freqs_f[2:end], dfs.n_f[2:end] .- sqrt.(dfsVar.n_f[2:end]), label="", linewidth=2, linestyle=:dashdot, color=2)
#
# xlims!((0, 0.2))
# ylims!((10^2, 2*10^3))
# xlabel!("variant frequency")
# ylabel!("number of variants")
# display(p2)
