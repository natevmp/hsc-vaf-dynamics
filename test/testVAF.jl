
include("../src/vafdyn.jl")


using .VAFDyn
using Plots, JLD
gr()

println("starting")

params = Dict(
    "N"=>450,
    "ρ"=>1.,
    "ϕ"=>3.,
    "μ"=>1.2
)
evolveTime = 100

# discrete evolve
println("running singles")
@time begin
    dfs = VAFDyn.DFreqspace(params["N"])
    VAFDyn.evolveVAF(dfs, params, evolveTime)
end
@time begin
    dfs = VAFDyn.DFreqspace(params["N"])
    VAFDyn.evolveVAF(dfs, params, evolveTime)
end
@time begin
    dfs = VAFDyn.DFreqspace(params["N"])
    VAFDyn.evolveVAF(dfs, params, evolveTime)
end
@time begin
    dfs = VAFDyn.DFreqspace(params["N"])
    VAFDyn.evolveVAF(dfs, params, evolveTime)
end
@time begin
    dfs = VAFDyn.DFreqspace(params["N"])
    VAFDyn.evolveVAF(dfs, params, evolveTime)
end

println("running loop 1")
@time for i in 1:50
    dfs = VAFDyn.DFreqspace(params["N"])
    VAFDyn.evolveVAF(dfs, params, evolveTime)
end

params = Dict(
    "N"=>352,
    "ρ"=>1.,
    "ϕ"=>3.,
    "μ"=>1.2
)

println("running loop 2")
@time for i in 1:50
    dfs = VAFDyn.DFreqspace(params["N"])
    VAFDyn.evolveVAF(dfs, params, evolveTime)
end

# # plotting
# h = plot(dfs.freqs_f[2:end], dfs.n_f[2:end], dpi=100, label="true", yaxis=:log10)
# # plot!(cfs.freqs_f[2:end-1], cfs.n_f[2:end-1]/params["N"], linestyle=:dash, label="PDE")
# # plot!(cfs.freqs_f[2:end-1], map(x -> params["μ"]/x, cfs.freqs_f[2:end-1]), linestyle=:dot, label="μ/f")
# xlims!((0, 1))
# ylims!((10^-1, 10^6))
# xlabel!("variant frequency")
# ylabel!("number of variants")
# title!("Discrete Freqspace")
# display(h)

sampledfs = VAFDyn.sampler(dfs, 80)

# plot!(h, sampledfs.freqs_f[2:end], sampledfs.n_f[2:end], label="sampled", yaxis=:log10)
# display(h)

# savefig("figures/VAF_400HSC_59y.png")

save("vafTest.jld", "params", params, "sampledfs", sampledfs, "evolveTime", evolveTime, "dfs", dfs)
