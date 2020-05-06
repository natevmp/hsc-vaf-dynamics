module TestVAF

include("../src/vafdyn.jl")

using .VAFDyn
using Plots
gr()

params = Dict(
    "ρ"=>1.,
    "μ"=>1.2,
    "N"=>400
)
evolveTime = 59

# discrete evolve
dfs = VAFDyn.DFreqspace(params["N"])
println("running Markov Chain")
VAFDyn.evolveVAF(dfs, params, evolveTime, 0.001)


# # diffusion evolve 
# cfs = VAFDyn.CFreqspace(2001)
# println("running PDE")
# VAFDyn.evolveVAF(cfs, params, evolveTime, 0.0001)

# plotting
h = plot(dfs.freqs_f[2:end], dfs.n_f[2:end], dpi=100, label="true", yaxis=:log10)
# plot!(cfs.freqs_f[2:end-1], cfs.n_f[2:end-1]/params["N"], linestyle=:dash, label="PDE")
# plot!(cfs.freqs_f[2:end-1], map(x -> params["μ"]/x, cfs.freqs_f[2:end-1]), linestyle=:dot, label="μ/f")
xlims!((0, 1))
ylims!((10^-1, 10^6))
xlabel!("variant frequency")
ylabel!("number of variants")
title!("Discrete Freqspace")
display(h)

sampledfs = VAFDyn.sampler(dfs, 80)[1]

plot!(h, sampledfs.freqs_f[2:end], sampledfs.n_f[2:end], label="sampled", yaxis=:log10)
display(h)

# savefig("figures/VAF_400HSC_59y.png")

end