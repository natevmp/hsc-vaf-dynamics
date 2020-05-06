include("../src/vafdyn.jl")

using .VAFDyn
using Plots
pyplot()


evolveTime = 59
poolSize = 400
ρ = 1.
μ = 1.2

params = Dict(
    "ρ"=>ρ,
    "μ"=>1.,
    "N"=>poolSize
)
# discrete evolve
dfs = VAFDyn.DFreqspace(params["N"])
println("running Markov Chain")
VAFDyn.evolveVAF(dfs, params, evolveTime, 0.001)

# plotting
h = plot(dfs.freqs_f[2:end], dfs.n_f[2:end], dpi=200, label="μ = "*string(1.0), yaxis=:log10)
xlims!((0, 1))
ylims!((10^-1, 10^6))
xlabel!("variant frequency")
ylabel!("number of variants")
title!("Discrete Freqspace")

# test mutation rate influence
for μ=4.:3.:13.
    params = Dict(
        "ρ"=>ρ,
        "μ"=>μ,
        "N"=>poolSize
    )
    # discrete evolve
    dfs = VAFDyn.DFreqspace(params["N"])
    println("running Markov Chain")
    VAFDyn.evolveVAF(dfs, params, evolveTime, 0.001)
    # plotting
    plot!(dfs.freqs_f[2:end], dfs.n_f[2:end], label="μ = "*string(μ), yaxis=:log10)
end

display(h)

# savefig("Figures/VAF_400HSC_59y_muCompare.png")

# ====== λ compare ======
evolveTime = 59
poolSize = 400
ρ = 1.
μ = 1.2

params = Dict(
    "ρ"=>ρ,
    "μ"=>μ,
    "N"=>poolSize
)
# discrete evolve
dfs = VAFDyn.DFreqspace(params["N"])
println("running Markov Chain")
VAFDyn.evolveVAF(dfs, params, evolveTime, 0.001)

# plotting
h = plot(dfs.freqs_f[2:end], dfs.n_f[2:end], dpi=200, label="ρ = "*string(ρ), yaxis=:log10)
xlims!((0, 1))
ylims!((10^-1, 10^6))
xlabel!("variant frequency")
ylabel!("number of variants")
title!("Discrete Freqspace")

# test mutation rate influence
for ρ=2.:1.:5.
    params = Dict(
        "ρ"=>ρ,
        "μ"=>μ,
        "N"=>poolSize
    )
    # discrete evolve
    dfs = VAFDyn.DFreqspace(params["N"])
    println("running Markov Chain")
    VAFDyn.evolveVAF(dfs, params, evolveTime, 0.001)
    # plotting
    plot!(dfs.freqs_f[2:end], dfs.n_f[2:end], label="ρ = "*string(ρ), yaxis=:log10)
end
display(h)

# savefig("Figures/VAF_400HSC_59y_rhoCompare.png")