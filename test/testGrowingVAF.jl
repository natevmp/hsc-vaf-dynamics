include("../src/vafdyn.jl")

using .VAFDyn
using Plots
pyplot()

par = Dict(
    "ρ"=>1.,
    "γ"=>0.06,
    "μ"=>1.2,
    "N0"=>200
)

evolveTime = 20
# final compartment size
nT = par["N0"] * exp(par["γ"]*evolveTime)
println("Final population size: "*string(nT))

cfs = VAFDyn.CFreqspace(2001)
println("running PDE")
VAFDyn.evolveGrowingVAF(cfs, par, evolveTime, 0.0001)

# plotting
h = plot(cfs.freqs_f[2:end-1], cfs.n_f[2:end-1]/par["N0"], dpi=200, 
        label="PDE; t="*string(evolveTime), yaxis=:log10)
plot!(cfs.freqs_f[2:end-1], map(x -> par["μ"]/x, cfs.freqs_f[2:end-1]), label="μ/x")
plot!(cfs.freqs_f[2:end-1], map(x -> par["μ"]/x^2, cfs.freqs_f[2:end-1]), label="μ/x^2")
xlims!((0, 1))
# ylims!((10^-10, 10^2))
xlabel!("variant frequency")
ylabel!("number of variants")
title!("Discrete Freqspace: t = " * string(evolveTime))
display(h)

# savefig("Figures/GrowingVAFtest1.png")

