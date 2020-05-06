module tester

include("../src/vafdyn.jl")

using .VAFDyn
using Plots
pyplot()

params = Dict(
    "ρ"=>1.,
    "μ"=>1.2,
    "N"=>400
)

# discrete evolve
dfs = VAFDyn.DFreqspace(params["N"])
VAFDyn.addClones(dfs, Integer(round(params["N"]*0.3)), 1.)
VAFDyn.evolveVAF(dfs, params, 10, 0.001, addClones=false)

# diffusion evolve 
cfs = VAFDyn.CFreqspace(1001)
VAFDyn.addClones(cfs, 0.3, 1.)
VAFDyn.evolveVAF(cfs, params, 10, 0.0001, addClones=false)

# kimura solution
probCloneKimura_f = VAFDyn.evolveCloneKimura(cfs.freqs_f, 0.3, 10, params["N"], 60)

# plotting
h = plot(dfs.freqs_f[2:end], dfs.n_f[2:end], dpi=200, label="MC")
plot!(cfs.freqs_f[2:end-1], cfs.n_f[2:end-1]/params["N"], linestyle=:dash, label="PDE")
plot!(cfs.freqs_f, probCloneKimura_f/params["N"], linestyle=:dot, label= "Kimura solution")
xlims!((0, 1))
xlabel!("variant frequency")
ylabel!("probability")
title!("Discrete Freqspace")
display(h)

# savefig("figures/singleCloneCompare_400HSC.png")




end