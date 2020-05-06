include("../src/vafdyn.jl")

using DelimitedFiles
using Plots
gr()
using .VAFDyn

# ===== Load Data =====
untypedM = readdlm("../data/Shearwater_calls_FDR0.95_all_muts.txt", '\t', Any; skipstart=1)
untypedM = untypedM[:, 5: end-1]
replace!(untypedM, "NA"=>0)

variants_var_col = Array{Int}(untypedM)

# ===== Order data =====
HSCMask_col = fill(false, size(variants_var_col, 2))
HSCMask_col[1:73] .= true
HSCMask_col[125:end] .= true
HPCMask_col = .!HSCMask_col

# mutational burden
mutBurden_col = sum(variants_var_col, dims=1)


# ===== variant allele frequencies =====
vaf_var = sum(variants_var_col, dims=2) / 140
vafHSC_var = sum(variants_var_col[:, HSCMask_col], dims=2) / sum(HSCMask_col)
vafHPC_var = sum(variants_var_col[:, HPCMask_col], dims=2) / sum(HPCMask_col)


# ===== plot VAF ===== #

gr()

p1 = histogram(vafHSC_var, yaxis=:log10, bins=sum(HSCMask_col), 
        ylims=(10^-0.3, 10^5), xlims=(0, 0.6), label="HSC", fillalpha=0.5, linealpha=0.4, dpi=200)
xlabel!("variant frequency")
ylabel!("number of variants")
title!("HSC")
display(p1)
println("number of HSC cells:", sum(HSCMask_col))

# ===== make VAF predictions =====

age = 59
dt = 0.0001
μ = 6.
ρ = 1.
params = Dict(
        "ρ"=>ρ,
        "μ"=>μ,
        "N"=>3000
)
dfs = vafMC(params, age, dt)
plot!(dfs.freqs_f, dfs.n_f, label="true: N = 3000", linewidth=2)
(sampfs, sampMutscalc) = sampler(dfs, sum(HSCMask_col))
plot!(sampfs.freqs_f, sampfs.n_f, label="sampled: N = 89", linewidth=2)
# title!("N="*string(params["N"])*"\nmu="*string(μ)*"\nrho="*string(ρ))
trueMutburden = mutBurden(dfs)
sampMutburden = mutBurden(sampfs)
println("true mutational burden: ",trueMutburden,"\n")
println("sample mutational burden measured: ",sampMutburden,"\n")
println("sample mutational burden calculated: ",sampMutburden,"\n")
display(p1)
# savefig("Figures/sampledDist")

# p3 = plot(dfs.freqs_f, dfs.n_f, label="true", linewidth=2)
# plot!(sampfs.freqs_f, sampfs.n_f, label="sampled", linewidth=2)
# display(p3)
