using Revise
using DelimitedFiles
using Statistics
using Distributions
using LaTeXStrings

##
using Plots
gr()

##
include("../src/vafdyn.jl")
using .VAFDyn

## ==== Load Data ====
untypedM = readdlm("./data/Shearwater_calls_FDR0.95_all_muts.txt", '\t', Any; skipstart=1)
untypedM = untypedM[:, 5: end-1]
replace!(untypedM, "NA"=>0)
variants_var_col = Array{Int}(untypedM)

HSCMask_col = fill(false, size(variants_var_col, 2))
HSCBMMask_col = fill(false, size(variants_var_col, 2))
HSCBMMask_col[1:73] .= true
HSCMask_col[1:73] .= true
HSCMask_col[125:end] .= true
HPCMask_col = .!HSCMask_col;
mutBurden_col = sum(variants_var_col, dims=1);

vaf_var = sum(variants_var_col, dims=2) / 140
vafHSC_var = sum(variants_var_col[:, HSCMask_col], dims=2);
vafHPC_var = sum(variants_var_col[:, HPCMask_col], dims=2) / sum(HPCMask_col);

vafHSCBM_var = sum(variants_var_col[:, HSCBMMask_col], dims=2) / sum(HSCBMMask_col)

vafAbsHSC_var = sum(variants_var_col[:, HSCMask_col], dims=2);
vafAbsHSCBM_var = sum(variants_var_col[:, HSCBMMask_col], dims=2);



## ==== Build VAF spectra ====
function makeHist(varPrev_var, N::Integer)
    n_f = zeros(Int64, N+1)
    for m in 1:N
        n_f[1+m] = length(findall(x -> x==m, varPrev_var))
    end
    return n_f
end
vafHSC_f = makeHist(vafAbsHSC_var, sum(HSCMask_col))
vafHSCBM_f = makeHist(vafAbsHSCBM_var, sum(HSCBMMask_col))

# println(vafHSC_f[2])



## ==== Fit VAF 1 data to equilibrium ====

# function eqVAF1(p::Real, μ::Real, N::Real)
#     (2 + p/(1-p)) * μ * N
# end

μ1 = 1.2
μ2 = 4.2
S = sum(HSCMask_col)
vRef = vafHSC_f[2]

pVals_i = 0:0.01:0.99
NVals_j = 1*10^4:10^2:5*10^4

# v1_p_N = Array{Float64, 2}(undef, length(pVals_i), length(NVals_j))
# vDist_p_N = Array{Float64, 2}(undef, length(pVals_i), length(NVals_j))
nSolved1_p = Array{Float64,1}(undef, length(pVals_i))
nSolved2_p = Array{Float64,1}(undef, length(pVals_i))
@time for (i,p) in enumerate(pVals_i)
    nSolved1_p[i] = vRef / ((2 + p/(1-p))*μ1)
    nSolved2_p[i] = vRef / ((2 + p/(1-p))*μ2)
end


## ===== Plotting =====

p2 = plot(pVals_i, nSolved1_p,
        linewidth=2,
#         dpi=300,
        label="μ = "*string(μ1))
plot!(pVals_i, nSolved2_p,
        linewidth=2,
        label="μ = "*string(μ2))
xlabel!("p")
ylabel!("N")
title!("f=1/N equilibrium fit")
display(p2)
# savefig(p2, "../figures/p-N_dependence.png")

