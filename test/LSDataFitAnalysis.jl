include("../src/vafdyn.jl")

using DelimitedFiles, JLD2
using ApproxBayes
using Distributions, Distances
using .VAFDyn

using Plots
gr()

## ============================== Load Lee-Six Data ==============================

untypedM = readdlm("data/Shearwater_calls_FDR0.95_all_muts.txt", '\t', Any; skipstart=1)
untypedM = untypedM[:, 5: end-1]
replace!(untypedM, "NA"=>0)

variants_var_col = Array{Int}(untypedM)

## ===== Order data =====
HSCMask_col = fill(false, size(variants_var_col, 2))
HSCMask_col[1:73] .= true   # bone marrow HSCs
# HSCMask_col[125:end] .= true  # Peripheral blood HSCs
HPCMask_col = .!HSCMask_col
sampleSize = sum(HSCMask_col)

# mutational burden
mutBurden_col = sum(variants_var_col, dims=1)

mutHSCBurdenAv = mean(mutBurden_col[HSCMask_col])
mutHSCBurdenVar = var(mutBurden_col[HSCMask_col])

# ===== variant allele frequencies =====
vaf_var = sum(variants_var_col, dims=2) / 140
prevHSC_var = sum(variants_var_col[:, HSCMask_col], dims=2)

function makeHist(varPrev_var, N::Integer)
    n_f = zeros(Int64, N+1)
    for m in 1:N
        n_f[1+m] = length(findall(x -> x==m, varPrev_var))
    end
    return n_f
end

vafHSCData_f = makeHist(prevHSC_var, sampleSize)
freqs_f = range(0, 1; length=sampleSize+1)

# dataParams = Dict(
#     "evolveTime" => 59,
#     "sampleSize" => sum(HSCMask_col),
#     "mutMean" => mutHSCBurdenAv,
#     "mutVar" => mutHSCBurdenVar,
#     "μ" => mutHSCBurdenVar/mutHSCBurdenAv - 1
# )

## ============================== Load fit data ==============================

# data = load("test/ABC_VAF1fit.jld")
@load "test/ABC_VAF1fit_mu4.29.jld2"
# @load "test/ABC_VAF1fit_mu1.2.jld2"

μString = string(round(dataParams["μ"],digits=2))

## ============================== Plot ABC results ==============================

# h1 = histogram(param1Accepted_part, label="", color=1, bins=20)
# xlabel!("N")
# ylabel!("accepted values")
# h2 = histogram(param2Accepted_part, label="", color=2, bins=20)
# xlabel!("p")
# # ylabel!("accepted values")
# p1 = plot(h1, h2, layout=2)

p2 = scatter(param2Accepted_part, param1Accepted_part,
    # size=(210,300),
    markersize=8,
    xtickfont=font(13), 
    ytickfont=font(13), 
    guidefont=font(16), 
    legendfont=font(13),
    label="accepted parameter pair")
xlabel!("p (between 0 and 1)")
ylabel!("N")
title!("μ = "*μString)

p3 = plot(p1, p2, layout=(2,1), size=(600, 800))
# display(p3)
# savefig(p3, "ABC_N-p_Sampled.pdf")

display(p2)
savefig(p2, "figures/ABCfitting/ABC_VAF1_mu"*μString*"_N-p_Sampled.pdf")

## ===== Params =====

display(dataParams)
# dataParams["λPar"] = dataParams["mutMean"]^2 / ((dataParams["mutVar"] - dataParams["mutMean"]))

evolveTime = 59
userParams = Dict(
    "p"=>0.8,
    "μ"=>dataParams["μ"],
    "N"=>13000
)
# userParams["λ"] = dataParams["mutMean"]^2 / ((dataParams["mutVar"] - dataParams["mutMean"])*(2-userParams["p"])*evolveTime)
userParams["λ"] = dataParams["λPar"] / ((2-userParams["p"])*evolveTime)

display(userParams)

## ============================== Compare LSdata and Fit ============================== 

evoParams = Dict(
    "ρ"=>userParams["λ"]*(1-userParams["p"]),
    "μ"=>userParams["μ"],
    "ϕ"=>userParams["λ"]*userParams["p"],
    "N"=>userParams["N"]
)

# ==== create fitted VAF ====

vfsFit = VAFDyn.VFreqspace(evoParams["N"],1001)
@time VAFDyn.evolveVAFfd(vfsFit, evoParams, evolveTime);
dfsFit = VAFDyn.makeDFSfromVFS(vfsFit, evoParams["N"])
@time dfsFitSampled = VAFDyn.sampler(dfsFit, dataParams["sampleSize"])

## ==== Plot fitted and measured VAFs ====
μString = string(round(dataParams["μ"],digits=2))
pString = string(userParams["p"])
NString = string(userParams["N"])


p0 = plot(freqs_f, vafHSCData_f, yaxis=:log10, ylims=(10^-0.3, 10^5),
    label="",
    linestyle=:solid,
    color=1,
    linewidth=2)
scatter!(freqs_f, vafHSCData_f,
    label="LS data",
    linestyle=:solid,
    color=1,
    linewidth=3)
plot!(dfsFitSampled.freqs_f, dfsFitSampled.n_f, 
    label="fit: N = "*NString*"; p = "*pString,       
    linestyle=:dash,
    color=2,
    linewidth=3)
xlims!(0,0.35)
xlabel!("variant frequency")
ylabel!("number of variants")
title!("μ = "*μString)
display(p0)

# savefig("figures/ABCfitting/ABC_VAF1_fitCompare_mu"*μString*"_p"*pString*".pdf")