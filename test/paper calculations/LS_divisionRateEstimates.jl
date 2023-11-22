using JLD2, LaTeXStrings, Glob, FileIO
include("../../src/inferencePipeline.jl")
using .InferencePipeline
include("../../src/theory.jl")
using .Theory

## load data
@load "HPC/LSDataStatsBM.jld2" sampleSize SCBurdenHSC_CID nVHSC_f

NOptInterpol_nP_P = Vector{Float64}[]


dataFolder = "data/LSData_NPSpace_21-08-05/"
filenameData = glob(dataFolder*"LSData_NPSpaceInference_tM5_lVFS500_pureGrowth100*")

@load filenameData[1] NOptInterpol_p paramsTot _p _N vaf1ErrorInterpol_p_N

##
μKnown = 1.2

paramsKnown = Dict{String, Real}(
    "evolve time" => 59,
    "N initial" => 1,
    "sample size" => sampleSize,
    "mature time" => paramsTot["mature time"],
    "pure births" => paramsTot["pure births"],
    "μ" => μKnown,
)
paramsEst = InferencePipeline.estimateRates(SCBurdenHSC_CID, μKnown)
paramsIn = merge(paramsKnown, paramsEst)
##
_λ = similar(_p)
_seλ = similar(_p)
for (i,p) in enumerate(_p)
    N = NOptInterpol_p[i]
    paramsParticle = InferencePipeline.createTestParams(paramsIn, N, p)
    _λ[i] = paramsParticle["λ"]
    _seλ[i] = paramsParticle["seλ"]
end

# λMean = _λ[1] + (_λ[end] - _λ[1])/2
λMean = round(mean(_λ),digits=1)
λLB = round(_λ[1] - _seλ[1],digits=1)
λUB = round(_λ[end] + _seλ[end],digits=1)

##
println("estimated division rate is: λ = $λMean ($λLB - $λUB) ")

