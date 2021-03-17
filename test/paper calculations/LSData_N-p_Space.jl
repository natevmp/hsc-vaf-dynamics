include("../../src/inferencePipeline.jl")
using .InferencePipeline
using JLD2
include("../../src/compoundPoisson.jl")
using .CompoundPoisson
include("../../src/theory.jl")
using .Theory

##
@load "data/LSDataStats.jld2"

##

paramsKnown = Dict{String, Real}(
    "evolve time" => 59,
    "N initial" => 1,
    "sample size" => sampleSize,
    "mature time" => 3
)


_NDisc = [501,1E5,3E5,7E5]
_pDisc = [0.1, 0.35, 0.65, 0.85]

##
using Plots
pyplot()
##

paramsEst = InferencePipeline.estimateRates(SCBurdenHSC_CID)

##
cpVals_id = CompoundPoisson.randComPois(paramsEst["divisions"], paramsEst["μ"], 50000)
##
figTest = histogram(SCBurdenHSC_CID, normalize=true, bins=25)
stephist!(cpVals_id, normalize=true)
display(figTest)

##
@time _N, _p, NOptInterpol_p, vaf1ErrorInterpol_p_N = calcNpSpace(paramsKnown, nVHSC_f, SCBurdenHSC_CID, _NDisc, _pDisc, 200; verbose=true)


