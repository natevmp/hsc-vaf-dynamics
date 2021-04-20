using JLD2, FileIO
using HypothesisTests
using Plots
pyplot()

include("../../src/theory.jl")
using .Theory

include("../../src/inferencePipeline.jl")
using .InferencePipeline

include("../../src/compoundPoisson.jl")
using .CompoundPoisson

##

# SCBurdenHSC_CID, SCBurdenHSCMean, SCBurdenHSCVar, freqs_f, nVHSC_f, sampleSize  =load("data/LSDataStats.jld2")
SCBurdenHSC_CID = load("data/LSDataStatsBM.jld2", "SCBurdenHSC_CID")

##
paramsEst = InferencePipeline.estimateRates(SCBurdenHSC_CID)

##
cpVals_id = CompoundPoisson.randComPois(paramsEst["divisions"], paramsEst["μ"], 50000)
##

figTest = histogram(SCBurdenHSC_CID, normalize=true, bins=25)
stephist!(cpVals_id, normalize=true)
display(figTest)

##

ADtest = KSampleADTest(SCBurdenHSC_CID, cpVals_id)



##

# paramsKnown = Dict{String, Real}(
#     "evolve time" => 59,
#     "N initial" => 1,
#     "sample size" => sampleSize,
#     "mature time" => 3
# )