using JLD2
include("../../src/inferencePipeline.jl")
using .InferencePipeline
include("../../src/theory.jl")
using .Theory
# include("../../src/compoundPoisson.jl")
# using .CompoundPoisson




# user params

tM = 10
lVfs = 500
_nPure = [0, 10, 30, 60, 100, 500]

##
@load "data/LSDataStatsBM.jld2" sampleSize SCBurdenHSC_CID nVHSC_f

μKnown = 1.2
paramsKnown = Dict{String, Real}(
    "evolve time" => 59,
    "N initial" => 1,
    "sample size" => sampleSize,
    "mature time" => tM,
    "pure births" => _nPure[6],
)
# paramsEst = InferencePipeline.estimateRates(SCBurdenHSC_CID)
paramsEst = InferencePipeline.estimateRates(SCBurdenHSC_CID, μKnown)
testPars = merge(paramsKnown,paramsEst)
testPars["N final"] = 2E4
testPars["p"] = 0.4
Theory.getλFromTotalDivisions(testPars)

##
# cpVals_id = CompoundPoisson.randComPois(paramsEst["divisions"], paramsEst["μ"], 50000)
# ##
# figTest = histogram(SCBurdenHSC_CID, normalize=true, bins=25)
# stephist!(cpVals_id, normalize=true)
# display(figTest)

##

VafFit = 1
_NDisc = range(1E4, 3.5E5, length=4)
_pDisc = range(0.1, 0.85, length=4)

##
@time _N, _p, NOptInterpol_p, vaf1ErrorInterpol_p_N = InferencePipeline.calcNpSpace(paramsKnown, paramsEst, nVHSC_f, SCBurdenHSC_CID, _NDisc, _pDisc, lVfs, VafFit; verbose=true)


##

# filename = "./LSData_NPSpaceInference_tM"*string(tM)*"_lVFS"*string(lVfs)*".jld2"
# @save filename _N _p NOptInterpol_p

##
using Plots
pyplot()

fig1 = plot(_p, NOptInterpol_p, label="vaf fit:"*string(VafFit))
xlabel!("p")
ylabel!("N")
title!("vfsL = "*string(lVfs))

display(fig1)
# savefig(fig1, "NP_fit"*string(VafFit)*".pdf")
