using JLD2
include("../src/inferencePipeline.jl")
using .InferencePipeline
include("../src/theory.jl")
using .Theory
# include("../../src/compoundPoisson.jl")
# using .CompoundPoisson




# user params


tM = ARGS[1]
lVfs = [1000, 5000, 10000][parse(Int, ARGS[2])]


##
@load "LSDataStatsBM.jld2" sampleSize SCBurdenHSC_CID nVHSC_f

μKnown = 1.2
paramsKnown = Dict{String, Real}(
    "evolve time" => 59,
    "N initial" => 1,
    "sample size" => sampleSize,
    "mature time" => tM
)
# paramsEst = InferencePipeline.estimateRates(SCBurdenHSC_CID)
paramsEst = InferencePipeline.estimateRates(SCBurdenHSC_CID, μKnown)

VafFit = 1
_NDisc = range(1E4, 2E5, length=5)
_pDisc = range(0.1, 0.85, length=5)

##
@time _N, _p, NOptInterpol_p, vaf1ErrorInterpol_p_N = InferencePipeline.calcNpSpace(paramsKnown, paramsEst, nVHSC_f, SCBurdenHSC_CID, _NDisc, _pDisc, lVfs, VafFit; verbose=true)


filename = "./LSData_NPSpaceInference_tM"*string(tM)*"_lVFS"*string(lVfs)*".jld2"
@save filename _N _p NOptInterpol_p

println("succes")
