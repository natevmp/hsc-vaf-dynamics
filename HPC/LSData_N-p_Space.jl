using JLD2
include("../src/inferencePipeline.jl")
using .InferencePipeline
include("../src/theory.jl")
using .Theory

# user params
paramsInput = [
    Dict(:tM => 3, :lVfs => 500, :Nmin => 2E4, :Nmax => 2E5),
    Dict(:tM => 5, :lVfs => 500, :Nmin => 4E4, :Nmax => 2.5E5),
    Dict(:tM => 8, :lVfs => 500, :Nmin => 4E4, :Nmax => 3.5E5),
    Dict(:tM => 10, :lVfs => 500, :Nmin => 5E4, :Nmax => 4E5)
]

# current sim params
simNum = parse(Int,ARGS[1])
tM = paramsInput[simNum][:tM]
lVfs = paramsInput[simNum][:lVfs]
Nmin = paramsInput[simNum][:Nmin]
Nmax = paramsInput[simNum][:Nmax]

# load data
@load "HPC/LSDataStatsBM.jld2" sampleSize SCBurdenHSC_CID nVHSC_f

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
_NDisc = range(Nmin, Nmax, length=20)
_pDisc = range(0.1, 0.9, length=20)

##
@time _N, _p, NOptInterpol_p, vaf1ErrorInterpol_p_N = InferencePipeline.calcNpSpace(paramsKnown, paramsEst, nVHSC_f, SCBurdenHSC_CID, _NDisc, _pDisc, lVfs, VafFit; verbose=true)


filename = "./LSData_NPSpaceInference_tM"*string(tM)*"_lVFS"*string(lVfs)*".jld2"
@save filename _N _p NOptInterpol_p

println("succes")
