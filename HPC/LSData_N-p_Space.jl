using JLD2
include("../src/inferencePipeline.jl")
using .InferencePipeline
include("../src/theory.jl")
using .Theory

# user params
paramsInput = [
    Dict(:nPure => 500, :tM => 10, :lVfs => 500, :Nmin => 5E3, :Nmax => 2E5),
    Dict(:nPure => 100, :tM => 10, :lVfs => 500, :Nmin => 5E3, :Nmax => 2E5),
    Dict(:nPure => 60, :tM => 10, :lVfs => 500, :Nmin => 5E3, :Nmax => 2E5),
    Dict(:nPure => 30, :tM => 10, :lVfs => 500, :Nmin => 5E3, :Nmax => 3E5),
    Dict(:nPure => 10, :tM => 10, :lVfs => 500, :Nmin => 5E3, :Nmax => 5E5),
    Dict(:nPure => 0, :tM => 10, :lVfs => 500, :Nmin => 5E3, :Nmax => 6E5),
]

# current sim params
simNum = parse(Int,ARGS[1])
nPure = paramsInput[simNum][:nPure]
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
    "mature time" => tM,
    "pure births" => _nPure[simNum]
)
# paramsEst = InferencePipeline.estimateRates(SCBurdenHSC_CID)
paramsEst = InferencePipeline.estimateRates(SCBurdenHSC_CID, μKnown)

VafFit = 1
_NDisc = range(Nmin, Nmax, length=40)
_pDisc = range(0.05, 0.99, length=40)

##
@time _N, _p, NOptInterpol_p, vaf1ErrorInterpol_p_N = InferencePipeline.calcNpSpace(paramsKnown, paramsEst, nVHSC_f, SCBurdenHSC_CID, _NDisc, _pDisc, lVfs, VafFit; verbose=true)


filename = "./LSData_NPSpaceInference_tM"*string(tM)*"_lVFS"*string(lVfs)*"_pureGrowth"*string(Npure)*".jld2"
paramsTot = merge(paramsKnown, paramsEst)
@save filename _N _p NOptInterpol_p paramsTot

println("succes")
