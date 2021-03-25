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


_NDisc = [1E3,1E4,6E4,9E4]
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

VafFit = 3

##
@time _N, _p, NOptInterpol_p, vaf1ErrorInterpol_p_N = InferencePipeline.calcNpSpace(paramsKnown, nVHSC_f, SCBurdenHSC_CID, _NDisc, _pDisc, 500, VafFit; verbose=true)


##

fig1 = plot(_p, NOptInterpol_p, label="vaf fit:"*string(VafFit))
xlabel!("p")
ylabel!("N")

display(fig1)
savefig(fig1, "NP_fit"*string(VafFit)*".pdf")
