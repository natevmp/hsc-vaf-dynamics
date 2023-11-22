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
using Statistics
##
μ = 1.2
eM = mean(SCBurdenHSC_CID)
vM = var(SCBurdenHSC_CID)

y = (vM - eM)/μ^2
l = eM - y*μ

λ(u) = (eM-(vM-eM)/u)/59
fig1 = plot(1:.1:20, λ.(1:.1:20),label="")
xlabel!("μ")
ylabel!("λ (years)")
ylims!(-10,20)
savefig(fig1, "figMutBurden.png")
display(fig1)

##
vY = (vM-eM)/ μ^2 / 59
λ2(γ) = (eM - γ*59 * μ)/59
fig2 = plot(0:.1:20, λ2.(0:.1:20),label="")
xlabel!("division rate γ = E(y)/t (years)")
ylabel!("external mutation rate λ (years)")
ylims!(-0,20)
annotate!([(15, 15, ("Var(y)/t = "*string(round((vM-eM)/ μ^2 / 59, digits=1)), 12, :center))])
# annotate!([(5, y[5], ("this is #5", 16, :red, :center)), (10, y[10], ("this is #10", :right, 20, "courier"))])
savefig(fig2, "figMutBurden2.png")
display(fig2)
# println((vM-eM)/ μ^2 / 59)
##
paramsEst = InferencePipeline.estimateRates(SCBurdenHSC_CID)
cpVals_id = CompoundPoisson.randComPois(paramsEst1["divisions"], paramsEst1["μ"], 50000)
##

figTest = histogram(SCBurdenHSC_CID, normalize=true, bins=25)
stephist!(cpVals1_id, normalize=true)
stephist!(cpVals2_id, normalize=true)
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