using JLD2
include("../src/inferencePipeline.jl")
using .InferencePipeline
include("../src/theory.jl")
using .Theory
using ProgressMeter


paramsInput = Dict(
    :nPure => 500,
    :tM => 10,
    :lVfs => 500,
    :N0 => 1,
    :Nmin => 5E3,
    :Nmax => 1.5E5
)

# load data

μKnown = 1.2
paramsKnown = Dict{String, Real}(
    "evolve time" => 59,
    "N initial" => paramsInput[:N0],
    "sample size" => sampleSize,
    "mature time" => paramsInput[:tM],
    "pure births" => paramsInput[:nPure]
)

@load "HPC/LSDataStatsBM.jld2" sampleSize SCBurdenHSC_CID nVHSC_f
paramsEst = InferencePipeline.estimateRates(SCBurdenHSC_CID, μKnown)

## Find optimal N-p pair

vafFit = 1
_NDisc = range(Nmin, Nmax, length=8)
p = 0.5

function getPopSize(paramsKnown, paramsEst, _N, p, lVfs::Int, vafFit::Int, nVafDataFit::Int)
    paramsIn = merge(paramsKnown, paramsEst)
    # nVafTheoryFit_p_N = calcExpectedVAFs(paramsIn, _N, _p, lVFS, vafFit)
    nVafFit_N = Vector{Float64}(undef, length(_N))
    @showprogress for (j,N) in enumerate(_N)
        paramsParticle = InferencePipeline.createTestParams(paramsIn, N, p)
        vfs = VAFDyn.VFreqspace(paramsParticle["N final"], lVfs)
        VAFDyn.evolveGrowingVAF(vfs, paramsParticle, paramsParticle["evolve time"])
        dfs = VAFDyn.makeDFSfromVFS(vfs, paramsParticle["N final"])
        dfsS = VAFDyn.sampler(dfs, paramsParticle["sample size"])
        nVafFit_N[j] = dfsS.n_f[1+vafFit]
    end
    InferencePipeline.findZeroErrorPoint(nVafFit_N, nVafDataFit, _N)
end

##
Nopt = getPopSize(paramsKnown, paramsEst, _N, p, 100, vafFit, nVHSC_f[1+vafFit])

##
paramsDef = InferencePipeline.createTestParams(merge(paramsEst, paramsKnown), Nopt, p)

vfs = VAFDyn.VFreqspace(paramsDef["N final"], 500)
VAFDyn.evolveGrowingVAF(vfs, paramsDef, paramsDef["evolve time"])
dfs = VAFDyn.makeDFSfromVFS(vfs, paramsDef["N final"])
dfsS = VAFDyn.sampler(dfs, paramsDef["sample size"])

## cumulative
nVTheoryCum_f = [sum(dfsS.n_f[k:end]) for k in 1:length(dfsS.n_f)]
nVHSCCum_f = [sum(nVHSC_f[k:end]) for k in 1:length(nVHSC_f)]

## plot

fig1 = plot(dfsS.freqs_f[2:end], nVHSCCum_f[2:end], 
    xaxis=:log10,
    yaxis=:log10, 
    color=:grey55,
    markershape=:circle,
    markersize=4,
    markerstrokecolor=:white,
    linewidth=1.5,
    ylims=(1E0, 8E4),
    label="HSC data",
)
plot!(dfsS.freqs_f[2:end], nVTheoryCum_f[2:end])

plot!(dfsS.freqs_f[2:end], 24 .*(1 ./dfsS.freqs_f[2:end] .- 1) )

xlabel!(L"Variant alelle frequency $f$")
ylabel!(L"Number of variants $< f$")
ylims!(10^-0.3, 10^5)

display(fig1)


##
fig2 = bar(freqs_f, nVHSC_f, yaxis=:log10,
    label="HSC data",
    # linestyle=:sticks,
    color=:grey70,
    linewidth=0)

plot!(dfsS.freqs_f[2:end], dfsS.n_f[2:end], 
    # label="pure birhts = ",
    linestyle=:dash,
    # color=2,
    # linewidth=3
)

xlims!(0,0.6)
ylims!(10^-0.3, 10^5)
xlabel!(L"Variant allele frequency $f$")
ylabel!("Number of variants")
# title!("μ = "*μString)
display(fig2)
