using JLD2, LaTeXStrings

## ============== Load Data ==============

pureGrowth=false
NOptInterpol_tM_P = Vector{Float64}[]
lVfs_ = Int[]
tM_ = Int[]

dataFolder = "data/LSData_NPSpace_21-06-15/"
let
    if pureGrowth
        filename = "LSData_NPSpaceInference_tM1_lVFS500_pureGrowthtrue.jld2"
    else
        # filename = "LSData_NPSpaceInference_tM1_lVFS500_pureGrowthfalse.jld2"
        filename = "LSData_NPSpaceInference_tM1_lVFS500_pureGrowth50.jld2"
    end
    @load dataFolder*filename _p NOptInterpol_p paramsTot
    push!(NOptInterpol_tM_P, NOptInterpol_p)
    push!(tM_, paramsTot["mature time"])
    global _p = _p
end
let
    if pureGrowth
        filename = "LSData_NPSpaceInference_tM3_lVFS500_pureGrowthtrue.jld2"
    else
        filename = "LSData_NPSpaceInference_tM3_lVFS500_pureGrowth50.jld2"
    end
    @load dataFolder*filename NOptInterpol_p paramsTot
    push!(NOptInterpol_tM_P, NOptInterpol_p)
    push!(tM_, paramsTot["mature time"])
end
let
    if pureGrowth
        filename = "LSData_NPSpaceInference_tM5_lVFS500_pureGrowthtrue.jld2"
    else
        filename = "LSData_NPSpaceInference_tM6_lVFS500_pureGrowth50.jld2"
    end
    @load dataFolder*filename NOptInterpol_p paramsTot
    push!(NOptInterpol_tM_P, NOptInterpol_p)
    push!(tM_, paramsTot["mature time"])
end
let
    if pureGrowth
        filename = "LSData_NPSpaceInference_tM8_lVFS500_pureGrowthtrue.jld2"
    else
        filename = "LSData_NPSpaceInference_tM10_lVFS500_pureGrowth50.jld2"
    end
    @load dataFolder*filename NOptInterpol_p paramsTot
    push!(NOptInterpol_tM_P, NOptInterpol_p)
    push!(tM_, paramsTot["mature time"])
end


## =========== Plotting ===============
using Plots
pyplot()
theme(:default,
    minorgrid=false,
    gridstyle=:dash,
    fontfamily="DejaVu Sans",
    showaxis=true,
    gridlinewidth=0.7,
    size=(0.9*500,0.9*400),
)



## --------- NP-space ------------
palette = cgrad(:tableau_red_blue, 4, categorical=true, rev=true)

fig1 = plot(size=(0.9*600,0.9*400))
for i in 1:length(NOptInterpol_tM_P)
    plot!(
        _p, NOptInterpol_tM_P[i], 
        label="Maturation time: "*string(tM_[i]), 
        # marker=:circle,
        color=palette[i],
    )
end
xlabel!(L"asymmetric division fraction $p$")
ylabel!(L"mature population size $N$")
xlims!(0, 1)
# ylims!(2E4, 4E5)
ylims!(0.5E4, 3E5)

display(fig1)
savefig(fig1, "Figures/LSData/NPSpace_mu1.2.pdf")


## --------- SC Burden -----------

@load "HPC/LSDataStatsBM.jld2" sampleSize SCBurdenHSC_CID nVHSC_f
paramsEst = InferencePipeline.estimateRates(SCBurdenHSC_CID)
cpVals_id = CompoundPoisson.randComPois(paramsEst["divisions"], paramsEst["μ"], Int(1E6))

##
fig2 = histogram(
    SCBurdenHSC_CID, normalize=true, bins=25,
    label="HSC data",
    color=:grey60,
    linecolor=:grey60,
    # linealpha=:0,
)
stephist!(
    cpVals_id, normalize=true,
    label="Compound Poisson\ndistribution",
    linestyle=:dash,
    color=:black,
)
xlabel!(L"Number of mutations $m$")
ylabel!("Density of cells")
display(fig2)

# figname = "Figures/LSData_Fitting/mutationalBurden.pdf"
# savefig(fig2, figname)

## --------- VAF spectrum -----------
include("../../src/inferencePipeline.jl")
using .InferencePipeline
include("../../src/vafdyn.jl")
using .VAFDyn

##

dfsFit1_tM = []
dfsFit2_tM = []
for tM in 1:length(NOptInterpol_tM_P)
    p = _p[40]
    N = NOptInterpol_tM_P[tM][40]

    @load "HPC/LSDataStatsBM.jld2" sampleSize SCBurdenHSC_CID nVHSC_f

    μKnown = 1.2
    paramsKnown = Dict{String, Real}(
        "evolve time" => 59,
        "N initial" => 1,
        "sample size" => sampleSize,
        "mature time" => tM
    )
    paramsEst = InferencePipeline.estimateRates(SCBurdenHSC_CID, μKnown)
    paramsComb = InferencePipeline.createTestParams(merge(paramsEst, paramsKnown), N, p; pureGrowth)
    println(paramsComb)

    vfsFit = VAFDyn.VFreqspace(Int(round(N)),500)
    @time VAFDyn.evolveGrowingVAF(vfsFit, paramsComb, paramsComb["evolve time"], pureGrowth=50)
    # @time VAFDyn.evolveGrowingVAFpureGrowth(vfsFit, paramsComb, paramsComb["evolve time"], tW=100)
    dfsFit = VAFDyn.makeDFSfromVFS(vfsFit, paramsComb["N final"])
    @time dfsFitSampled = VAFDyn.sampler(dfsFit, paramsComb["sample size"])
    push!(dfsFit1_tM, dfsFitSampled)
end
##

## ==== Plot fitted and measured VAFs ====
# μString = "1.2"
# pString = string(userParams["p"])
# NString = string(userParams["N"])
@load "HPC/LSDataStatsBM.jld2" nVHSC_f freqs_f

##

fig2 = bar(freqs_f, nVHSC_f, yaxis=:log10,
    label="HSC data",
    # linestyle=:sticks,
    color=:grey70,
    linewidth=0)

for (i,tM) in enumerate(tM_)
    plot!(dfsFit1_tM[i].freqs_f[2:end], dfsFit1_tM[i].n_f[2:end], 
        label="maturation time = "*string(tM),       
        linestyle=:dash,
        # color=2,
        # linewidth=3
        )
end

xlims!(0,0.6)
ylims!(10^-0.3, 10^5)
xlabel!(L"Variant allele frequency $f$")
ylabel!("Number of variants")
# title!("μ = "*μString)
display(fig2)

figname = "Figures/LSData_Fitting/VAFSpectrum_density.pdf"
savefig(fig2, figname)

## ----- cumulative -----
dfsCum_tM_F = [
    [
        sum(dfs.n_f[k:end]) for k in 1:length(dfs.n_f)
    ] 
    for dfs in dfsFit1_tM
]
nVHSCCum_f = [sum(nVHSC_f[k:end]) for k in 1:length(nVHSC_f)]
##

fig3 = plot(freqs_f[2:end], nVHSCCum_f[2:end], 
    # xaxis=:log10,
    yaxis=:log10, 
    color=:grey55,
    markershape=:circle,
    markersize=4,
    markerstrokecolor=:white,
    linewidth=1.5,
    ylims=(1E0, 8E4),
    label="HSC data",
)
for (i,tM) in enumerate(tM_)
    plot!(freqs_f[2:end], dfsCum_tM_F[i][2:end],
        label="maturation time = "*string(tM),
    )
end
xlabel!(L"Variant alelle frequency $f$")
ylabel!(L"Number of variants $< f$")
ylims!(10^-0.3, 10^5)

display(fig3)
figname = "Figures/LSData_Fitting/VAFSpectrum_cumulative.pdf"
savefig(fig3, figname)