using JLD2, LaTeXStrings, Glob, FileIO
include("../../src/inferencePipeline.jl")
using .InferencePipeline
include("../../src/compoundPoisson.jl")
using .CompoundPoisson
include("../../src/vafdyn.jl")
using .VAFDyn

SAVEFIGS = true
## ============== Load Data ==============

NOptInterpol_nP_P = Vector{Float64}[]
vaf1Er_nP_P_n = Array{Float64,2}[]
lVfs_ = Int[]
# _tM = Int[]
_nP = Int[]

dataFolder = "data/LSData_NPSpace_21-08-05/"
filenames_ = glob(dataFolder*"LSData_NPSpaceInference_tM5_lVFS500_pureGrowth*")

for filename in filenames_

    @load filename NOptInterpol_p paramsTot _p _N vaf1ErrorInterpol_p_N
    global _p = _p
    global _N = _N
    push!(NOptInterpol_nP_P, NOptInterpol_p)
    println(paramsTot)
    # push!(_tM, paramsTot["mature time"])
    push!(_nP, paramsTot["pure births"])
    push!(vaf1Er_nP_P_n, vaf1ErrorInterpol_p_N)

end

@load "HPC/LSDataStatsBM.jld2" nVHSC_f freqs_f sampleSize SCBurdenHSC_CID


## ============== Calculate VAF spectra =============

dfsFit1_nP = VAFDyn.DFreqspace[]
paramsFit = Dict{String, Real}
for nP in 1:length(NOptInterpol_nP_P)
    p = _p[40]
    N = NOptInterpol_nP_P[nP][40]

    S, SCBurden_CID = load("HPC/LSDataStatsBM.jld2", "sampleSize", "SCBurdenHSC_CID")
    μKnown = 1.2
    global paramsKnown = Dict{String, Real}(
        "evolve time" => 59,
        "N initial" => 1,
        # "N initial" => 4,
        "sample size" => S,
        "mature time" => 5,
        "pure births" => nP,
    )
    paramsEst = InferencePipeline.estimateRates(SCBurden_CID, μKnown)
    paramsComb = InferencePipeline.createTestParams(merge(paramsEst, paramsKnown), N, p)
    paramsFit = paramsComb
    println(paramsComb)
    println(paramsComb["ρ"] + paramsComb["ϕ"])

    vfsFit = VAFDyn.VFreqspace(Int(round(N)),500)
    @time VAFDyn.evolveGrowingVAF(vfsFit, paramsComb, paramsComb["evolve time"])
    # @time VAFDyn.evolveGrowingVAFpureGrowth(vfsFit, paramsComb, paramsComb["evolve time"], tW=100)
    dfsFit = VAFDyn.makeDFSfromVFS(vfsFit, paramsComb["N final"])
    @time dfsFitSampled = VAFDyn.sampler(dfsFit, paramsComb["sample size"])
    push!(dfsFit1_nP, dfsFitSampled)
end

## -------- Calculate Compound Poisson Distribution -----------
@load "HPC/LSDataStatsBM.jld2" sampleSize SCBurdenHSC_CID nVHSC_f
paramsEst = InferencePipeline.estimateRates(SCBurdenHSC_CID)
cpVals_id = CompoundPoisson.randComPois(paramsEst["divisions"], paramsEst["μ"], Int(1E6))

## -------- get Regular Poisson distribution -------------
using Distributions

burdenDistPoisson = Poisson(mean(SCBurdenHSC_CID))




## ===================== Plotting ============================
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


## --------- SC Burden -----------
fig2 = histogram(
    SCBurdenHSC_CID, normalize=true, bins=25,
    label="HSC data",
    color=:grey60,
    linecolor=:grey60,
    # linealpha=:0,
    title = "a)",
    titleloc = :left,
    titlefont=font(20, "DejaVu Sans"),
)
plot!(
    700:1400, (x -> pdf(burdenDistPoisson, x)).(700:1400),
    color=:grey35,
    linestyle=:solid,
    label="Poisson distribution",
)
stephist!(
    cpVals_id, normalize=true, bins=300,
    label="Comp. Poisson\ndistribution",
    linestyle=:dash,
    color=:black,
)
xlabel!(L"Number of mutations $m$")
ylabel!("Density of cells")
ylims!(0, 0.008)
display(fig2)
# figname = "Figures/LSData_Fitting/mutationalBurden.pdf"
figname="Figures/Paper/4a.pdf"
SAVEFIGS && savefig(fig2, figname)

## --------- VAF spectrum -----------
fig2 = plot()
plot!(freqs_f, nVHSC_f, yaxis=:log10,
    label="HSC data",
    # linetype=:bar,
    # linetype=:stepmid,
    # fillrange=0,
    markershape=:diamond,
    markersize=4,
    color=:grey70,
    linewidth=0,
    title = "b)",
    titleloc = :left,
    titlefont=font(20, "DejaVu Sans"),
    )
# for (i,nP) in enumerate(_nP)
#     plot!(dfsFit1_nP[i].freqs_f[2:end], dfsFit1_nP[i].n_f[2:end], 
#         label="pure birhts = "*string(nP),       
#         linestyle=:dash,
#         # color=2,
#         # linewidth=3
#         )
# end
plot!(dfsFit1_nP[1].freqs_f[2:end], dfsFit1_nP[1].n_f[2:end], 
    label=L"$v(f,t)$ fit",       
    linestyle=:dash,
    color=:black,
    # linewidth=3
    )
xlims!(0,0.5)
ylims!(10^-0.3, 10^5)
xlabel!(L"Variant allele frequency $f$")
ylabel!(L"Number of variants $v_f$")
# title!("μ = "*μString)
display(fig2)
# figname = "Figures/LSData_Fitting/VAFSpectrum_density.pdf"
figname = "Figures/Paper/4b.pdf"
SAVEFIGS && savefig(fig2, figname)


## --------- NP-space ------------
# plotInds = [1,2,5,6,7]
plotInds = 1:length(_nP)

palette = cgrad(:tableau_red_blue, length(plotInds), categorical=true, rev=true)


fig1 = plot(size=(0.9*500,0.9*350))
for i in plotInds
    plot!(
        _p, NOptInterpol_nP_P[i], 
        label=latexstring("N_H = "*string(_nP[i])), 
        # marker=:circle,
        color=palette[i],
    )
end
xlabel!(L"Asymmetric division fraction $p$")
ylabel!(L"Mature population size $N$")
xlims!(0, 1)
# ylims!(2E4, 4E5)
ylims!(0.5E4, 4E5)

display(fig1)
# SAVEFIGS && savefig(fig1, "Figures/LSData/NPSpace_mu1.2.pdf")

## --------- NP-space with error -----------
pyplot()
# gr()
theme(:default,
    minorgrid=false,
    gridstyle=:dash,
    fontfamily="DejaVu Sans",
    showaxis=true,
    gridlinewidth=0.7,
    size=(0.9*500,0.9*400),
)

##

colorgrad = cgrad(:grayC, rev = false, alpha = nothing, scale = nothing, categorical = nothing)
colorgrad=:vik
figNP = plot(
    size=(0.9*500,0.9*350),
    title = "c)",
    titleloc = :left,
    titlefont=font(20, "DejaVu Sans"),
)
heatmap!(
    _p, _N, transpose(vaf1Er_nP_P_n[2]),
    clims=(-.3,.3),
    c=colorgrad,
    colorbar_title=L"Relative error $ \Delta v_{1/S} $",
)
# contour!(_p, _N, transpose(vaf1Er_nP_P_n[1]), clims=(-.3,.3), fill=true, c=colorgrad)
plot!(
    _p, NOptInterpol_nP_P[2], 
    label="Optimal fit",
    # marker=:circle,
    color=:black,
)
xlabel!(L"Asymmetric divisions fraction $p$")
ylabel!(L"Population size at maturity $N$")

# ylims!(0,3E5)
display(figNP)
SAVEFIGS && savefig(figNP, "Figures/Paper/4c.pdf")


## ----- cumulative -----
dfsCum_nP_F = [
    [
        sum(dfs.n_f[k:end]) for k in 1:length(dfs.n_f)
    ] 
    for dfs in dfsFit1_nP
]
nVHSCCum_f = [sum(nVHSC_f[k:end]) for k in 1:length(nVHSC_f)]
##

fig3 = plot(freqs_f[2:end], nVHSCCum_f[2:end], 
    xaxis=:log10,
    yaxis=:log10, 
    color=:grey55,
    markershape=:circle,
    markersize=4,
    markerstrokecolor=:white,
    linewidth=1.5,
    # ylims=(1E0, 8E4),
    label="HSC data",
)
for (i,nP) in enumerate(_nP)
    plot!(freqs_f[2:end], dfsCum_nP_F[i][2:end],
        label="pure births = "*string(nP),
    )
end
xlabel!(L"Variant alelle frequency $f$")
ylabel!(L"Number of variants $< f$")
xlims!(1/sampleSize, 1)
ylims!(10^-0.3, 10^5)
display(fig3)
figname = "Figures/LSData_Fitting/VAFSpectrum_cumulative.pdf"
SAVEFIGS && savefig(fig3, figname)