## Load packages
# using Gadfly
using Plots, LaTeXStrings, ColorSchemes
using JLD2, FileIO
using Glob
using Loess
include("../../src/vafdyn.jl")
using .VAFDyn
include("../../src/theory.jl")
using .Theory
using Statistics

SAVEPLOTS = true

## =============== Load data ===============
# --------- sim VAF data ---------
# fileNames_ = glob("singlePatientFullSim_Ni50_Nf*NH100*.jld2", "./data/Simulations/Nf1000")
fileNames_ = glob("singlePatientFullSim_Ni1_Nf*.jld2", "./data/Simulations/Nf10000")

nSims = length(fileNames_)
println(nSims)
@load fileNames_[1] paramsTrue timesSim_ nCellSim_t
nV_sim_f = Array{Float64, 2}(undef, nSims, paramsTrue["N final"]+1)
nVS_sim_f = Array{Float64, 2}(undef, nSims, paramsTrue["sample size"]+1)
for (i,fName) in enumerate(fileNames_)
    @load fName nVarSim_f nVarSimS_f
    nV_sim_f[i, :] .= nVarSim_f
    nVS_sim_f[i, :] .= nVarSimS_f
end
# Get averaged vaf spectra
nVAv_f = sum(nV_sim_f, dims=1)/nSims
nVSAv_f = sum(nVS_sim_f, dims=1)/nSims

## --------- pred VAF data ---------
if !haskey(paramsTrue, "pure birth")
    paramsTrue["pure births"] = 0   # must be added as it did not exist in previous version
end
vfs = VAFDyn.VFreqspace(paramsTrue["N final"], 500)
VAFDyn.evolveGrowingVAF(vfs, paramsTrue, paramsTrue["evolve time"])
dfs = VAFDyn.makeDFSfromVFS(vfs, paramsTrue["N final"])
dfsS = VAFDyn.sampler(dfs, paramsTrue["sample size"])


## ======================= Plotting =======================
pyplot()
theme(:default,
    minorgrid=false,
    gridstyle=:dash,
    fontfamily="DejaVu Sans",
    showaxis=true,
    gridlinewidth=0.7,
    # size=(0.9*500,0.9*400),
    size=(500,400),
    legendfontsize=10,
    guidefontsize=12,
    tickfontsize=10,
)
simId = 7

freqs_f = (0:paramsTrue["N final"]) / paramsTrue["N final"]
freqsS_f = (0:paramsTrue["sample size"]) / paramsTrue["sample size"]

## ---------- Sim vs EV - Complete - Log ---------
figSimVFP = plot(yscale=:log10, xscale=:log10)
plot!(freqs_f[2:end], nVAv_f[2:end],
    color=:grey45,
    linealpha=1,
    label="Simulations average",
    # title = "c)",
    # titleloc = :left,
    # titlefont=font(20, "DejaVu Sans"),
)
plot!(freqs_f[2:end], nV_sim_f[simId, 2:end],
    # linewidth=0.3,
    linewidth=0,
    # linetype=:sticks,
    # linetype=:bar,
    # linealpha=0.8,
    # fillalpha=0.8,
    # fillrange=0,
    # fillcolor=:match,
    label="Single simulation",
    markershape=:diamond,
    markersize=2,
    markerstrokewidth=0,
    markeralpha=0.7,
    color=1,
)
# for sim in 2:5
#     plot!(freqs_f[2:end], nV_sim_f[sim, 2:end],
#     # linewidth=0,
#     linewidth=0.2,
#     linealpha=0.5,
#     fillalpha=0.5,
#     fillrange=0,
#     fillcolor=:match,
#     label=""
#     )
# end
plot!(freqs_f[2:end], dfs.n_f[2:end],
    color=:black,
    linestyle=:dash,
    label=L"$v(f,t)$ prediction",
)
plot!(freqs_f[2:end], 10^(-.88)*(nV_sim_f[simId, 2:end].==0),
    markershape=:diamond,
    markersize=1.5,
    markerstrokewidth=0,
    color=1,
    markeralpha=0.7,
    label="",
)
# xlims!(0,0.25)
ylims!(10^(-.9),10^5)
xlims!(1/paramsTrue["N final"],1)
xlabel!(L"Variant allele frequency $f$")
ylabel!(L"Number of variants $v_f$")
display(figSimVFP)

SAVEPLOTS && savefig(figSimVFP, "figures/paper/VafSpectrumSimFPAlt.pdf")

##

## ----------  Sim vs EV - Sample - Log  ---------
# # fig2 = plot(yscale=:log10)
# fig2 = plot(yscale=:log10, xscale=:log10)
# plot!(freqsS_f[2:end], nVS_sim_f[1, 2:end],
# color=2,
# # linetype=:steppre,
# # linetype=:steppost,
# # linewidth=0.2,
# linetype=:sticks,
# linewidth=2,
# linealpha=0.9,
# fillrange=0,
# fillalpha=0.9,
# fillcolor=:match,
# label="single simulation"
# )
# # for sim in 2:5
# #     bar!(freqsS_f[2:end], nVS_sim_f[sim, 2:end],
# #     linewidth=0,
# #     fillalpha=0.5,
# #     label=""
# #     )
# # end
# plot!(freqsS_f[2:end], nVSAv_f[2:end],
# color=:grey35,
# linealpha=1,
# label="simulations average")
# plot!(freqsS_f[2:end], dfsS.n_f[2:end],
# color=:black,
# linestyle=:dash,
# label="predicted average")
# ylims!(10^0,10^5)
# xlims!(1/(1.1paramsTrue["sample size"]),1)
# xlabel!("VAF")
# ylabel!("number of variants")
# title!("Sample")
# display(fig2)

# fig3 = plot(fig1, fig2, layout=2, size=(900, 400))
# display(fig3)
# # SAVEPLOTS && savefig(fig3, "figures/paper/simsVAFSpectrumTrueVsSample.pdf")
# # SAVEPLOTS && savefig(fig3, "figures/paper/simsVAFSpectrumTrueVsSampleLogLog.pdf")


## -------------sampled vs true --------------------
# colorSDark = cgrad(:seaborn_dark, categorical=true)
colorSDark = ColorSchemes.seaborn_dark
colorSLight = ColorSchemes.seaborn_muted
figPopVSample = plot(
    yscale=:log10, xscale=:log10,
    legend=:bottomleft,
)
plot!(freqs_f[2:end], nV_sim_f[1, 2:end],
    color=1,
    # linetype=:sticks,
    # linewidth=0.3,
    linewidth=0,
    # fillrange=0,
    # fillalpha=0.6,
    # linealpha=0.6,
    # fillcolor=:match,
    # label="single simulation",
    label="",
    markershape=:diamond,
    markeralpha=0.6,
    markersize=2.5,
    markerstrokewidth=0,
    title = "c)",
    titleloc = :left,
    titlefont=font(20, "DejaVu Sans"),
)
plot!(freqsS_f[2:end], nVS_sim_f[1, 2:end],
    color=2,
    # linetype=:sticks,
    # linewidth=0.3,
    linewidth=0,
    # fillrange=0,
    # fillalpha=0.6,
    # linealpha=0.6,
    # fillcolor=:match,
    # label="single simulation,"
    label="",
    markershape=:diamond,
    markeralpha=0.6,
    markersize=4,
    markerstrokewidth=0,
)
plot!(freqs_f[2:end], nVAv_f[2:end],
    color=colorSLight[1],
    # label="population:\nsim average",
    # label="population",
    label="",
)
plot!(freqs_f[2:end], dfs.n_f[2:end],
    color=colorSDark[1],
    linestyle=:dash,
    # label="population:\nmodel",
    label="",
)
plot!(freqsS_f[2:end], nVSAv_f[2:end],
    color=colorSLight[2],
    # label="sample:\nsim average",
    # label="sample",
    label="",
)
plot!(freqsS_f[2:end], dfsS.n_f[2:end],
    color=colorSDark[2],
    linestyle=:dash,
    # label="sample:\nmodel",
    label="",
)

# legend labels
plot!(0:0.1:1, (0:0.1:1)*1E-5,
    color=colorSLight[1],
    markershape=:diamond,
    markersize=4,
    markerstrokewidth=0,
    label="Population",
)
plot!(0:0.1:1, (0:0.1:1)*1E-5,
    color=colorSLight[2],
    markershape=:diamond,
    markersize=4,
    markerstrokewidth=0,
    label="Sample",
)
plot!(0:0.1:1, (0:0.1:1)*1E-5,
    color=:grey70,
    # fillrange=1E-6,
    markershape=:diamond,
    markerstrokewidth=0,
    linewidth=0,
    label="Single simulation",
)
plot!(0:0.1:1, (0:0.1:1)*1E-5,
    color=:grey45,
    label="Sim average",
)
plot!(0:0.1:1, (0:0.1:1)*1E-5,
    color=:black,
    linestyle=:dash,
    label=L"$v(f,t)$ prediction",
)
# annotate!((0.001,0.001), Plots.text("a)"))

ylims!(10^-.5,10^5)
xlabel!(L"Variant allele frequency $f$")
ylabel!(L"Number of variants $v_f$")
xlims!(1/paramsTrue["N final"],1)
display(figPopVSample)

SAVEPLOTS && savefig(figPopVSample, "Figures/Paper/3c.pdf")

# fig4 = plot(
#     fig3,
#     xscale=:linear,
#     legend=:topright,
# )
# display(fig4)

## ============== Cumulative VAF ======================

nVCum_sim_f = zeros(size(nV_sim_f))
for simID in 1:size(nV_sim_f,1)
    nVCum_sim_f[simID, 2:end-1] = [sum(@view nV_sim_f[simID, i:end-1]) for i in 2:size(nV_sim_f,2)-1]
end
##
nVCumVar_f = zeros(size(nV_sim_f,2))
nVCumAv_f = zeros(size(nV_sim_f,2))
for fI in 2:(size(nV_sim_f,2)-1)
    nVCumVar_f[fI] = var(nVCum_sim_f[:,fI])
    nVCumAv_f[fI] = mean(nVCum_sim_f[:,fI])
end
##

nVCum_f = [sum(@view nV_sim_f[simId, i:end-1]) for i in 2:size(nV_sim_f,2)-1]
nVAvCum_f = [sum(@view nVAv_f[i:end-1]) for i in 2:length(nVAv_f)-1]
nVExp_f = [sum(dfs.n_f[i:end-1]) for i in 2:length(dfs.n_f)-1]
nVCumStdUp = nVCumAv_f .+ sqrt.(nVCumVar_f)
nVCumStdDw = nVCumAv_f .- sqrt.(nVCumVar_f)
# nVExpMC_f = [sum(dfsMC.n_f[i:end-1]) for i in 2:length(dfsMC.n_f)-1]

## ----------  Sim vs EV - Cumulative - Log  ---------
figCumVaf = plot(
    xscale=:log10,
    yscale=:log10,
    # legend=:topright,
    legend=:bottomleft,
    ylims=(10^-0,10^7)
)
plot!(freqs_f[2:end-1], nVCum_f,
    # markerstrokewidth=0,
    # linetype=:stepmid,
    # fillcolor=:match,
    # fillrange=0,
    # linealpha=0.7,
    # fillalpha=0.8,
    linewidth=0,
    markershape=:hexagon,
    # markershape=:diamond,
    markersize=2,
    markerstrokewidth=0,
    # markeralpha=0.7,
    label="Single simulation",
)
plot!(freqs_f[2:end-1], nVAvCum_f,
    color=:grey45,
    linealpha=1,
    label="Simulations average"
)
# standard deviations
plot!(freqs_f[2:end-1], nVCumStdUp[2:end-1],
    linestyle=:dashdot,
    color=:grey45,
    label="Simulations std"
)
plot!(freqs_f[2:end-1], nVCumStdDw[2:end-1],
    linestyle=:dashdot,
    color=:grey45,
    label=""
)
plot!(freqs_f[2:end-1], nVExp_f,
    color=:black,
    linestyle=:dash,
    label=L"$w(f,t)$ prediction"
)

# title!("growing population (1 -> 10'000)")
# xlims!(0,0.25)
# ylims!(0,200)
# ylims!(10^-0,10^7)
# xlims!(0,0.1)
xlims!(1/paramsTrue["N final"],1)
ylims!(10^(-.5),10^7)
xlabel!(L"Variant allele frequency $f$")
ylabel!(L"Number of variants < $f$")
# xlims!(1/paramsTrue["N final"],1)
display(figCumVaf)

SAVEPLOTS ? savefig(figCumVaf, "Figures/Paper/CumVafSpectrumSimFP.pdf") : 0

## ============= Sampled ================
# 0
# simId = 7
# nVSCum_f = [sum(@view nVS_sim_f[simId, i:end-1]) for i in 2:size(nVS_sim_f,2)-1]
# nVSAvCum_f = [sum(@view nVSAv_f[i:end-1]) for i in 2:length(nVAv_f)-1]
# nVSExp_f = [sum(dfsS.n_f[i:end-1]) for i in 2:length(dfsS.n_f)-1]
# # nVExpMC_f = [sum(dfsMC.n_f[i:end-1]) for i in 2:length(dfsMC.n_f)-1]

# ## ---------- sample ---------
# # fig1 = plot(yscale=:log10)
# fig1 = plot(yscale=:log10)
# # fig1 = plot(yscale=:log10, xscale=:log10, legend=:bottomleft)
# # fig1 = plot()

# for simId in 10:15
#     nVCum_f = [sum(@view nVS_sim_f[simId, i:end-1]) for i in 2:size(nVS_sim_f,2)-1]
#     plot!(freqsS_f[2:end-1], nVSCum_f,
#         markerstrokewidth=0,
#         fillcolor=:match,
#         label="single simulation"
#     )
# end
# plot!(freqsS_f[2:end-1], nVSAvCum_f,
# color=:grey35,
# linealpha=1,
# label="simulations average")
# plot!(freqsS_f[2:end-1], nVSExp_f,
# color=:black,
# linestyle=:dash,
# label="predicted average"
# )
# # plot!(freqs_f[2:end-1], nVExpMC_f,
# # color=:black,
# # linestyle=:dashdot,
# # label="predicted average"
# # )
# title!("growing population (1 -> 10'000)")
# xlims!(0,0.25)
# ylims!(10^-0,10^7)
# xlabel!("Variant allele frequency")
# ylabel!("Number of variants < f")
# xlims!(1/paramsTrue["N final"],1)
# ylims!(1E-0.5, 1E5)
# display(fig1)

# ## ============ Sample v True =============

# fig3 = plot(yscale=:log10)
# fig3 = plot(yscale=:log10, xscale=:log10, legend=:bottomleft)
# # fig3 = plot()

# # for simId in 10:15
# #     nVCum_f = [sum(@view nVS_sim_f[simId, i:end-1]) for i in 2:size(nVS_sim_f,2)-1]
# #     plot!(freqsS_f[2:end-1], nVSCum_f,
# #         markerstrokewidth=0,
# #         fillcolor=:match,
# #         label="single simulation"
# #     )
# # end
# plot!(freqs_f[2:end-1], nVAvCum_f,
# color=1,
# linealpha=1,
# label="true population average")
# # plot!(freqs_f[2:end-1], nVExp_f,
# # color=1,
# # linestyle=:dash,
# # label="predicted average"
# # )

# plot!(freqsS_f[2:end-1], nVSAvCum_f,
# color=2,
# linealpha=1,
# label="sampled population average")
# # plot!(freqsS_f[2:end-1], nVSExp_f,
# # color=2,
# # linestyle=:dash,
# # label="predicted average"
# # )

# # plot!(freqs_f[2:end-1], nVExpMC_f,
# # color=:black,
# # linestyle=:dashdot,
# # label="predicted average"
# # )
# title!("growing population (1 -> 10'000)")
# xlims!(0,0.25)
# ylims!(10^-0,10^7)
# xlabel!("Variant allele frequency")
# ylabel!("Number of variants < f")
# xlims!(1/paramsTrue["N final"],1)
# display(fig3)

# SAVEPLOTS ? figname = "cumulativeVAF_TrueVSample.pdf" : 0
# savefig(fig3, figname)
# ##

