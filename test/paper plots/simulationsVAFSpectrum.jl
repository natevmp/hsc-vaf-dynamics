## Load packages
# using Gadfly
using Plots
using JLD2, FileIO
using Glob
using Loess
include("../../src/vafdyn.jl")
using .VAFDyn
include("../../src/theory.jl")
using .Theory

SAVEPLOTS = false

## =============== Load data ===============
# --------- sim VAF data ---------

# fileNames_ = glob("singlePatientFullSim_Ni10000_Nf*.jld2", "./data/Simulations/Nf10000")
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
# filenameVAFDyn = "data/Simulations/Nf10000/vafDynGrowth_N"*string(paramsTrue["N final"])*".jld2"
# @load filenameVAFDyn dfs dfsS dfsPDE dfsPDES

# dfsMC = VAFDyn.DFreqspace(paramsTrue["N final"])
# VAFDyn.evolveGrowingVAF(dfsMC, paramsTrue, paramsTrue["evolve time"])

##
vfs = VAFDyn.VFreqspace(paramsTrue["N final"], 500)
VAFDyn.evolveGrowingVAF(vfs, paramsTrue, paramsTrue["evolve time"])
dfs = VAFDyn.makeDFSfromVFS(vfs, paramsTrue["N final"])
dfsS = VAFDyn.sampler(dfs, paramsTrue["sample size"])




## ======================= Plotting =======================
pyplot(legendfontsize=10, guidefontsize=12, tickfontsize=10,  size=(500,400))
simTestID = 5

freqs_f = (0:paramsTrue["N final"]) / paramsTrue["N final"]
freqsS_f = (0:paramsTrue["sample size"]) / paramsTrue["sample size"]

## ---------- complete ---------
fig1 = plot(yscale=:log10)
# fig1 = plot(yscale=:log10, xscale=:log10)
plot!(freqs_f[2:end], nV_sim_f[simTestID, 2:end],
# linewidth=0,
# linetype=:steppre,
# linetype=:steppost,
# linewidth=0.2,
linestyle=:sticks,
linewidth=0.5,
linealpha=0.9,
fillrange=0,
fillalpha=0.9,
fillcolor=:match,
label="single simulation"
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
# model = loess(freqs_f[2:end], nVAv_f[2:end])
# vs = predict(model, freqs_f[2:end])
# plot!(freqs_f[2:end], vs,
plot!(freqs_f[2:end], nVAv_f[2:end],
color=:grey35,
linealpha=1,
label="simulations average")
plot!(freqs_f[2:end], dfs.n_f[2:end],
color=:black,
linestyle=:dash,
label="predicted average"
)
# xlims!(0,0.25)
ylims!(10^0,10^5)
xlabel!("VAF")
ylabel!("number of variants")
xlims!(1/paramsTrue["N final"],1)
title!("complete")
display(fig1)

## ---------- sample ---------
# fig2 = plot(yscale=:log10)
fig2 = plot(yscale=:log10, xscale=:log10)
plot!(freqsS_f[2:end], nVS_sim_f[1, 2:end],
color=2,
# linetype=:steppre,
# linetype=:steppost,
# linewidth=0.2,
linetype=:sticks,
linewidth=2,
linealpha=0.9,
fillrange=0,
fillalpha=0.9,
fillcolor=:match,
label="single simulation"
)
# for sim in 2:5
#     bar!(freqsS_f[2:end], nVS_sim_f[sim, 2:end],
#     linewidth=0,
#     fillalpha=0.5,
#     label=""
#     )
# end
plot!(freqsS_f[2:end], nVSAv_f[2:end],
color=:grey35,
linealpha=1,
label="simulations average")
plot!(freqsS_f[2:end], dfsS.n_f[2:end],
color=:black,
linestyle=:dash,
label="predicted average")
ylims!(10^0,10^5)
xlims!(1/(1.1paramsTrue["sample size"]),1)
xlabel!("VAF")
ylabel!("number of variants")
title!("sample")
display(fig2)

## ------- save fig --------
fig3 = plot(fig1, fig2, layout=2, size=(900, 400))
display(fig3)
# SAVEPLOTS && savefig(fig3, "figures/paper/simsVAFSpectrumTrueVsSample.pdf")
# SAVEPLOTS && savefig(fig3, "figures/paper/simsVAFSpectrumTrueVsSampleLogLog.pdf")



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

simId = 7
nVCum_f = [sum(@view nV_sim_f[simId, i:end-1]) for i in 2:size(nV_sim_f,2)-1]
nVAvCum_f = [sum(@view nVAv_f[i:end-1]) for i in 2:length(nVAv_f)-1]
nVExp_f = [sum(dfs.n_f[i:end-1]) for i in 2:length(dfs.n_f)-1]
nVCumStdUp = nVCumAv_f .+ sqrt.(nVCumVar_f)
nVCumStdDw = nVCumAv_f .- sqrt.(nVCumVar_f)
# nVExpMC_f = [sum(dfsMC.n_f[i:end-1]) for i in 2:length(dfsMC.n_f)-1]

## ---------- complete ---------
# fig1 = plot(yscale=:log10)
# fig1 = plot(yscale=:log10)
fig1 = plot(yscale=:log10, legend=:topright,ylims=(10^-0,10^7))
# fig1 = plot(legend=:topright)

# for simId in 10:15
#     nVCum_f = [sum(@view nV_sim_f[simId, i:end-1]) for i in 2:size(nV_sim_f,2)-1]
#     plot!(freqs_f[2:end-1], nVCum_f,
#         markerstrokewidth=0,
#         fillcolor=:match,
#         label="single simulation"
#     )
# end
plot!(freqs_f[2:end-1], nVAvCum_f,
color=:grey35,
linealpha=1,
label="simulations average")
plot!(freqs_f[2:end-1], nVExp_f,
color=:black,
linestyle=:dash,
label="predicted average"
)
plot!(freqs_f[2:end-1], nVCumStdUp[2:end-1],
linestyle=:dash,
label="standard deviation"
)
plot!(freqs_f[2:end-1], nVCumStdUp[2:end-1],
linestyle=:dashdot,
color=:grey55,
label="standard deviation"
)
plot!(freqs_f[2:end-1], nVCumStdDw[2:end-1],
linestyle=:dashdot,
color=:grey55,
# label="standard deviation"
)
# plot!(freqs_f[2:end-1], nVExpMC_f,
# color=:black,
# linestyle=:dashdot,
# label="predicted average"
# )
title!("growing population (1 -> 10'000)")
# xlims!(0,0.25)
# ylims!(0,200)
# ylims!(10^-0,10^7)
# xlims!(0,0.1)
xlabel!("Variant allele frequency")
ylabel!("Number of variants < f")
# xlims!(1/paramsTrue["N final"],1)
display(fig1)

figname = "cumulativeVAF_growingN_linlin.pdf"
savefig(fig1, figname)

## ============= Sampled ================

simId = 7
nVSCum_f = [sum(@view nVS_sim_f[simId, i:end-1]) for i in 2:size(nVS_sim_f,2)-1]
nVSAvCum_f = [sum(@view nVSAv_f[i:end-1]) for i in 2:length(nVAv_f)-1]
nVSExp_f = [sum(dfsS.n_f[i:end-1]) for i in 2:length(dfsS.n_f)-1]
# nVExpMC_f = [sum(dfsMC.n_f[i:end-1]) for i in 2:length(dfsMC.n_f)-1]

## ---------- sample ---------
# fig1 = plot(yscale=:log10)
fig1 = plot(yscale=:log10)
# fig1 = plot(yscale=:log10, xscale=:log10, legend=:bottomleft)
# fig1 = plot()

for simId in 10:15
    nVCum_f = [sum(@view nVS_sim_f[simId, i:end-1]) for i in 2:size(nVS_sim_f,2)-1]
    plot!(freqsS_f[2:end-1], nVSCum_f,
        markerstrokewidth=0,
        fillcolor=:match,
        label="single simulation"
    )
end
plot!(freqsS_f[2:end-1], nVSAvCum_f,
color=:grey35,
linealpha=1,
label="simulations average")
plot!(freqsS_f[2:end-1], nVSExp_f,
color=:black,
linestyle=:dash,
label="predicted average"
)
# plot!(freqs_f[2:end-1], nVExpMC_f,
# color=:black,
# linestyle=:dashdot,
# label="predicted average"
# )
title!("growing population (1 -> 10'000)")
xlims!(0,0.25)
ylims!(10^-0,10^7)
xlabel!("Variant allele frequency")
ylabel!("Number of variants < f")
xlims!(1/paramsTrue["N final"],1)
display(fig1)

## ============ Sample v True =============

fig3 = plot(yscale=:log10)
fig3 = plot(yscale=:log10, xscale=:log10, legend=:bottomleft)
# fig3 = plot()

# for simId in 10:15
#     nVCum_f = [sum(@view nVS_sim_f[simId, i:end-1]) for i in 2:size(nVS_sim_f,2)-1]
#     plot!(freqsS_f[2:end-1], nVSCum_f,
#         markerstrokewidth=0,
#         fillcolor=:match,
#         label="single simulation"
#     )
# end
plot!(freqs_f[2:end-1], nVAvCum_f,
color=1,
linealpha=1,
label="true population average")
# plot!(freqs_f[2:end-1], nVExp_f,
# color=1,
# linestyle=:dash,
# label="predicted average"
# )

plot!(freqsS_f[2:end-1], nVSAvCum_f,
color=2,
linealpha=1,
label="sampled population average")
# plot!(freqsS_f[2:end-1], nVSExp_f,
# color=2,
# linestyle=:dash,
# label="predicted average"
# )

# plot!(freqs_f[2:end-1], nVExpMC_f,
# color=:black,
# linestyle=:dashdot,
# label="predicted average"
# )
title!("growing population (1 -> 10'000)")
xlims!(0,0.25)
ylims!(10^-0,10^7)
xlabel!("Variant allele frequency")
ylabel!("Number of variants < f")
xlims!(1/paramsTrue["N final"],1)
display(fig3)

figname = "cumulativeVAF_TrueVSample.pdf"
savefig(fig3, figname)
##

