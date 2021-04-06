##
using Plots
pyplot()
# gr()
# plotlyjs()
# plotly()
##
using JLD2
using Statistics
using Distributions
using ProgressMeter
using Dierckx, Roots, Interpolations
using LaTeXStrings
##
include("../../src/vafdyn.jl")
using .VAFDyn
include("../../src/vafSim.jl")
using .VAFSim
include("../../src/theory.jl")
using .Theory
include("../../src/inferencePipeline.jl")
using .InferencePipeline
##
SAVEPLOTS=false
LOADDATA=true

## ================================== get simulation data ==================================

if LOADDATA
    # @load "./data/vafSim_50sims_Ni1_Nf2000.jld2"
    @load "data/SinglePatientPipeline/Nf10000/singlePatientFullSim_Nf10000_sim10.jld2"
    # @load "data/SinglePatientPipeline/singlePatientFullSim_Nf1000.jld2"
    N_N = 500:10000:80000
    p_p = 0.1:0.2:0.9
else

    # User params
    paramsTrue = Dict{String,Real}(
        "N initial" => 1,
        "N final" => 500,
        "μ" => 2.,
        "λ" => 5.,
        "p" => 0.5,
        "sample size" => 100,
        "mature time" => 15,
        "evolve time" => 70
    )
    Theory.extendParams!(paramsTrue)

    # Single Patient simulation
    @time timesSim_t, nCellSim_t, nVarSim_f, nVarSimS_f, nCellSim_m, nCellSimS_m, mSim_cid, mSimS_cid = VAFSim.birthDeathFixedGrowth(paramsTrue, paramsTrue["evolve time"], 0.1, showprogress=true)
    # filename = "data/SinglePatientPipeline/singlePatientFullSim_Nf"*string(paramsTrue["N final"])*".jld2"
    # @save filename paramsTrue timesSim_t nCellSim_t nVarSim_f nVarSimS_f nCellSim_m nCellSimS_m mSim_cid mSimS_cid

    N_N = 100:500:3000
    p_p = 0.1:0.2:0.9
end

## ============================== create VAF spectra ==============================

paramsKnown = Dict{String,Real}(
    "N initial" => 1,
    "sample size" => paramsTrue["sample size"],
    "mature time" => paramsTrue["mature time"],
    "evolve time" => paramsTrue["evolve time"]
)

paramsKnown["mature time"] = 3

fFit_ = [1,2,3,4,5]
# fFit_ = [1,]

NOptInterpol_FFit_p = Vector{Float64}[]
for fFit in fFit_
    _N, _p, NOptInterpol_p, vaf1ErrorInterpol_p_N = InferencePipeline.calcNpSpace(paramsKnown, nVarSimS_f, mSimS_cid, N_N, p_p, 500, fFit, verbose=true)
    push!(NOptInterpol_FFit_p, NOptInterpol_p)
    global _p = _p
end

## ==================================== Plotting =======================================

fig1 = plot()
for (i,fFit) in enumerate(fFit_)
    plot!(_p, NOptInterpol_FFit_p[i], label="vaf fit:"*string(fFit), linestyle=:dash)
end
xlabel!("p")
ylabel!("N")

display(fig1)


##
# dfsTrue = VAFDyn.DFreqspace(paramsTrue["N final"])
# @time VAFDyn.evolveGrowingVAF(dfsTrue, paramsTrue, paramsTrue["evolve time"])
# dfsTrueS = VAFDyn.sampler(dfsTrue, paramsTrue["sample size"])

##
fig2 = plot((1:paramsTrue["sample size"])/paramsTrue["sample size"], nVarSimS_f,
    seriestype=:sticks,
    yscale=:log10,
)
plot!(dfsTrueS.freqs_f[2:end], dfsTrueS.n_f[2:end])

xlims!(0,0.15)
# ylims!(5E-1, maximum(nVarSimS_f))
ylims!(1E3, 1E6)

display(fig2)





# NMature = paramsTrue["N final"]
# μTrue = paramsTrue["μ"]
# # nPlots = 50

# fig1 = histogram(mSimS_cid,
#     bins=30,
#     label=L"\tilde{\mu}="*string(round(paramsIn["μ"], digits=2)),
#     legendfontsize=11,
#     # linewidth=2,
#     # legend=:topleft
#     )
# xlabel!("# mutations")
# ylabel!("# cells")
# title!("μ true = $μTrue, N mature = $NMature")

# # mMax = Int(round(2*paramsTrue["λ"]*paramsTrue["evolve time"]*paramsTrue["μ"]*3))
# # xlims!(mMax/8,3*mMax/8)
# # xlims!(300,550)
# display(fig1)

# ##
# NMature = paramsTrue["N final"]
# μ = paramsTrue["μ"]
# # gr()
# # fig = heatmap(pCont_p, NCont_N, transpose(vaf1ErrorInterpol_p_N),
# #     fillalpha=0.6)
# # # plot!(pCont_p, NCont_N, transpose(vaf1ErrorInterpol_p_N))
# # # plot!(p_p, Nopt_p,
# # #     linewidth=2,
# # #     linestyle=:dash,
# # #     label="best fit",
# # #     color=1)
# # xlims!(0,1)
# # # ylims!(10^4,2*10^5)
# # xlabel!("p (fraction of asymmetric divisions)")
# # ylabel!("N")
# # # title!("μ = 1.2")
# # title!("error (number of variants at 1/S)")
# # display(fig)

# fig2 = heatmap(pCont_p, NCont_N, transpose(abs.(vaf1ErrorInterpol_p_N)),
# fillalpha=0.5)
# plot!(pCont_p, NCont_N, transpose(abs.(vaf1ErrorInterpol_p_N)))
# plot!(p_p, Nopt_p,
# linewidth=2,
# linestyle=:dash,
# label="best fit",
# colorbar_title="relative error (number of variants at 1/S)",
# color=1)
# scatter!([paramsTrue["p"],], [paramsTrue["N final"],],
#     # markerstyle=
#     label="True values")
# xlims!(0,1)
# # ylims!(10^4,2*10^5)
# xlabel!("p (fraction of asymmetric divisions)")
# ylabel!("N")
# # title!("relative error (number of variants at 1/S)")
# # title!("μ = 1.2")
# # savefig("figures/ABCfitting/N-p_spaceAbs_mu"*μString*".pdf")
# display(fig2)

# SAVEPLOTS ? savefig(fig2, "figures/ABCfitting/N-p_spaceAbs_mu$μ _N$NMature"*".pdf") : nothing

# ##
# # pyplot()
# fig3 = bar((1:paramsTrue["sample size"])/paramsTrue["sample size"], nVarSimS_f[2:end-1],
#     # linewidth=0.8,
#     linewidth=0,
#     color=:black,
#     fillalpha=0.4,
#     # fillrange=0
#     label="data"
#     )
# pTrue = paramsTrue["p"]
# pTest_p = 0.1:0.2:0.9
# for pTest in pTest_p
#     paramsEvo = createTestParams(paramsIn, NoptSpl_p(pTest), pTest)
#     vfs = VAFDyn.VFreqspace(Int(round(NoptSpl_p(pTest))), 300)
#     VAFDyn.evolveGrowingVAF(vfs, paramsEvo, paramsEvo["evolve time"])
#     dfsPDE = VAFDyn.makeDFSfromVFS(vfs, paramsEvo["N final"])
#     dfsPDES = VAFDyn.sampler(dfsPDE, paramsEvo["sample size"])
#     plot!(dfsPDES.freqs_f[2:end-1], dfsPDES.n_f[2:end-1],
#     yscale=:log10,
#     label="PDE: p = $pTest")
# end
# ylims!(1E-1, 3E4)
# # title!("true p = $pTrue")
# display(fig3)
# SAVEPLOTS ? savefig(fig3, "figures/ABCfitting/VAF1Fit_mu$μ _N$NMature"*".pdf") : nothing

# ##
# fig4 = plot(fig1, fig2, fig3, layout=3, size=(1000,800))
# SAVEPLOTS ? savefig(fig4, "figures/SinglePatientPipeline/singlePatientResults_mu$μ _N$NMature"*".pdf") : nothing
# display(fig4)

##

# pTestA = 0.4
# paramsEvoA = createTestParams(paramsIn, paramsTrue["N final"], pTestA)
# vfsA = VAFDyn.VFreqspace(paramsEvoA["N final"], 500)
# VAFDyn.evolveGrowingVAF(vfsA, paramsEvoA, paramsEvoA["evolve time"])
# dfsPDEA = VAFDyn.makeDFSfromVFS(vfsA, paramsEvoA["N final"])
# dfsPDEAS = VAFDyn.sampler(dfsPDEA, paramsEvoA["sample size"])

# println(dfsPDEAS.n_f[2])

