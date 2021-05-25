##
using JLD2, FileIO
using Statistics
using LaTeXStrings
using Glob
include("../../src/vafdyn.jl")
using .VAFDyn
include("../../src/theory.jl")
using .Theory

##
SAVEPLOTS=false

# ========== Load data ==========

fileNames_ = glob("singlePatient*.jld2", "./data/SinglePatientPipeline/Nf10000")
nSims = length(fileNames_)
@load fileNames_[1] paramsTrue timesSim_ nCellSim_t
nV_sim_f = Array{Float64, 2}(undef, nSims, paramsTrue["N final"]+1)
nVS_sim_f = Array{Float64, 2}(undef, nSims, paramsTrue["sample size"]+1)
for (i,fName) in enumerate(fileNames_)
    @load fName nVarSim_f nVarSimS_f
    nV_sim_f[i, :] .= nVarSim_f
    nVS_sim_f[i, :] .= nVarSimS_f
end

filenameRes = "data/SinglePatientPipeline/Nf"*string(paramsTrue["N final"])*"/npSpaceDistribution_N"*string(paramsTrue["N final"])*".jld2"
@load filenameRes _p _N NOpt_sim_p

NOptAv_p = vec(mean(NOpt_sim_p, dims=1))
NOptStd_p = vec(std(NOpt_sim_p, dims=1))
NOpt95L_p = NOptAv_p .- 2*NOptStd_p
NOpt95U_p = NOptAv_p .+ 2*NOptStd_p

## ================================ Plotting ===================================

using Plots
pyplot(legendfontsize=10, guidefontsize=12, tickfontsize=10,  size=(600,400))
using ColorSchemes
#3
myGrad = cgrad(:lapaz, rev = true)

fig1 = histogram2d(pOpts_id, NOpts_id,
normalize=true,
bins=25,
# c=:lapaz
# c=myGrad
c=:bilbao
)
scatter!([paramsTrue["p"]], [paramsTrue["N final"]],
color=:black,
markersize=6,
label="true value")
xlabel!("p")
ylabel!("N")
display(fig1)

figname = "figures/paper/simsNpInference"
savefig(fig1, figname*"_alt.pdf")

##
# fig1 = plot(pCont_p, NOptAv_p, color=:black)
# plot!(pCont_p, NOpt95L_p,
# linestyle=:dash,
# color=:black,
# label="95% confidence interval")
# plot!(pCont_p, NOpt95U_p,
# color=:black,
# linestyle=:dash,
# label=""
# )
# xlabel!("p")
# ylabel!("N")
# display(fig1)

# ##

# fig2 = histogram(NOpt_sim_p[:, 40], bins=15, normalize=true,
# color=:grey65,
# label="")
# xlabel!("N")
# ylabel!("density of results")
# title!("p="*string(paramsTrue["p"]))

# display(fig2)


# ##
# freqsS_f = (1:paramsTrue["sample size"])/paramsTrue["sample size"]


# fig3 = plot(yscale=:log10)
# for sim in 1:5
#     bar!(freqsS_f, nVS_Sim_f[sim][2:end],
#     color=2+sim,
#     # linewidth=0.8,
#     linewidth=0,
#     fillalpha=0.5,
#     # fillrange=0
#     label="sim "*string(sim)
#     )
# end
# pTrue = paramsTrue["p"]
# # pTest_p = 0.1:0.2:0.9
# # pTest_p = [0.5,]
# # for pTest in pTest_p
# #     paramsEvo = createTestParams(paramsIn, NoptSpl_p(pTest), pTest)
# #     vfs = VAFDyn.VFreqspace(Int(round(NoptSpl_p(pTest))), 300)
# #     VAFDyn.evolveGrowingVAF(vfs, paramsEvo, paramsEvo["evolve time"])
# #     dfsPDE = VAFDyn.makeDFSfromVFS(vfs, paramsEvo["N final"])
# #     dfsPDES = VAFDyn.sampler(dfsPDE, paramsEvo["sample size"])
# #     plot!(dfsPDES.freqs_f[2:end-1], dfsPDES.n_f[2:end-1],
# #     color=:black,
# #     yscale=:log10,
# #     label="PDE: p = $pTest")
# # end
# ylims!(1E-1, 3E4)
# # title!("true p = $pTrue")
# display(fig3)
# SAVEPLOTS ? savefig(fig3, "figures/ABCfitting/VAF1Fit_mu$μ _N$NMature"*".pdf") : nothing

# ##
# fig4 = plot(fig1, fig2, fig3, layout=3, size=(1000,800))
# SAVEPLOTS ? savefig(fig4, "figures/SinglePatientPipeline/singlePatientResults_mu$μ _N$NMature"*".pdf") : nothing
# display(fig4)

# ##

# # pTestA = 0.4
# # paramsEvoA = createTestParams(paramsIn, paramsTrue["N final"], pTestA)
# # vfsA = VAFDyn.VFreqspace(paramsEvoA["N final"], 500)
# # VAFDyn.evolveGrowingVAF(vfsA, paramsEvoA, paramsEvoA["evolve time"])
# # dfsPDEA = VAFDyn.makeDFSfromVFS(vfsA, paramsEvoA["N final"])
# # dfsPDEAS = VAFDyn.sampler(dfsPDEA, paramsEvoA["sample size"])

# # println(dfsPDEAS.n_f[2])

