##
using Plots
pyplot()
# gr()
# plotlyjs()
# plotly()
##
using JLD2, FileIO
using Statistics
using Distributions
using ProgressMeter
using Dierckx, Roots, Interpolations
using LaTeXStrings
using Glob
using StatsBase
##
include("../../src/vafdyn.jl")
using .VAFDyn
include("../../src/theory.jl")
using .Theory
include("../../src/inferencePipeline.jl")
using .InferencePipeline
##
SAVEPLOTS=false
LOADDATA = true

# ==================== Load data ====================

if LOADDATA
    @load "data/Simulations/npSpaceDistribution_N10000.jld2"
else

    # ----- Load sim data -----
    fileNames_ = glob("singlePatientFullSim_Ni1_*.jld2", "./data/Simulations/Nf10000")
    # fileNames_ = fileNames_[1:4]
    nFiles = length(fileNames_)
    @load fileNames_[1] paramsTrue timesSim_ nCellSim_t

    # nV_Sim_f = Vector{Vector{Int64}}(undef, nFiles)
    nVS_Sim_f = Vector{Vector{Int64}}(undef, nFiles)
    nVS_Sim_cid = Vector{Vector{Int64}}(undef, nFiles)
    for (i,fName) in enumerate(fileNames_)
        @load fName nVarSim_f nVarSimS_f mSimS_cid
        # nV_Sim_f[i] = nVarSim_f
        nVS_Sim_f[i] = nVarSimS_f
        nVS_Sim_cid[i] = mSimS_cid
    end

    ## ----- Build N-p space -----
    N_N = range(500, 40000, length=8)
    p_p = range(0.1, 0.9, length=8)
    lVfs = 400

    paramsKnown = Dict{String,Real}(
        "N initial" => 1,
        "sample size" => paramsTrue["sample size"],
        "mature time" => paramsTrue["mature time"],
        "evolve time" => paramsTrue["evolve time"],
        "μ" => paramsTrue["μ"]
    )

    NCont_N, pCont_p, NOpt_sim_p = InferencePipeline.calcNpSpaceMultPopDistribution(N_N, p_p, paramsTrue, nVS_Sim_f, nVS_Sim_cid, lVfs, verbose=true)

    filename = "data/Simulations/npSpaceDistribution_N"*string(paramsTrue["N final"])*".jld2"
    save(filename, 
        Dict(
            "paramsTrue" => paramsTrue,
            "_p" => pCont_p,
            "_N" => NCont_N,
            "NOpt_sim_p" => NOpt_sim_p,
        )
    )
end

NOptAv_p = vec(mean(NOpt_sim_p, dims=1))
NOptStd_p = vec(std(NOpt_sim_p, dims=1))
NOpt95L_p = NOptAv_p .- 2*NOptStd_p
NOpt95U_p = NOptAv_p .+ 2*NOptStd_p


## ================================ Plotting ===================================
pyplot()
theme(:default, minorgrid=false, gridstyle=:dash, fontfamily="DejaVu Sans")

##

fig1 = plot(_p, NOptAv_p, color=:black)
plot!(_p, NOpt95L_p,
linestyle=:dash,
color=:black,
label="95% confidence interval")
plot!(_p, NOpt95U_p,
color=:black,
linestyle=:dash,
label=""
)
xlabel!("p")
ylabel!("N")
display(fig1)

##

fig2 = histogram(NOpt_sim_p[:, 40], bins=15, normalize=true,
color=:grey65,
label="")
xlabel!("N")
ylabel!("density of results")
title!("p="*string(paramsTrue["p"]))

display(fig2)


##
freqsS_f = (1:paramsTrue["sample size"])/paramsTrue["sample size"]


fig3 = plot(yscale=:log10)
for sim in 1:5
    bar!(freqsS_f, nVS_Sim_f[sim][2:end],
    color=2+sim,
    # linewidth=0.8,
    linewidth=0,
    fillalpha=0.5,
    # fillrange=0
    label="sim "*string(sim)
    )
end
pTrue = paramsTrue["p"]
# pTest_p = 0.1:0.2:0.9
# pTest_p = [0.5,]
# for pTest in pTest_p
#     paramsEvo = createTestParams(paramsIn, NoptSpl_p(pTest), pTest)
#     vfs = VAFDyn.VFreqspace(Int(round(NoptSpl_p(pTest))), 300)
#     VAFDyn.evolveGrowingVAF(vfs, paramsEvo, paramsEvo["evolve time"])
#     dfsPDE = VAFDyn.makeDFSfromVFS(vfs, paramsEvo["N final"])
#     dfsPDES = VAFDyn.sampler(dfsPDE, paramsEvo["sample size"])
#     plot!(dfsPDES.freqs_f[2:end-1], dfsPDES.n_f[2:end-1],
#     color=:black,
#     yscale=:log10,
#     label="PDE: p = $pTest")
# end
ylims!(1E-1, 3E4)
# title!("true p = $pTrue")
display(fig3)
SAVEPLOTS && savefig(fig3, "figures/ABCfitting/VAF1Fit_mu$μ _N$NMature"*".pdf")

##
fig4 = plot(fig1, fig2, fig3, layout=3, size=(1000,800))
SAVEPLOTS && savefig(fig4, "figures/SinglePatientPipeline/singlePatientResults_mu$μ _N$NMature"*".pdf")
display(fig4)

##

# pTestA = 0.4
# paramsEvoA = createTestParams(paramsIn, paramsTrue["N final"], pTestA)
# vfsA = VAFDyn.VFreqspace(paramsEvoA["N final"], 500)
# VAFDyn.evolveGrowingVAF(vfsA, paramsEvoA, paramsEvoA["evolve time"])
# dfsPDEA = VAFDyn.makeDFSfromVFS(vfsA, paramsEvoA["N final"])
# dfsPDEAS = VAFDyn.sampler(dfsPDEA, paramsEvoA["sample size"])

# println(dfsPDEAS.n_f[2])

## ===== heatmap =====

# histogram2d(NOpt_sim_p)
using StatsBase

_NDist = 1e3:1e3:2e4
NOptDistribution_N_p = Array{Float64,2}(undef, length(_NDist)-1, length(_p))
for (i,p) in enumerate(_p)
    h = fit(Histogram, NOpt_sim_p[:,i], _NDist, closed=:right)
    h = StatsBase.normalize(h, mode=:density)
    NOptDistribution_N_p[:, i] = h.weights
end


##
using Colors, ColorSchemes

    # mygrays = ColorScheme(append!([RGB{Float64}(1, 1, 1),],[RGB{Float64}(i, i, i) for i in 1:-0.1:0]))

colorgrad = cgrad(:grayC, rev = false, alpha = nothing, scale = nothing, categorical = nothing)

contour(_p, _NDist[2:end], NOptDistribution_N_p, fill = true, c=colorgrad)