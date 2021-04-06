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
##
include("../../src/vafdyn.jl")
using .VAFDyn
include("../../src/theory.jl")
using .Theory

##
SAVEPLOTS=false


## ============= create analysis functions =============

# --------- Estimate division and mutation rate ---------
function estimateRates(paramsTrue, mVars_cid)
    paramsKnown = Dict{String, Real}(
        "evolve time" => paramsTrue["evolve time"],
        "N initial" => paramsTrue["N initial"],
        "sample size" => paramsTrue["sample size"],
        "mature time" => paramsTrue["mature time"]
    )
    varBurdenS = var(mVars_cid)
    meanBurdenS = mean(mVars_cid)
    r, μ = Theory.getCPRateParamsFromBurdenStats(meanBurdenS, varBurdenS)
    paramsEst = Dict{String, Real}(
        "divisions" => r,
        "μ" => μ
    )

    paramsIn = merge(paramsKnown, paramsEst)
    return paramsIn
end

function createTestParams(params::Dict, N::Real, p::Real)
    r = params["divisions"]
    t = params["evolve time"]
    tM = params["mature time"]
    Nf = Int(round(N))
    γ = log(Nf)/tM

    paramsTest = deepcopy(params)
    paramsTest["N final"] = Nf
    paramsTest["p"] = p
    paramsTest["growth rate"] = γ
    λ = Theory.getλFromTotalDivisions(paramsTest)
    paramsTest["λ"] = λ
    Theory.extendParams!(paramsTest)
    return paramsTest
end

# ---------- create VAF spectra ----------
function calcExpectedVAFs(paramsIn::Dict, _N::AbstractArray, _p::AbstractArray, lengthVFS::Int)
    nVAF1_p_N = Array{Float64,2}(undef, length(_p), length(_N))
    for (i,p) in enumerate(_p)
        for (j,N) in enumerate(_N)
            paramsParticle = createTestParams(paramsIn, N, p)
            vfs = VAFDyn.VFreqspace(paramsParticle["N final"], lengthVFS)
            VAFDyn.evolveGrowingVAF(vfs, paramsParticle, paramsParticle["evolve time"])
            dfs = VAFDyn.makeDFSfromVFS(vfs, paramsParticle["N final"])
            dfsS = VAFDyn.sampler(dfs, paramsParticle["sample size"])
            nVAF1_p_N[i,j] = dfsS.n_f[2]
        end
    end
    return nVAF1_p_N
end

# ---------- get interpolated solutions ----------
function findZeroErrorLine(nVAF1_p_N, nData_f, _N, _p; l::Int=100, verbose=false)
    
    vaf1Error_p_N = (nVAF1_p_N .- nData_f[2])/nData_f[2]
    vaf1ErrorSpl_p_N = Spline2D(_p, _N, vaf1Error_p_N)
    
    verbose ? display(vaf1Error_p_N) : nothing

    Nopt_p = Vector{Float64}(undef, length(_p))
    for i in 1:length(_p)
        vaf1ErrorInterpol_N = Spline1D(_N, vaf1Error_p_N[i,:])
        Nopt_p[i] = try
            find_zero(n->vaf1ErrorInterpol_N(n), (_N[1],_N[end]), Bisection())
        catch e
            println("error: no zero found for p="*string(_p[i])*" between "*string(_N[1])*" and "*string(_N[end])*".")
            if vaf1ErrorInterpol_N(_N[end]) < 0
                _N[end]
            else
                _N[1]
            end
        end
        verbose ? println("Nopt_p = ", Nopt_p[i]) : nothing
    end
    NoptSpl_p = Spline1D(_p, Nopt_p)
    
    pCont_p = range(_p[1], _p[end], length=l)
    NCont_N = range(_N[1], _N[end], length=l)

    NoptInterpol_p = NoptSpl_p.(pCont_p)
    vaf1ErrorInterpol_p_N = Array{Float64}(undef, l, l)
    for (i,p) in enumerate(pCont_p)
        for (j,N) in enumerate(NCont_N)
            vaf1ErrorInterpol_p_N[i,j] = vaf1ErrorSpl_p_N(p, N)
        end
    end
    # return Nopt_p, NoptSpl_p, pOpt_N, pOptSpl_N, pCont_p, NCont_N, NoptInterpol_p, vaf1ErrorInterpol_p_N
    return Nopt_p, NoptSpl_p, pCont_p, NCont_N, NoptInterpol_p, vaf1ErrorInterpol_p_N
end








## =================== Full Pipeline ===================

function calcSingleSimNpSpace(paramsTrue::Dict, nVarS_f, mVarsS_cid::Vector{Int}, _N, _p, lVFS::Int; lCont::Int=100, verbose=false)
    paramsIn = estimateRates(paramsTrue, mVarsS_cid)
    verbose ? display(paramsIn) : nothing
    nVAF1_p_N = calcExpectedVAFs(paramsIn, _N, _p, lVFS)
    Nopt_p, NoptSpl_p, pCont_p, NCont_N, NOptInterpol_p, vaf1ErrorInterpol_p_N = 
        findZeroErrorLine(nVAF1_p_N, nVarS_f, _N, _p; l=lCont, verbose=verbose)
    return NCont_N, pCont_p, NOptInterpol_p, vaf1ErrorInterpol_p_N
end
                                                                                                   
##
function calcNpSpaceDistribution(_N, _p, paramsTrue, nVS_Sim_f, nVS_Sim_cid, lVFS; verbose=false)
    nSims = length(nVS_Sim_f)
    lCont = 100
    NOpt_sim_p = Array{Float64}(undef, nSims, lCont)
    NCont_N, pCont_p, NOptInterpol_p, ___ = calcSingleSimNpSpace(paramsTrue, nVS_Sim_f[1], nVS_Sim_cid[1], _N, _p, lVFS; lCont, verbose=verbose)
    NOpt_sim_p[1, :] = NOptInterpol_p
    @showprogress for sim in 2:nSims
        NCont_N, pCont_p, NOptInterpol_p, ___ = calcSingleSimNpSpace(paramsTrue, nVS_Sim_f[sim], nVS_Sim_cid[sim], _N, _p, lVFS; lCont, verbose=verbose)
        NOpt_sim_p[sim, :] = NOptInterpol_p
    end
    # NOptAv_p = vec(mean(NOpt_sim_p, dims=1))
    # NOptStd_p = vec(std(NOpt_sim_p, dims=1))
    # NOpt95L_p = NOptAv_p .- 2*NOptStd_p
    # NOpt95U_p = NOptAv_p .+ 2*NOptStd_p

    return NCont_N, pCont_p, NOpt_sim_p
end



## ============================================================================================== ==================================== implement pipeline =================================== ==============================================================================================




# ========== Load data ==========

fileNames_ = glob("singlePatient*.jld2", "./data/SinglePatientPipeline/Nf10000")
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


## =============== Build N-p space ===============
# N_N = 1E2:5E2:9E3
N_N = 499:2500:26000
p_p = 0.05:0.2:0.85

##
NCont_N, pCont_p, NOpt_sim_p = calcNpSpaceDistribution(N_N, p_p, paramsTrue, nVS_Sim_f, nVS_Sim_cid, 500, verbose=true)

##
filename = "data/SinglePatientPipeline/npSpaceDistribution_N"*string(paramsTrue["N final"])*".jld2"
save(filename, Dict(
    "_p" => _p,
    "_N" => _N,
    "NOpt_sim_p" => NOpt_sim_p
)
)


##

NOptAv_p = vec(mean(NOpt_sim_p, dims=1))
NOptStd_p = vec(std(NOpt_sim_p, dims=1))
NOpt95L_p = NOptAv_p .- 2*NOptStd_p
NOpt95U_p = NOptAv_p .+ 2*NOptStd_p


## ================================ Plotting ===================================



fig1 = plot(pCont_p, NOptAv_p, color=:black)
plot!(pCont_p, NOpt95L_p,
linestyle=:dash,
color=:black,
label="95% confidence interval")
plot!(pCont_p, NOpt95U_p,
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

