module InferencePipeline


using Statistics, Distributions, ProgressMeter, Dierckx, Roots, Interpolations
export calcNpSpace, calcNpSpaceMultPopDistribution

include("./vafdyn.jl")
using .VAFDyn
include("./theory.jl")
using .Theory


function estimateRates(mVars_cid)
    varBurdenS = var(mVars_cid)
    meanBurdenS = mean(mVars_cid)
    # standardErrorMean = std(mVars_cid)/sqrt(length(mVars_cid))
    stdMuts = std(mVars_cid)
    r, μ = Theory.getCPRateParamsFromBurdenStats(meanBurdenS, varBurdenS)
    paramsEst = Dict{String, Real}(
        "divisions" => r,
        "μ" => μ,
        # "ste mutations" => standardErrorMean,
        "σ mutations" => stdMuts
    )
    return paramsEst
end

function estimateRates(mVars_cid::Vector{<:Integer}, μKnown::Real)
    varBurdenS = var(mVars_cid)
    meanBurdenS = mean(mVars_cid)
    # standardErrorMean = std(mVars_cid)/sqrt(length(mVars_cid))
    stdMuts = std(mVars_cid)
    r, μ = Theory.getCPRateParamsFromBurdenStats(meanBurdenS, varBurdenS, μKnown=μKnown)
    paramsEst = Dict{String, Real}(
        "divisions" => r,
        "μ" => μ,
        # "ste mutations" => standardErrorMean
        "σ mutations" => stdMuts,
    )
    return paramsEst
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
    seλ = Theory.getSeλFromTotalMutations(paramsTest)
    paramsTest["seλ"] = seλ
    Theory.extendParams!(paramsTest)
    return paramsTest
end

# ---------- create VAF spectra ----------
function calcExpectedVAFs(paramsIn::Dict, _N::AbstractArray, _p::AbstractArray, lengthVFS::Int, vafFit::Int=1; cumulative::Bool=false)
    nVafFit_p_N = Array{Float64,2}(undef, length(_p), length(_N))
    for (i,p) in enumerate(_p)
        for (j,N) in enumerate(_N)
            paramsParticle = createTestParams(paramsIn, N, p)
            vfs = VAFDyn.VFreqspace(paramsParticle["N final"], lengthVFS)
            VAFDyn.evolveGrowingVAF(vfs, paramsParticle, paramsParticle["evolve time"])
            dfs = VAFDyn.makeDFSfromVFS(vfs, paramsParticle["N final"])
            dfsS = VAFDyn.sampler(dfs, paramsParticle["sample size"])
            if cumulative
                nVafFit_p_N[i,j] = sum(dfsS.n_f[(1+vafFit):end])
            else
                nVafFit_p_N[i,j] = dfsS.n_f[1+vafFit]
            end
        end
    end
    return nVafFit_p_N
end

function calcExpectedVAFs(paramsIn::Dict, _N::AbstractArray, p::Real, lengthVFS::Int, vafFit::Int=1; cumulative::Bool=false)
    nVafFit_N = Vector{Float64}(undef, length(_N))
    for (j,N) in enumerate(_N)
        paramsParticle = createTestParams(paramsIn, N, p)
        vfs = VAFDyn.VFreqspace(paramsParticle["N final"], lengthVFS)
        VAFDyn.evolveGrowingVAF(vfs, paramsParticle, paramsParticle["evolve time"])
        dfs = VAFDyn.makeDFSfromVFS(vfs, paramsParticle["N final"])
        dfsS = VAFDyn.sampler(dfs, paramsParticle["sample size"])
        if cumulative
            nVafFit_N[j] = sum(dfsS.n_f[(1+vafFit):end])
        else
            nVafFit_N[j] = dfsS.n_f[1+vafFit]
        end
    end
    return nVafFit_N
end

function findZeroErrorPoint(nVafTheoryfit_N, nVafDataFit::Int, _N)
    vaf1Error_N = (nVafTheoryfit_N .- nVafDataFit)/nVafDataFit
    vaf1ErrorInterpol_N = Spline1D(_N, vaf1Error_N)
    Nopt = try
        find_zero(n->vaf1ErrorInterpol_N(n), (_N[1],_N[end]), Bisection())
    catch e
        println("error: no zero found for p="*string(_p[i])*" between "*string(_N[1])*" and "*string(_N[end])*".")
        if vaf1ErrorInterpol_N(_N[end]) < 0
            _N[end]
        else
            _N[1]
        end
    end
    return Nopt
end

# ---------- get interpolated solutions ----------
function findZeroErrorLine(nVarTheoryfit_p_N, nVarDataFit::Int, _N, _p; l::Int=100, verbose=false)
    
    nVarError_p_N = (nVarTheoryfit_p_N .- nVarDataFit)/nVarDataFit
    nVarErrorSpl_p_N = Spline2D(_p, _N, nVarError_p_N)
    
    verbose ? display(nVarError_p_N) : nothing

    Nopt_p = Vector{Float64}(undef, length(_p))
    for i in 1:length(_p)
        vaf1ErrorInterpol_N = Spline1D(_N, nVarError_p_N[i,:])
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
    nVarErrorInterpol_p_N = Array{Float64}(undef, l, l)
    for (i,p) in enumerate(pCont_p)
        for (j,N) in enumerate(NCont_N)
            nVarErrorInterpol_p_N[i,j] = nVarErrorSpl_p_N(p, N)
        end
    end
    # return Nopt_p, NoptSpl_p, pOpt_N, pOptSpl_N, pCont_p, NCont_N, NoptInterpol_p, vaf1ErrorInterpol_p_N
    return Nopt_p, NoptSpl_p, pCont_p, NCont_N, NoptInterpol_p, nVarErrorInterpol_p_N
end

# =================== Full Pipeline ===================

function calcNpSpace(paramsKnown::Dict, nVarS_f::Vector{Int}, mVarsS_cid::Vector{Int}, _N, _p, lVFS::Int, vafFit::Int=1; lCont::Int=100, verbose=false, cumulative::Bool=false)
    if haskey(paramsKnown, "μ")
        paramsEst = estimateRates(mVarsS_cid, paramsKnown["μ"])
    else
        paramsEst = estimateRates(mVarsS_cid)
    end
    calcNpSpace(paramsKnown, paramsEst, nVarS_f, _N, _p, lVFS, vafFit; lCont, verbose, cumulative)
end

function calcNpSpace(paramsKnown::Dict, paramsEst::Dict, nVarS_f, _N, _p, lVFS::Int, vafFit::Int=1; lCont::Int=100, verbose=false, cumulative::Bool=false)
    paramsIn = merge(paramsKnown, paramsEst)
    verbose ? display(paramsIn) : nothing
    nVarTheoryFit_p_N = calcExpectedVAFs(paramsIn, _N, _p, lVFS, vafFit; cumulative)
    nVarRef = cumulative ? sum(nVarS_f[1+vafFit:end]) : nVarS_f[1+vafFit]
    _, _, pCont_p, NCont_N, NOptInterpol_p, vaf1ErrorInterpol_p_N = 
        findZeroErrorLine(nVarTheoryFit_p_N, nVarRef, _N, _p; l=lCont, verbose=verbose)
    return NCont_N, pCont_p, NOptInterpol_p, vaf1ErrorInterpol_p_N
end

function calcNpSpaceMultPopDistribution(_N, _p, paramsKnown, nVS_Pop_f, nVS_Pop_cid, lVFS, vafFit::Int=1; verbose=false)
    nPops = length(nVS_Pop_f)
    lCont = 100
    NOpt_sim_p = Array{Float64}(undef, nPops, lCont)
    NCont_N, pCont_p, NOptInterpol_p, ___ = calcNpSpace(paramsKnown, nVS_Pop_f[1], nVS_Pop_cid[1], _N, _p, lVFS, vafFit; lCont=lCont, verbose=verbose)
    NOpt_sim_p[1, :] = NOptInterpol_p
    @showprogress for sim in 2:nPops
        NCont_N, pCont_p, NOptInterpol_p, ___ = calcNpSpace(paramsKnown, nVS_Pop_f[sim], nVS_Pop_cid[sim], _N, _p, lVFS, vafFit; lCont=lCont, verbose=verbose)
        NOpt_sim_p[sim, :] = NOptInterpol_p
    end
    # NOptAv_p = vec(mean(NOpt_sim_p, dims=1))
    # NOptStd_p = vec(std(NOpt_sim_p, dims=1))
    # NOpt95L_p = NOptAv_p .- 2*NOptStd_p
    # NOpt95U_p = NOptAv_p .+ 2*NOptStd_p

    return NCont_N, pCont_p, NOpt_sim_p
end

# ================== Max pop size t_M - N_p =======================

function calcTmNhSpace(paramsKnown::Dict, nVarS_f::Vector{Int}, mVarsS_cid::Vector{Int}, _N, p, _tM, _Nh, lVFS::Int, vafFit::Int=1; lCont::Int=100, verbose=false)
    if haskey(paramsKnown, "μ")
        paramsEst = estimateRates(mVarsS_cid, paramsKnown["μ"])
    else
        paramsEst = estimateRates(mVarsS_cid)
    end
    calcTmNhSpace(paramsKnown, paramsEst, nVarS_f, _N, p, _tM, _Nh, lVFS, vafFit; lCont, verbose)
end

function calcTmNhSpace(paramsKnown::Dict, paramsEst::Dict, nVarS_f, _N, p, _tM, _Nh, lVFS::Int, vafFit::Int=1; lCont::Int=100, verbose=false)
    paramsIn = merge(paramsKnown, paramsEst)
    verbose ? display(paramsIn) : nothing

    NmaxOpt_tM_nh = Array{Float64, 2}(undef, length(_tM), length(_Nh))
    for (i,tM) in enumerate(_tM)
        paramsIn["mature time"] = tM
        for (j,Nh) in enumerate(_Nh)
            verbose && println("testing tM="*string(tM)*"; Nh="*string(Nh))
            paramsIn["pure births"] = Nh
            nVarTheoryFit_N = calcExpectedVAFs(paramsIn, _N, p, lVFS, vafFit)
            nVarRef = nVarS_f[1+vafFit]
            NmaxOpt_tM_nh[i,j] = findZeroErrorPoint(nVarTheoryFit_N, nVarRef, _N)
            verbose && println("N max for tM="*string(tM)*" and Nh="*string(Nh)*" is: ", NmaxOpt_tM_nh[i,j])
        end
    end
    NmaxSpl_tM_nh = Spline2D(_tM, _Nh, NmaxOpt_tM_nh)
    _tMCont = range(_tM[1], _tM[end], length=lCont)
    _NhCont = range(_Nh[1], _Nh[end], length=lCont)
    NmaxInterpol_tM_nh = Array{Float64}(undef, lCont, lCont)
    for (i,tM) in enumerate(_tMCont)
        for (j,Nh) in enumerate(_NhCont)
            NmaxInterpol_tM_nh[i,j] = NmaxSpl_tM_nh(tM, Nh)
        end
    end
    return _tMCont, _NhCont, NmaxInterpol_tM_nh, NmaxOpt_tM_nh
end




end