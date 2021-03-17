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
    r, μ = Theory.getCPRateParamsFromBurdenStats(meanBurdenS, varBurdenS)
    paramsEst = Dict{String, Real}(
        "divisions" => r,
        "μ" => μ
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

# =================== Full Pipeline ===================

function calcNpSpace(paramsKnown::Dict, nVarS_f, mVarsS_cid::Vector{Int}, _N, _p, lVFS::Int; lCont::Int=100, verbose=false)
    paramsEst = estimateRates(mVarsS_cid)
    paramsIn = merge(paramsKnown, paramsEst)
    verbose ? display(paramsIn) : nothing
    nVAF1_p_N = calcExpectedVAFs(paramsIn, _N, _p, lVFS)
    Nopt_p, NoptSpl_p, pCont_p, NCont_N, NOptInterpol_p, vaf1ErrorInterpol_p_N = 
        findZeroErrorLine(nVAF1_p_N, nVarS_f, _N, _p; l=lCont, verbose=verbose)
    return NCont_N, pCont_p, NOptInterpol_p, vaf1ErrorInterpol_p_N
end


function calcNpSpaceMultPopDistribution(_N, _p, paramsTrue, nVS_Pop_f, nVS_Pop_cid, lVFS; verbose=false)
    nPops = length(nVS_Pop_f)
    lCont = 100
    NOpt_sim_p = Array{Float64}(undef, nPops, lCont)
    NCont_N, pCont_p, NOptInterpol_p, ___ = calcNpSpace(paramsTrue, nVS_Pop_f[1], nVS_Pop_cid[1], _N, _p, lVFS; lCont, verbose=verbose)
    NOpt_sim_p[1, :] = NOptInterpol_p
    @showprogress for sim in 2:nPops
        NCont_N, pCont_p, NOptInterpol_p, ___ = calcNpSpace(paramsTrue, nVS_Pop_f[sim], nVS_Pop_cid[sim], _N, _p, lVFS; lCont, verbose=verbose)
        NOpt_sim_p[sim, :] = NOptInterpol_p
    end
    # NOptAv_p = vec(mean(NOpt_sim_p, dims=1))
    # NOptStd_p = vec(std(NOpt_sim_p, dims=1))
    # NOpt95L_p = NOptAv_p .- 2*NOptStd_p
    # NOpt95U_p = NOptAv_p .+ 2*NOptStd_p

    return NCont_N, pCont_p, NOpt_sim_p
end





end