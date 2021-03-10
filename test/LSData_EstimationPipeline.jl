include("../src/vafdyn.jl")
using .VAFDyn
include("../src/theory.jl")
using .Theory
include("../src/compoundPoisson.jl")
using .CompoundPoisson
include("../src/burdenDyn.jl")
using .BurdenDyn

using DelimitedFiles, FileIO
using ProgressMeter
using ApproxBayes
using Distributions, Distances
using Dierckx, Roots, Interpolations

using Plots
gr()
# using Gadfly

## =========================== pipeline ==============================

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
    extendParams!(paramsTest)
    return paramsTest
end

# ---------- create VAF spectra ----------
function calcExpectedVAFs(paramsIn::Dict, _N::AbstractArray, _p::AbstractArray, lengthVFS::Int)
    nVAF1_p_N = Array{Float64,2}(undef, length(_p), length(_N))
    @showprogress for (i,p) in enumerate(_p)
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

function calcSingleSimNpSpace(paramsIn::Dict, nVarS_f, _ND, _pD, lVFS::Int; lCont::Int=100, verbose=false)
    verbose ? display(paramsIn) : nothing
    nVAF1_p_N = calcExpectedVAFs(paramsIn, _ND, _pD, lVFS)
    Nopt_p, NoptSpl_p, _p, _N, NOptInterpol_p, vaf1ErrorInterpol_p_N = 
        findZeroErrorLine(nVAF1_p_N, nVarS_f, _ND, _pD; l=lCont, verbose=verbose)
    return _N, _p, NOptInterpol_p, vaf1ErrorInterpol_p_N
end


## ============================== Load Lee-Six Data ==============================

let
    nVData_cid, SCBurdenMean, SCBurdenVar, freqsData_f, nVData_f, sampleSize = load("data/LSDataStats.jld2", "SCBurdenHSC_CID","SCBurdenHSCMean", "SCBurdenHSCVar", "freqs_f", "nVHSC_f", "sampleSize")
    global paramsData = Dict(
        "mutMean" => SCBurdenMean,
        "mutVar" => SCBurdenVar,
        "sample size" => sampleSize,
        "t" => 59
    )
    global freqsData_f = freqsData_f
    global nVData_f = nVData_f
    global nVData_cid = nVData_cid
end

## ========================== Parameters for evolution ===========================
μKnown = nothing
# μKnown = 1.2
paramsIn = Dict{String,Number}(
    "evolve time" => paramsData["t"],
    "sample size" => paramsData["sample size"],
    "N initial" => 1
)
paramsIn["divisions"], paramsIn["μ"] = Theory.getCPRateParamsFromBurdenStats(paramsData["mutMean"], paramsData["mutVar"], μKnown=μKnown)


_mT = [2,7,15]
# paramsIn_mT = Dict[]
# for mT in  _mT
# end

N_N = 1*10^3:1*10^4:7*10^4
p_p = 0.05:0.05:0.95

## ================= call NP space pipeline =================

lCont = 100
NOpt_mT_p = Array{Float64, 2}(undef, length(_mT), lCont)
for (i,mT) in enumerate(_mT)
    params = deepcopy(paramsIn)
    params["mature time"] = mT
    # Theory.extendParamsFromData!(params)
    _N, _p, NOpt_p, vaf1Error_p_N = calcSingleSimNpSpace(params, nVData_f, N_N, p_p, 500, lCont=lCont)
    global _N = _N
    global _p = _p
    global vaf1Error_p_N = vaf1Error_p_N
    NOpt_mT_p[i, :] = NOpt_p
end

##
nVEV_mT_f = Array{Float64, 2}(undef, length(_mT), paramsData["sample size"]+1)
nCEV_MT_m = Vector{Float64}[]

nCP = 200000
cpData_mT_id = Array{Int, 2}(undef, length(_mT), nCP)
let
    pTest = 0.5
    for i in 1:length(_mT)
        NOptSpl_p = Spline1D(_p, NOpt_mT_p[i, :])
        NTest = NOptSpl_p(pTest)
        paramsAlt = copy(paramsIn); paramsAlt["mature time"] = _mT[i]
        paramsEvo = createTestParams(paramsAlt, NTest, pTest)
        vfs = VAFDyn.VFreqspace(Int(round(NTest)), 1001)
        @time VAFDyn.evolveGrowingVAF(vfs, paramsEvo, paramsEvo["evolve time"])
        dfsPDE = VAFDyn.makeDFSfromVFS(vfs, paramsEvo["N final"])
        dfsPDES = VAFDyn.sampler(dfsPDE, paramsEvo["sample size"])
        nVEV_mT_f[i, :] = dfsPDES.n_f
        
        cpData_id = CompoundPoisson.randComPois(Theory.getEffectiveDivisions(paramsEvo), paramsEvo["μ"], nCP)
        cpData_mT_id[i, :] = cpData_id
    end
end

##
let
    λ, μ = Theory.getCPRateParamsFromBurdenStats(paramsData["mutMean"], paramsData["mutVar"])
    global cpData_id = CompoundPoisson.randComPois(λ, μ, nCP)
end

## ========== Plotting ==========
pyplot()

## --------- SC burden ---------
fig1 = histogram(nVData_cid, dpi=400,
    bins=25,
    linewidth=0,
    linealpha=0,
    bar_width=15.0,
    color=:grey60,
    normalize=true,
    label="LS data"
)
stephist!(cpData_mT_id[1,:],
    color=:black,
    normalize=true,
    bins=250,
    linewidth=1.5,
    linestyle=:dash,
    label="Comp Pois"
)
# stephist!(cpData_id,
#     color=:blue,
#     linealpha=0.6,
#     normalize=true,
#     bins=250,
#     linewidth=1.,
#     linestyle=:dash,
#     label="Compound Poisson"
# )
xlabel!("number of mutations")
ylabel!("density of cells")
xlims!(700, 1400)
display(fig1)

# savefig(fig1, "figures/LSData/SCBurden.pdf")
#
# cpData_id = CompoundPoisson.randComPois(Theory.getEffectiveDivisions(paramsTrue), paramsTrue["μ"], 200000)
# stephist!(cpData_id,
#     color=:black,
#     linealpha=0.6,
#     normalize=true,
#     bins=250,
#     linewidth=1.6,
#     linestyle=:dash,
#     label=L"Comp. Poisson: $\mu=1.20$")

## --------- VAF spectrum ---------
fig2 = bar((1:paramsData["sample size"])/paramsData["sample size"], nVData_f[2:end-1],
dpi=400,
yscale=:log10,
# linewidth=0.8,
linewidth=0,
linealpha=0,
color=:grey60,
# fillalpha=0.4,
# fillrange=0
label="LS data"
)
for i in 1:length(_mT)
    plot!(freqsData_f[2:end], nVEV_mT_f[i,2:end],
        label="mature time: "*string(_mT[i])*" years",
        color=2+i,
        linewidth=2
    )
end
xlims!(0, 0.6)
ylims!(8*1E-1, 3E4)
xlabel!("VAF")
ylabel!("number of variants")
display(fig2)
# savefig(fig2, "figures/LSData/VAFSpectrum.pdf")


## --------- N-p space ---------
fig3 = plot(dpi=300)
for i in 1:length(_mT)
    plot!(_p, NOpt_mT_p[i,:],
    color=2+i,
    linewidth=2,
    label="mature time = "*string(_mT[i])*" years",
    )
end

ylims!(10^4,6*10^4)
xlabel!("p")
ylabel!("N")
display(fig3)

# savefig(fig3, "figures/LSData/NPSpace.pdf")

##
fscaler = 0.9
fig4 = plot(fig1, fig2, fig3, layout=(1,3), size=(fscaler*1400,fscaler*400))

display(fig4)
savefig(fig4, "figures/paper/LSdata_analysisPipeline.pdf" )

##

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


## ========== GADFLY Plotting ==========
# fig1 = plot(x=_p, y=NOpt_p, Geom.line, Scale.x_continuous(minvalue=0, maxvalue=3))
# display(fig1)

# fig1 = plot(x=_p, y=_N, z=vaf1Error_p_N, Geom.contour
# # Coord.cartesian(xmin=0, xmax=1, ymin=10^4, ymax=6*10^4)
# # Geom.rectbin,  
# # Scale.x_continuous(minvalue=0, maxvalue=1),
# # Scale.y_continuous(minvalue=1*10^4, maxvalue=6*10^4),
# # Scale.color_continuous()
# )
# # push!(fig1, layer(x=_p, y=NOpt_p, Geom.line))
# display(fig1)

# coord = Coord.cartesian(xmin=0, xmax=1, ymin=10^4, ymax=6*10^4)