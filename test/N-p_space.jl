include("../src/vafdyn.jl")
using .VAFDyn
include("../src/theory.jl")
using .Theory

using DelimitedFiles, FileIO
using ProgressMeter
using ApproxBayes
using Distributions, Distances

using Plots
gr()

## ========================== Parameters for evolution ===========================
knownMu = false
paramsIn = Dict{String,Number}(
    "evolve time" => paramsData["t"],
    "sample size" => paramsData["sample size"],
    "N initial" => 1,
    "mature time" => 15
    )

if knownMu
    paramsIn["μ"] = 1.2
    paramsIn["division rate"] = let 
        m = paramsData["mutMean"]
        t = paramsData["t"]
        μ = paramsIn["μ"]
        m / (μ*t)
    end
else
    paramsIn["μ"] = paramsData["mutVar"]/paramsData["mutMean"] - 1.
    paramsIn["division rate"] = let 
        m = paramsData["mutMean"]
        v = paramsData["mutVar"]
        t = paramsData["t"]
        m^2/(v-m)/t
    end
end

N_N = 10^4:5*10^3:14*10^4
p_p = 0.05:0.025:0.95


## ============================== create VAF spectra ==============================

function calcExpectedVAFs()

    vaf1_p_N = Array{Float64,2}(undef, length(pVals_p), length(popSizes_N))
    @showprogress for (i,p) in enumerate(pVals_p)
        for (j,N) in enumerate(popSizes_N)

            if knownMu
                λ = mutMean / (evoParams["μ"]*(2-p)*evoParams["t"])
                # λ = 10.
                divRateInferred = meanBurdenS^2/(mutVar-mutMean)/paramsKnown["evolve time"]
            else
                λ = mutMean^2/(mutVar - mutMean) / ( (2-p)*evoParams["t"] )
            end

            # r = params["division rate"]
            # t = params["evolve time"]
            # tf = params["mature time"]
            # Nf = Int(round(N))
            # γ = log(Nf)/tf
            # λ = (r-log(Nf)/t)/(2-p)
            # paramsTest = deepcopy(params)
            # paramsTest["N final"] = Nf
            # paramsTest["p"] = p
            # paramsTest["growth rate"] = γ
            # paramsTest["λ"] = λ

            evoParams["N final"] = N
            evoParams["ρ"] = λ*(1-p)
            evoParams["ϕ"] = λ*p
            evoParams["N initial"] = 1
            evoParams["N mature"] = 15

            vfs = VAFDyn.VFreqspace(evoParams["N"],1001)
            # VAFDyn.evolveVAFfd(vfs, evoParams, evoParams["t"]);
            VAFDyn.evolveVAFfd(vfs, evoParams, evoParams["t"]);
            dfs = VAFDyn.makeDFSfromVFS(vfs, evoParams["N"])
            dfsSampled = VAFDyn.sampler(dfs, evoParams["sampleSize"])

            vaf1_p_N[i,j] = dfsSampled.n_f[2]

        end
    end

    return vaf1_p_N
end

function checkMinmax(λ)

    vaf1InOut_p_N = Array{Float64,2}(undef, 2, 2)
    @showprogress for (i,p) in enumerate( [ pVals_p[1],pVals_p[end] ] )
        for (j,N) in enumerate( [ popSizes_N[1],popSizes_N[end] ] )

            evoParams["ρ"] = λ*(1-p)
            evoParams["ϕ"] = λ*p
            evoParams["N"] = N

            vfs = VAFDyn.VFreqspace(evoParams["N"],1001)
            VAFDyn.evolveVAFfd(vfs, evoParams, evoParams["t"]);
            dfs = VAFDyn.makeDFSfromVFS(vfs, evoParams["N"])
            dfsSampled = VAFDyn.sampler(dfs, evoParams["sampleSize"])

            vaf1InOut_p_N[i,j] = dfsSampled.n_f[2]

        end
    end

    return vaf1InOut_p_N
end

##

vaf1InOut_p_N = checkMinmax(λ)
display(vaf1InOut_p_N)
vaf1InOutError_p_N = vaf1InOut_p_N .- vafData_f[2]
display(vaf1InOutError_p_N)

##

function bestPars(vafError_p_N, NVals_N)

    pLen = size(vafError_p_N)[1]
    vafBestN_p = Vector{Float64}(undef, pLen)
    for i in 1:pLen
        bestNID = argmin(abs.(vafError_p_N[i,:]))
        vafBestN_p[i] = NVals_N[bestNID]
        # println(bestNID)
    end

    return vafBestN_p
end

##
@time vaf1_p_N = calcExpectedVAFs()

##

vaf1Error_p_N = vaf1_p_N .- vafData_f[2]
vafBestNStep_p = bestPars(vaf1Error_p_N, popSizes_N)


## ========== get interpolated solutions ===========
using Dierckx, Roots, Interpolations

# vaf1ErrorInterpol_p_N = interpolate(vaf1Error_p_N, BSpline(Linear()))

vafBestN_p = Vector{Float64}(undef, length(pVals_p))
for i in 1:length(pVals_p)
    vaf1ErrorInterpol_N = Spline1D(popSizes_N, vaf1Error_p_N[i,:])
    vafBestN_p[i] = find_zero(x->vaf1ErrorInterpol_N(x), (popSizes_N[1],popSizes_N[end]), Bisection())
end


## ========== Plotting ==========

μString = string(round(evoParams["μ"],digits=2))

vaf1ErrorSpline_p_N = Spline2D(pVals_p, popSizes_N, vaf1Error_p_N)
pValsCont_p = 0.05:0.01:0.95
NValsCont_N = 10^4:10^3:13*10^4
vaf1ErrorSpline_p_N(0.05, 10000)

vaf1ErrorInterpol_p_N = Array{Float64}(undef, length(pValsCont_p), length(NValsCont_N))
for (i,p) in enumerate(pValsCont_p)
    for (j,N) in enumerate(NValsCont_N)
        vaf1ErrorInterpol_p_N[i,j] = vaf1ErrorSpline_p_N(p, N)
    end
end

# fig = heatmap(transpose(vaf1Error_p_N), clims=[-40000.,40000.])
fig = heatmap(pValsCont_p, NValsCont_N, transpose(vaf1ErrorInterpol_p_N),
    fillalpha=0.6)
plot!(pValsCont_p, NValsCont_N, transpose(vaf1ErrorInterpol_p_N))
plot!(pVals_p, vafBestN_p,
    linewidth=2,
    linestyle=:dash,
    label="best fit",
    color=1)
xlims!(0,1)
# ylims!(10^4,2*10^5)
xlabel!("p (fraction of asymmetric divisions)")
ylabel!("N")
# title!("μ = 1.2")
title!("error (number of variants at 1/S)")
# savefig("figures/ABCfitting/N-p_space_mu"*μString*".pdf")
display(fig)

# fig2 = heatmap(pValsCont_p, NValsCont_N, transpose(vaf1ErrorInterpol_p_N),
#     fillalpha=0.5)
# plot!(pVals_p, vafBestN_p,
#     linewidth=2,
#     linestyle=:dash,
#     label="best fit",
#     color=1)
# # plot!(pVals_p, popSizes_N, transpose(vaf1Error_p_N))
# xlims!(0,1)
# # ylims!(10^4,2*10^5)
# xlabel!("p (fraction of asymmetric divisions)")
# ylabel!("N")
# # title!("μ = 1.2")
# title!("error (number of variants at 1/S)")
# # savefig("figures/ABCfitting/N-p_space_heatmap_mu"*μString*".pdf")
# display(fig2)

##
fig = heatmap(pValsCont_p, NValsCont_N, transpose(abs.(vaf1ErrorInterpol_p_N)),
fillalpha=0.6)
plot!(pValsCont_p, NValsCont_N, transpose(abs.(vaf1ErrorInterpol_p_N)))
plot!(pVals_p, vafBestN_p,
linewidth=2,
linestyle=:dash,
label="best fit",
color=1)
xlims!(0,1)
# ylims!(10^4,2*10^5)
xlabel!("p (fraction of asymmetric divisions)")
ylabel!("N")
title!("absolute value of error (number of variants at 1/S)")
# title!("μ = 1.2")
# savefig("figures/ABCfitting/N-p_spaceAbs_mu"*μString*".pdf")
display(fig)

# fig2 = heatmap(pVals_p, popSizes_N, transpose(abs.(vaf1Error_p_N)))
# plot!(pVals_p, vafBestN_p,
#     linewidth=2,
#     linestyle=:dash,
#     label="data best fit",
#     color=1)
# # plot!(pVals_p, popSizes_N, transpose(vaf1Error_p_N))
# xlims!(0,1)
# # ylims!(10^4,2*10^5)
# xlabel!("p (fraction of asymmetric divisions)")
# ylabel!("N")
# title!("absolute value of error (number of variants at 1/S)")
# # title!("μ = 1.2")
# savefig("figures/ABCfitting/N-p_spaceAbs_heatmap_mu"*μString*".pdf")
# display(fig2)