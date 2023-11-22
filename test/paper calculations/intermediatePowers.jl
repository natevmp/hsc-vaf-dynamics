include("../../src/vafdyn.jl")
using .VAFDyn
# using DifferentialEquations
# using JLD2, FileIO, Glob

LOADDATA = false
SAVEDATA = false
##

if !LOADDATA
    params = Dict(
        "ρ"=>365/7/4 *1/4,
        "ϕ"=>365/7/4 *3/4,
        "μ"=>1.,
        "N initial"=>1,
        "N final"=>1000,
        "pure births"=>1000,
        "growth rate"=>10.,
        "sample size"=>50,
        # "evolve time"=>60.,
    )
    VAFDyn.completeParams!(params)
    tMat = log(params["N final"])/params["growth rate"]
    println("mature time: ", tMat)

    _tP = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    # _tP = [5, 100]
    _t = _tP .+ tMat

    n_f_T = Vector{Float64}[]
    for (i,t) in enumerate(_t)
        vfs = VAFDyn.VFreqspace(params["N final"], 500)
        # @time VAFDyn.evolveGrowingVAF(vfs, params, t; alg=TRBDF2())
        @time VAFDyn.evolveGrowingVAF(vfs, params, t)
        dfs = VAFDyn.makeDFSfromVFS(vfs, params["N final"])
        # dfs = VAFDyn.DFreqspace(params["N final"])
        # @time VAFDyn.evolveGrowingVAF(dfs, params, t)
        push!(n_f_T, deepcopy(dfs.n_f))
    end
    SAVEDATA && jldsave("data/intermediatePowers.jld2"; _tP, _t, n_f_T, params, tMat)
else
    jldopen("data/intermediatePowers.jld2")
end

##

##
nExpCap(t) = exp(params["growth rate"]*t) < params["N final"] ? exp(params["growth rate"]*t) : params["N final"]
sfsExpT(k) = 2/(k+k^2)
sfsFixedT(k) = 1/k

# n_t = Int.(floor.(nExpCap.(_t)))
n_t = nExpCap.(_t)
nInt_t = Int.(floor.(n_t))
_f_T = [range(0,1, length=n+1) for n in nInt_t]
_f = range(0,1, length=params["N final"]+1)


##
using CairoMakie, LaTeXStrings
using ColorSchemes, Colors
##
fIn = findfirst(_f.>=0.005)
resScale=0.9
cGreys = ColorScheme([RGB{Float64}(i, i, i) for i in range(0.8,0.1,length=length(_t))])
f1 = Figure(
    resolution=(resScale*500,resScale*400),
    fontsize=14,
    font="DejaVu Sans"
)
ax1 = Axis(
    f1[1:2,1:2], xscale=log10, yscale=log10, xlabel="variant allelle frequency", ylabel="rescaled number of variants",
)
hidespines!(ax1, :t, :r)

lines!(_f[fIn:end-1]/2, sfsExpT.(params["N final"]*_f[fIn:end-1])/sfsExpT(params["N final"]*_f[fIn]), linestyle=:solid, label="growing population", color=:purple)
lines!(_f[fIn:end-1]/2, sfsFixedT.(params["N final"]*_f[fIn:end-1])/sfsFixedT(params["N final"]*_f[fIn]), linestyle=:solid, label="fixed population", color=:green)

for i in 1:length(_t)
    lines!(
        _f_T[i][fIn:end-1]/2, n_f_T[i][fIn:nInt_t[i]]/n_f_T[i][fIn],
        # label="t = "*string(_tP[i]),
        color=cGreys[i],
    )
end

# gridLegend = GridLayout()
# f1[1,2] = gridLegend
# gridLegend[1,1] = Legend(f1, ax1, framevisible=false, tellheight=true)
f1[1,2] = Legend(f1, ax1, framevisible=false, tellheight=true)
cbLayout = GridLayout()
f1[2,1] = cbLayout
colsize!(cbLayout, 1, Relative(0.))
Colorbar(
    cbLayout[1,1], limits = (_tP[1]-5,_tP[end]+5), colormap=cgrad(cGreys, 10, categorical = true),
    width=25,
    alignmode = Mixed(right = 15),
    # label="time after maturity",
)
Label(cbLayout[1,2], "years after maturity", rotation = pi/2,alignmode = Mixed(right = 30))
rowsize!(f1.layout, 2, Relative(0.80))
colsize!(f1.layout, 1, Relative(0.45))
CairoMakie.xlims!(ax1, _f[fIn], 1.)
CairoMakie.ylims!(10^-4,1.)
display(f1)
figname ="growthTransition_N"*string(params["N final"])*"_rho"*string(round(params["ρ"],digits=1))*".png"
# save("Figures/GrowthTransition/"*figname, f1)
# save("Figures/GrowthTransition/"*"plot2.pdf", f1)



## ================= Varying population size ========================

params = Dict(
    "ρ"=>1.,
    "ϕ"=>0,
    "μ"=>1.,
    "N initial"=>1,
    "growth rate"=>1.,
)
sampleSize = 500
_N = [600, 1000, 2000]

n_f_N = Vector{Float64}[]
for (i,N) in enumerate(_N)
    params["N final"] = N
    params["pure births"] = params["N final"]
    tMat = log(params["N final"])/params["growth rate"]
    println("mature time: ", tMat)
    tFinal = tMat+20
    dfs = VAFDyn.DFreqspace(params["N final"])
    @time VAFDyn.evolveGrowingVAF(dfs, params, tFinal)
    dfsS = VAFDyn.sampler(dfs, sampleSize)
    push!(n_f_N, deepcopy(dfsS.n_f))
end

##

sfsExpT(k) = 2/(k+k^2)
sfsF2(k) = 1/k^2
sfsFixedT(k) = 1/k
# _f_N = [range(0,1, length=N+1) for N in _N]
_f = range(0,1,length=sampleSize+1)

sfsExpTh_f = 1000*sfsExpT.(range(0,1001,step=1))
sfsExpSTh_f = VAFDyn.sampler(sfsExpTh_f, 1000, sampleSize)
##
using CairoMakie, LaTeXStrings
##
f2 = Figure(resolution=(800,600), fontsize=16)
ax2 = Axis(f2[1,1], xscale=log10, yscale=log10, xlabel="f", ylabel="rescaled number of variants")
# lines!(_f[2:end-1], sfsExpSTh_f[2:end-1]/sfsExpSTh_f[2])
# lines!(_f[2:end-1], sfsF2.(sampleSize*_f[2:end-1]), linestyle=:dash, label="1/f^2")
lines!(_f[2:end-1], sfsExpT.(sampleSize*_f[2:end-1]), linestyle=:dash, label="1/(k+k^2)")
lines!(_f[2:end-1], sfsFixedT.(sampleSize*_f[2:end-1]), linestyle=:dashdot, label="1/k")
for i in 1:length(_N)
    lines!(_f[2:end-1], n_f_N[i][2:end-1]/n_f_N[i][2], label="N = "*string(_N[i]))
end
axislegend(ax2, position=:lb)
display(f2)
