include("../../src/vafdyn.jl")
using .VAFDyn
using CairoMakie

##

measureTime1 = 25
measureTime2 = 50
# logMatFrac = 9//10
lenVfs = 1001
logMatFrac = 99//100
expGrowthRateFromNT(Nf, tF) = log(Nf)/tF
logisticGrowthRateFromNT(Nf, tF) = log(numerator(logMatFrac)*(Nf-1))/tF

sfsExpT(k) = 2/(k+k^2)
sfsFixedT(k) = 1/k

params = Dict(
    "ρ"=>1.,
    "ϕ"=>1.,
    "μ"=>1.2,
    "N initial"=>1,
    "N final"=>10000,
    "pure births"=>0,
    "sample size"=>50,
    "mature time"=>20.,
    "evolve time"=>measureTime2,
)
params["growth rate exp"] = expGrowthRateFromNT(params["N final"], params["mature time"])
params["growth rate logistic"] = logisticGrowthRateFromNT(params["N final"], params["mature time"])

##
paramsExp = deepcopy(params)
paramsExp["growth rate"] = params["growth rate exp"]
vfs = VAFDyn.VFreqspace(paramsExp["N final"], lenVfs)
nExpT1_f = zeros(length(vfs)-2)
_f = similar(nExpT1_f)
for attempt in 1:100
    @time VAFDyn.evolveGrowingVAF(vfs, paramsExp, measureTime1, growthType="cappedExp")
    # dfs = VAFDyn.makeDFSfromVFS(vfs, paramsExp["N final"])
    if vfs.n_f[2] > 0
        nExpT1_f .= vfs.n_f[2:end-1]
        _f .= vfs.freqs_f[2:end-1]
        break
    end
end
vfs = VAFDyn.VFreqspace(paramsExp["N final"], lenVfs)
nExpT2_f = zeros(length(vfs)-2)
for attempt in 1:100
    @time VAFDyn.evolveGrowingVAF(vfs, paramsExp, measureTime2, growthType="cappedExp")
    # dfs = VAFDyn.makeDFSfromVFS(vfs, paramsExp["N final"])
    if vfs.n_f[2] > 0
        nExpT2_f .= vfs.n_f[2:end-1]
        _f .= vfs.freqs_f[2:end-1]
        break
    end
end
paramsLog = deepcopy(params)
paramsLog["growth rate"] = params["growth rate logistic"]
vfs = VAFDyn.VFreqspace(paramsLog["N final"], lenVfs)
nLogT1_f = zeros(length(vfs)-2)
_f = similar(nExpT2_f)
for attempt in 1:100
    @time VAFDyn.evolveGrowingVAF(vfs, paramsLog, measureTime1, growthType="logistic")
    # dfs = VAFDyn.makeDFSfromVFS(vfs, paramsLog["N final"])
    if vfs.n_f[2] > 0
        nLogT1_f .= vfs.n_f[2:end-1]
        _f .= vfs.freqs_f[2:end-1]
        break
    end
end
vfs = VAFDyn.VFreqspace(paramsLog["N final"], lenVfs)
nLogT2_f = zeros(length(vfs)-2)
_f = similar(nExpT2_f)
for attempt in 1:100
    @time VAFDyn.evolveGrowingVAF(vfs, paramsLog, measureTime2, growthType="logistic")
    # dfs = VAFDyn.makeDFSfromVFS(vfs, paramsLog["N final"])
    if vfs.n_f[2] > 0
        nLogT2_f .= vfs.n_f[2:end-1]
        _f .= vfs.freqs_f[2:end-1]
        break
    end
end

##

using JLD2

jldsave("data/expVsLogisticGrowth.jld2"; params, paramsExp, paramsLog, _f, nExpT1_f ,nExpT2_f, nLogT1_f, nLogT2_f, measureTime1, measureTime2)


##
resScale=1.0
fig1 = Figure(
    resolution=(resScale*800,resScale*600),
    fontsize=24,
    font="DejaVu Sans"
)
Axis(fig1[1,1], xlabel="frequency", ylabel="number of variants", xscale=log10, yscale=log10)
lines!(_f[1:end-50], nExpT1_f[1:end-50] ./ nExp_f[1], label="exponential")
lines!(_f[1:end-50], nLogT1_f[1:end-50] ./ nLog_f[1], label="logistic")
lines!(_f, sfsFixedT.(params["N final"]*_f)/sfsFixedT(params["N final"]*_f[1]), linestyle=:solid, color=:green, label="theory fixed")
lines!(_f, sfsExpT.(params["N final"]*_f)/sfsExpT(params["N final"]*_f[1]), linestyle=:solid, color=:purple, label="theory growth")
axislegend(position=:lb)
display(fig1)

##

fig3  = Figure(
    resolution=(1400, 600),
    fontsize=24,
    font="DejaVu Sans"
)

# --------- growth plot ---------
ax1 = Axis(fig3[1,1], xlabel="time", ylabel="population size")
_t = range(0, params["evolve time"]+1,length=100)
lines!(
    _t, (t->VAFDyn.cappedExponentialGrowth(1, params["N final"], params["growth rate exp"], t)).(_t),
    label="exponential maturation",
    color=:black,
    linestyle=:dot,
    linewidth=1.5,
)
lines!(
    _t, (t->VAFDyn.logisticGrowth(params["N final"], params["growth rate logistic"], t)).(_t),
    label="logistic maturation",
    color=:black,
    linestyle=:dash,
    # linewidth=2,
)
text!(
    [measureTime1,measureTime2],[params["N final"],params["N final"]];
    text=["t₁", "t₂"],
    align=(:left, :top),
    offset=(7,-5),
    textsize=24,)
# text!(measureTime2, params["N final"]; text="t₂", align=(:left, :top), offset=(7,-5))
vlines!(ax1, measureTime1, color=Cycled(1))
vlines!(ax1, measureTime2, color=Cycled(2))
axislegend(position=:rb)
xlims!(0,params["evolve time"])
ylims!(0,nothing)
#fig label

# --------- VAF spectrum plot ---------
ax2 = Axis(fig3[1,2], xlabel="variant allele frequency", ylabel="rescaled number of variants", xscale=log10, yscale=log10)
# measure time 1
lines!(
    _f[1:end-200], nExpT1_f[1:end-200] ./ nExpT1_f[1], 
    color=Cycled(1),
    label="exponential: t₁ = $measureTime1",
    linestyle=:dot,
    linewidth=2.5,
)
lines!(
    _f[1:end-200], nLogT1_f[1:end-200] ./ nLogT1_f[1], 
    color=Cycled(1),
    label="logistic: t₁ = $measureTime1",
    linestyle=:dash,
    linewidth=2,
)
# measure time 2
lines!(
    _f[1:end-200], nExpT2_f[1:end-200] ./ nExpT2_f[1], 
    color=Cycled(2),
    label="exponential: t₂ = $measureTime2",
    linestyle=:dot,
    linewidth=2.5,
)
lines!(
    _f[1:end-200], nLogT2_f[1:end-200] ./ nLogT2_f[1], 
    color=Cycled(2),
    label="logistic: t₂ = $measureTime2",
    linestyle=:dash,
    linewidth=2,
)
# theory predictions
lines!(_f, sfsFixedT.(params["N final"]*_f)/sfsFixedT(params["N final"]*_f[1]), linestyle=:solid, color=:green, label="theory fixed")
lines!(_f, sfsExpT.(params["N final"]*_f)/sfsExpT(params["N final"]*_f[1]), linestyle=:solid, color=:purple, label="theory growth")

axislegend(position=:lb)

Label(
    fig3[1, 1], "a)",
    textsize = 38,
    font = "TeX Gyre Heros Bold",
    padding = (10, 0, 0, 10),
    tellheight = false,
    tellwidth = false,
    halign = :left,
    valign = :top,
)

Label(
    fig3[1, 2], "b)",
    textsize = 38,
    font = "TeX Gyre Heros Bold",
    padding = (10, 0, 0, 10),
    tellheight = false,
    tellwidth = false,
    halign = :left,
    valign = :top,
)

display(fig3)
##
figname ="S1-ExpVsLogisticGrowth.pdf"
save("../Figures/GrowthTransition/"*figname, fig3)
