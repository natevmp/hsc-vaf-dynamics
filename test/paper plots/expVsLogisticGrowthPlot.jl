using CairoMakie
resScale=0.9
##
using JLD2
@load "data/expVsLogisticGrowth.jld2"

##

expGrowthRateFromNT(Nf, tF) = log(Nf)/tF
logisticGrowthRateFromNT(Nf, tF) = log(numerator(logMatFrac)*(Nf-1))/tF

function cappedExponentialGrowth(Ni, K, r, t)
	return Ni*exp(r*t)<K ? Ni*exp(r*t) : K
end
function logisticGrowth(K, r, t)
    K / ( 1 + (K-1)*exp(-r*t) )
end

sfsExpT(k) = 2/(k+k^2)
sfsFixedT(k) = 1/k

##

fig3  = Figure(
    resolution=(950, 400),
    fontsize=14,
    font="DejaVu Sans"
)

# --------- growth plot ---------
ax1 = Axis(fig3[1,1], xlabel="time", ylabel="population size")
_t = range(0, params["evolve time"]+1,length=100)
lines!(
    _t, (t->cappedExponentialGrowth(1, params["N final"], params["growth rate exp"], t)).(_t),
    label="exponential maturation",
    color=:black,
    linestyle=:dot,
    linewidth=1.4,
)
lines!(
    _t, (t->logisticGrowth(params["N final"], params["growth rate logistic"], t)).(_t),
    label="logistic maturation",
    color=:black,
    linestyle=:dash,
    linewidth=1.2,
)
text!(
    [measureTime1,measureTime2],[params["N final"],params["N final"]];
    text=["t₁", "t₂"],
    align=(:left, :top),
    offset=(7,-5),
    textsize=14,
)
# text!(measureTime2, params["N final"]; text="t₂", align=(:left, :top), offset=(7,-5))
vlines!(ax1, measureTime1, color=Cycled(1))
vlines!(ax1, measureTime2, color=Cycled(2))
axislegend(position=:rb)
xlims!(0,params["evolve time"]+5)
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
    linewidth=1.4,
)
lines!(
    _f[1:end-200], nLogT1_f[1:end-200] ./ nLogT1_f[1], 
    color=Cycled(1),
    label="logistic: t₁ = $measureTime1",
    linestyle=:dash,
    linewidth=1.2,
)
# measure time 2
lines!(
    _f[1:end-200], nExpT2_f[1:end-200] ./ nExpT2_f[1], 
    color=Cycled(2),
    label="exponential: t₂ = $measureTime2",
    linestyle=:dot,
    linewidth=1.4,
)
lines!(
    _f[1:end-200], nLogT2_f[1:end-200] ./ nLogT2_f[1], 
    color=Cycled(2),
    label="logistic: t₂ = $measureTime2",
    linestyle=:dash,
    linewidth=1.2,
)
# theory predictions
lines!(_f, sfsFixedT.(params["N final"]*_f)/sfsFixedT(params["N final"]*_f[1]), linestyle=:solid, color=:green, label="theory fixed")
lines!(_f, sfsExpT.(params["N final"]*_f)/sfsExpT(params["N final"]*_f[1]), linestyle=:solid, color=:purple, label="theory growth")

axislegend(position=:lb)

Label(
    fig3[1, 1], "A",
    textsize = 24,
    font = "TeX Gyre Heros Bold",
    padding = (10, 0, 0, 10),
    tellheight = false,
    tellwidth = false,
    halign = :left,
    valign = :top,
)

Label(
    fig3[1, 2], "B",
    textsize = 24,
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
save("Figures/Paper/"*figname, fig3)
