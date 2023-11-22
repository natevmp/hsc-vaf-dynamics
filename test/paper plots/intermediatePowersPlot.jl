
using JLD2
# jldopen("data/intermediatePowers.jld2")
@load "data/intermediatePowers.jld2" _tP _t n_f_T params tMat
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


# sizeMilimeters = (90, 75)
# sizePoints = 2.8346456692913 .* sizeMilimeters

f1 = Figure(
    # resolution=(resScale*500,resScale*400),
    resolution=(700,600),
    fontsize=24,
    font="DejaVu Sans"
)

ax1 = Axis(
    f1[1:2,1:2],
    xscale=log10,yscale=log10,
    xlabel="Variant allele frequency",
    ylabel="Rescaled number of variants",
)
hidespines!(ax1, :t, :r)

lines!(_f[fIn:end-1]/2, sfsExpT.(params["N final"]*_f[fIn:end-1])/sfsExpT(params["N final"]*_f[fIn]), linestyle=:solid, label="growing population", color=:purple)
lines!(_f[fIn:end-1]/2, sfsFixedT.(params["N final"]*_f[fIn:end-1])/sfsFixedT(params["N final"]*_f[fIn]), linestyle=:solid, label="fixed population", color=:green)
for i in 1:length(_t)
    lines!(
        _f_T[i][fIn:end-1]/2, n_f_T[i][fIn:nInt_t[i]]/n_f_T[i][fIn],
        color=cGreys[i],
    )
end
text!(
    _f[end]/2,
    sfsExpT.(params["N final"]*_f[end-1])/sfsExpT(params["N final"]*_f[fIn]),
    text=L"$f^{-2}$",
    align = (:left, :center),
    textsize=30,
    markerspace = :pixel,
)
text!(
    _f[end]/2,
    sfsFixedT.(params["N final"]*_f[end-1])/sfsFixedT(params["N final"]*_f[fIn]),
    textsize=30,
    text=L"$f^{-1}$",
    align = (:left, :center),
)

f1[1,2] = Legend(f1, ax1, framevisible=false, tellheight=true)
cbLayout = GridLayout()
f1[2,1] = cbLayout
colsize!(cbLayout, 1, Relative(0.))
Colorbar(
    cbLayout[1,1], limits = (_tP[1]-5,_tP[end]+5), colormap=cgrad(cGreys, 10, categorical = true),
    width=25,
    alignmode = Mixed(right = 30),
)
Label(cbLayout[1,2], "years after maturity", rotation = pi/2,alignmode = Mixed(right = 70))
rowsize!(f1.layout, 2, Relative(0.80))
colsize!(f1.layout, 1, Relative(0.45))
CairoMakie.xlims!(ax1, _f[fIn]/2, 1.)
# CairoMakie.ylims!(10^-4,1.)
CairoMakie.ylims!(10^-(5.5),1.)

Label(
    f1[1, 1, TopLeft()], "a)",
    textsize = 40,
    font = "DejaVu Sans",
    padding = (0, 0, 10, 0),
    halign = :left,
    valign = :top,
)

display(f1, )

figname ="growthTransition_N"*string(params["N final"])*"_rho"*string(round(params["ρ"],digits=1))*".png"
save("Figures/GrowthTransition/"*figname, f1)

save("Figures/Paper/"*"2a.pdf", f1)