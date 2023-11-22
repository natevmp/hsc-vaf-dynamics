## ====== Load Data ======
using JLD2
@load "data/Simulations/npSpaceDistribution_N10000.jld2" _p NOpt_sim_p paramsTrue _N

## ====== Get histograms for each point ======
using StatsBase

_NDist = 1e3:1e3:2e4
NOptDistribution_N_p = Array{Float64,2}(undef, length(_NDist)-1, length(_p))
for (i,p) in enumerate(_p)
    h = fit(Histogram, NOpt_sim_p[:,i], _NDist, closed=:right)
    h = StatsBase.normalize(h, mode=:density)
    NOptDistribution_N_p[:, i] = h.weights
end


## ====== Plot ======
using Plots, LaTeXStrings, Colors, ColorSchemes
pyplot()
theme(:default,
    minorgrid=false, 
    gridstyle=:dash, 
    fontfamily="DejaVu Sans", 
    legendfontsize=9, 
    size=(500,400),
)

mygrays = ColorScheme(append!([RGB{Float64}(1, 1, 1),RGB{Float64}(1, 1, 1),RGB{Float64}(1, 1, 1),RGB{Float64}(1, 1, 1)],[RGB{Float64}(i, i, i) for i in 1:-0.01:0]))

colorgrad = cgrad(:grayC, rev = false, alpha = nothing, scale = nothing, categorical = nothing)
# colorgrad = cgrad(mygrays, rev = false, alpha = nothing, scale = nothing, categorical = nothing)


fig = contour(
    _p, _NDist[2:end], NOptDistribution_N_p,
    fill=true,
    c=colorgrad,
    title = "d)",
    titleloc = :left,
    titlefont=font(20, "DejaVu Sans"),
)
scatter!(
    [paramsTrue["p"],], [paramsTrue["N final"],],
    markersize=8,
    markershape=:diamond,
    color=:white,
    label="True value",
)
xlabel!(L"Asymmetric divisions fraction $p$")
ylabel!(L"Population size $N$")
display(fig)
filename = "3d.pdf"
savefig(fig, "Figures/Paper/"*filename)