using JLD2, Glob
using Plots, LaTeXStrings


## ================== Load data ====================

@load "data/LSData_TmNhSpace/LSData_TmNhSpaceInference_lN10.jld2" NmaxInterpol_tM_nh _tM _Nh NmaxOpt_tM_nh paramsTot

_tMC = range(_tM[1], _tM[end], length=size(NmaxInterpol_tM_nh,1))
_NhC = range(_Nh[1], _Nh[end], length=size(NmaxInterpol_tM_nh,2))

## ================== Plotting ====================

pyplot()
theme(:default,
    minorgrid=false,
    gridstyle=:dash,
    fontfamily="DejaVu Sans",
    showaxis=true,
    gridlinewidth=0.7,
    # size=(0.9*500,0.9*400),
)

##

colorgrad = cgrad(:grayC, rev = false, alpha = nothing, scale = nothing, categorical = nothing)

fig1 = contour(
    _tMC, _NhC, transpose(NmaxInterpol_tM_nh),
    c=colorgrad,
    fill=true,
    colorbar_title=L"Maximum population size at maturity $N$",
    size=(0.9*500,0.9*350),
    title = "d)",
    titleloc = :left,
    titlefont=font(20, "DejaVu Sans"),
)
xlims!(2,10)
xlabel!(L"Maturation time $t_M$ (years)")
ylabel!(L"Population size at maintenance $N_H$")
display(fig1)

savefig(fig1, "Figures/Paper/4d.pdf")