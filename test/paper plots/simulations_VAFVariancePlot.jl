using JLD2, Glob, DataFrames, FileIO
using Statistics
using DataFrames

##

_fNames = glob("data/Simulations/Nf10000/singlePatientFullSim_Ni1_*.jld2")
paramsTrue, nVarSim_f, nVarSimS_f = load(_fNames[1], "paramsTrue", "nVarSim_f", "nVarSimS_f")
lenF = length(nVarSim_f)
n_sim_f = Array{Int,2}(undef, length(_fNames), lenF)
for (i,fName) in enumerate(_fNames)
    paramsTrue, nVarSim_f, nVarSimS_f = load(fName, "paramsTrue", "nVarSim_f", "nVarSimS_f")
    n_sim_f[i,:] = nVarSim_f
end
nAv_f = vec(mean(n_sim_f,dims=1))
nStd_f = vec(std(n_sim_f,dims=1))

##
using CairoMakie, LaTeXStrings
using ColorSchemes, Colors
resScale=0.9
##


function smoother(u_x, w=9; d::Int=0)
    if w%2==0
        println("error: please take uneven value for windowsize")
    end
    wSide = Int(floor(w/2))
    s_x = Vector{Float64}(undef, length(u_x))
    s_x[1:d+wSide] .= u_x[1:d+wSide]
    for i in eachindex(u_x)[d+wSide+1:end-wSide]
        s_x[i] = mean(u_x[i-wSide:i+wSide])
    end
    s_x[end-wSide:end] .= u_x[end-wSide:end]
    return s_x
end

# smoother(nAv_f .+ nStd_f)


##


fig1 = Figure(
    resolution=(resScale*500,resScale*400),
    fontsize=14,
    font="DejaVu Sans"
)
ax1 = Axis(
    fig1[1,1],
    yscale=log10,
    xscale=log10,
    xlabel="variant allele frequency",
    ylabel="number of variants",
)

lines!(
    range(0,1,length=lenF), nAv_f .+ nStd_f,
    linestyle=:dash,
    label=L"V_f \pm \sigma_f",
    color=Cycled(2),
)
lines!(
    range(0,1,length=lenF), (x-> x>0 ? x : 1E-10).(nAv_f .- nStd_f),
    linestyle=:dash,
    color=Cycled(2),
)
lines!(
    range(0,1,length=lenF), nAv_f,
    label=L"V_f",
    color=Cycled(1),
)
xlims!(1/lenF,1E-1)
ylims!(0.1, maximum(nAv_f))
axislegend(position=:rt)
display(fig1)

##

fig2 = Figure(
    resolution=(resScale*500,resScale*400),
    fontsize=14,
    font="DejaVu Sans"
)

ax2 = Axis(
    fig2[1,1],
    yscale=log10,
    xscale=log10,
    xlabel="variant allele frequency",
    ylabel="relative standard deviation",
)

lines!(
    range(0,1,length=lenF)[nAv_f.!=0], nStd_f[nAv_f.!=0] ./ nAv_f[nAv_f.!=0],
    color=:grey15,
    # label=L"\sigma_f / V_f",
)

xlims!(1/lenF,1E-1)
ylims!(1E-2, 1E2)

display(fig2)

##



fig3 = Figure(
    resolution=(resScale*950,resScale*400),
    fontsize=14,
    font="DejaVu Sans"
)

ax1 = Axis(
    fig3[1,1],
    yscale=log10,
    xscale=log10,
    xlabel="variant allele frequency",
    ylabel="number of variants",
)
# lines!(
#     # range(0,1,length=lenF), nAv_f .+ nStd_f,
#     range(0,1,length=lenF), smoother(nAv_f .+ nStd_f, 15, d=20),
#     linestyle=:dash,
#     label=L"V_f \pm \sigma_{f}",
#     color=Cycled(2),
# )
# lines!(
#     # range(0,1,length=lenF), (x-> x>0 ? x : 1E-10).(nAv_f .- nStd_f),
#     range(0,1,length=lenF), (x-> x>0 ? x : 1E-10).(smoother(nAv_f .- nStd_f, 15; d=20)),
#     linestyle=:dash,
#     color=Cycled(2),
# )
band!(
    range(0,1,length=lenF),
    (x->x>0 ? x : 1E-10).(smoother(nAv_f .+ nStd_f, 9; d=20)),
    (x->x>0 ? x : 1E-10).(smoother(nAv_f .- nStd_f, 5; d=20)),
    color=Cycled(2),
    label=L"V_f \pm \sigma_{f}",
)
lines!(
    range(0,1,length=lenF), nAv_f,
    label=L"V_f",
    color=Cycled(1),
)
xlims!(1/lenF,1E-1)
ylims!(0.1, maximum(nAv_f))
axislegend(position=:rt)

ax2 = Axis(
    fig3[1,2],
    yscale=log10,
    xscale=log10,
    xlabel="variant allele frequency",
    ylabel="relative standard deviation",
)
lines!(
    range(0,1,length=lenF)[nAv_f.!=0], nStd_f[nAv_f.!=0] ./ nAv_f[nAv_f.!=0],
    color=:grey15,
    label=L"\sigma_f / V_f",
)
xlims!(1/lenF,1E-1)
ylims!(1E-2, 1E2)
axislegend(position=:rt)

Label(
    fig3[1, 1], "a)",
    textsize = 24,
    font = "TeX Gyre Heros Bold",
    padding = (10, 0, 0, 10),
    tellheight = false,
    tellwidth = false,
    halign = :left,
    valign = :top,
)

Label(
    fig3[1, 2], "b)",
    textsize = 24,
    font = "TeX Gyre Heros Bold",
    padding = (10, 0, 0, 10),
    tellheight = false,
    tellwidth = false,
    halign = :left,
    valign = :top,
)

display(fig3)

figname ="S2-vafSpectrumStd.pdf"
save("Figures/Paper/"*figname, fig3)