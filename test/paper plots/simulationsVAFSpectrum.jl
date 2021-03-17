## Load packages
# using Gadfly
using Plots
using JLD2, FileIO
using Glob
using Loess
include("../../src/vafdyn.jl")
using .VAFDyn
include("../../src/theory.jl")
using .Theory


## =============== Load data ===============
# --------- sim VAF data ---------

fileNames_ = glob("singlePatient*.jld2", "./data/SinglePatientPipeline/Nf10000")
nSims = length(fileNames_)
@load fileNames_[1] paramsTrue timesSim_ nCellSim_t
nV_sim_f = Array{Float64, 2}(undef, nSims, paramsTrue["N final"]+1)
nVS_sim_f = Array{Float64, 2}(undef, nSims, paramsTrue["sample size"]+1)
for (i,fName) in enumerate(fileNames_)
    @load fName nVarSim_f nVarSimS_f
    nV_sim_f[i, :] .= nVarSim_f
    nVS_sim_f[i, :] .= nVarSimS_f
end

# --------- pred VAF data ---------
filenameVAFDyn = "data/SinglePatientPipeline/Nf"*string(paramsTrue["N final"])*"/vafDynGrowth_N"*string(paramsTrue["N final"])*".jld2"
@load filenameVAFDyn dfs dfsS dfsPDE dfsPDES

# Get averaged vaf spectra
nVAv_f = sum(nV_sim_f, dims=1)/nSims
nVSAv_f = sum(nVS_sim_f, dims=1)/nSims


## ======================= Plotting =======================
pyplot(legendfontsize=10, guidefontsize=12, tickfontsize=10,  size=(500,400))
simTestID = 3

freqs_f = (0:paramsTrue["N final"]) / paramsTrue["N final"]
freqsS_f = (0:paramsTrue["sample size"]) / paramsTrue["sample size"]

## ---------- complete ---------
fig1 = plot(yscale=:log10)
plot!(freqs_f[2:end], nV_sim_f[1, 2:end],
# linewidth=0,
linetype=:steppre,
linewidth=0.2,
linealpha=0.9,
fillrange=0,
fillalpha=0.9,
fillcolor=:match,
label="single simulation"
)
# for sim in 2:5
#     plot!(freqs_f[2:end], nV_sim_f[sim, 2:end],
#     # linewidth=0,
#     linewidth=0.2,
#     linealpha=0.5,
#     fillalpha=0.5,
#     fillrange=0,
#     fillcolor=:match,
#     label=""
#     )
# end
# model = loess(freqs_f[2:end], nVAv_f[2:end])
# vs = predict(model, freqs_f[2:end])
# plot!(freqs_f[2:end], vs,
plot!(freqs_f[2:end], nVAv_f[2:end],
color=:grey35,
linealpha=1,
label="simulations average")
plot!(freqs_f[2:end], dfs.n_f[2:end],
color=:black,
linestyle=:dash,
label="predicted average"
)
# xlims!(0,0.25)
ylims!(10^0,10^5)
xlabel!("VAF")
ylabel!("number of variants")
xlims!(1/paramsTrue["N final"],1)
title!("complete")
display(fig1)

## ---------- sample ---------
fig2 = plot(yscale=:log10)
plot!(freqsS_f[2:end], nVS_sim_f[1, 2:end],
linetype=:steppre,
color=2,
linewidth=0.2,
linealpha=0.9,
fillrange=0,
fillalpha=0.9,
fillcolor=:match,
label="single simulation"
)
# for sim in 2:5
#     bar!(freqsS_f[2:end], nVS_sim_f[sim, 2:end],
#     linewidth=0,
#     fillalpha=0.5,
#     label=""
#     )
# end
plot!(freqsS_f[2:end], nVSAv_f[2:end],
color=:grey35,
linealpha=1,
label="simulations average")
plot!(freqsS_f[2:end], dfsS.n_f[2:end],
color=:black,
linestyle=:dash,
label="predicted average")
ylims!(10^0,10^5)
xlims!(1/paramsTrue["sample size"],1)
xlabel!("VAF")
ylabel!("number of variants")
title!("sample")
display(fig2)

## ------- save fig --------
fig3 = plot(fig1, fig2, layout=2, size=(900, 400))
display(fig3)
savefig(fig3, "figures/paper/simsVAFSpectrumTrueVsSample.pdf")






## ===== GADFLY TESTING =====
# fig1data = layer(x=freqs_f[2:end], y=nV_sim_f[simTestID, 2:end],
# Geom.bar,
# Theme(default_color=color("deepskyblue"))
# )
# fig1dataAv = layer(x=freqs_f[2:end], y=nVAv_f[2:end],
# Geom.smooth(method=:loess,smoothing=0.001),
# Theme(default_color=color("grey40"))
# )
# fig1ev = layer(x=freqs_f[2:end], y=dfs.n_f[2:end],
# Geom.line,
# Theme(default_color=color("black"), line_style=[:dash])
# )

# fig1 = plot(fig1ev, fig1dataAv, fig1data,
# Scale.y_log10,
# Coord.cartesian(xmin=0, xmax=1, ymin=0.1, ymax=5),
# Guide.xlabel("VAF"),
# Guide.ylabel("number of variants"),
# Guide.manual_color_key("", ["I'm blue l1", "I'm red l2", "I'm green l3"], ["deepskyblue", "grey40", "black"])
# )
# display(fig1)