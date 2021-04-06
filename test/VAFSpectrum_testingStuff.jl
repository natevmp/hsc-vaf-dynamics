##
using Plots
pyplot()
using JLD2, FileIO
using Glob
using Loess
##
include("../src/vafdyn.jl")
using .VAFDyn
# include("../src/vafSim.jl")
# using .VAFSim
include("../src/theory.jl")
using .Theory


##

LOADDATA = true


if !LOADDATA

    params = Dict(
        "N initial" => 1,
        "N final" => 1000,
        "μ" => 1.2,
        "λ" => 4.2,
        "p" => 0.5,
        "growth rate" => 0.3,
        "exp growth rate" => 0.3,
        "lin growth rate" => 0.01,
        "sample size" => 89,
        "evolve time" => 59
    )
    Theory.extendParams!(params)

    # Multiple Single Cell simulations

    nSims = 1
    @time times_t, nLive_t, vaf_n, vafS_n, burden_m, burdenB_m = VAFSim.birthDeathFixedGrowth(params, params["evolve time"], 0.1, showprogress=true)

    nV_Sim_f = Vector{Int64}[]
    nVS_Sim_f = Vector{Int64}[]
    for i in 1:nSims
        @time times_t, nLive_t, nV_f, nVS_f, burden_m, burdenS_m = VAFSim.birthDeathFixedGrowth(params, params["evolve time"], 0.1, showprogress=true)
        push!(nV_Sim_f, nV_f)
        push!(nVS_Sim_f, nVS_f)
    end

    # Save data

    name = "./data/vafSim_"*string(nSims)*"sims_Ni"*string(params["N initial"])*"_Nf"*string(params["N final"])*".jld2"
    @save name nV_Sim_f nVS_Sim_f

else
    ## Load data
    fileNames_ = glob("singlePatient*.jld2", "./data/SinglePatientPipeline/Nf10000")
    nSims = length(fileNames_)
    @load fileNames_[1] paramsTrue timesSim_ nCellSim_t

    # nV_Sim_f = Vector{Int64}[]
    # nVS_Sim_f = Vector{Int64}[]
    nV_sim_f = Array{Float64, 2}(undef, nSims, paramsTrue["N final"]+1)
    nVS_sim_f = Array{Float64, 2}(undef, nSims, paramsTrue["sample size"]+1)
    for (i,fName) in enumerate(fileNames_)
        @load fName nVarSim_f nVarSimS_f
        # push!(nV_Sim_f, nVarSim_f)
        # push!(nVS_Sim_f, nVarSimS_f)
        nV_sim_f[i, :] .= nVarSim_f
        nVS_sim_f[i, :] .= nVarSimS_f
    end

end

## VAF spectrum Growth ODE Evolve

vfs = VAFDyn.VFreqspace(paramsTrue["N final"], 200)
@time VAFDyn.evolveGrowingVAF(vfs, paramsTrue, paramsTrue["evolve time"])
dfsPDE = VAFDyn.makeDFSfromVFS(vfs, paramsTrue["N final"])
dfsPDES = VAFDyn.sampler(dfsPDE, paramsTrue["sample size"])

dfsGrowth = VAFDyn.DFreqspace(paramsTrue["N final"])
@time VAFDyn.evolveGrowingVAF(dfsGrowth, paramsTrue, paramsTrue["evolve time"])
dfsGrowthS = VAFDyn.sampler(dfsGrowth, paramsTrue["sample size"])

filename = "data/SinglePatientPipeline/vafDynGrowth_N"*string(paramsTrue["N final"])*".jld2"
save(filename, Dict(
    "dfs" => dfsGrowth,
    "dfsS" => dfsGrowthS,
    "dfsPDE" => dfsPDE,
    "dfsPDES" => dfsPDES
))

##
# p0a = Plots.plot(dfsPDE.freqs_f[2:end-1], dfsPDE.n_f[2:end-1],
# yscale=:log10, 
#     linewidth=1.5,
#     label="PDE")
# plot!(dfsGrowth.freqs_f[2:end-1], dfsGrowth.n_f[2:end-1], 
#     linewidth=1.5,
#     linestyle=:dash,
#     label="MC")
# plot!(dfsPDES.freqs_f[2:end-1], dfsPDES.n_f[2:end-1], 
#     linewidth=1.5,
#     label="PDE S")
# plot!(dfsGrowthS.freqs_f[2:end-1], dfsGrowthS.n_f[2:end-1],
#     linewidth=1.5,
#     linestyle=:dash,
#     label="MC S")

# display(p0a)
# p0b = Plots.plot(dfsPDES.freqs_f[2:end-1], dfsPDES.n_f[2:end-1], yscale=:log10, label="PDE S")
# plot!(dfsGrowthS.freqs_f[2:end-1], dfsGrowthS.n_f[2:end-1], label="MC S")

# p0 = Plots.plot(p0a, p0b, layout=2)



## Get averaged vaf spectrum

nVAv_f = sum(nV_sim_f, dims=1)/nSims
nVSAv_f = sum(nVS_sim_f, dims=1)/nSims

## VAF Spectrum ODE evolve
timeODE1 = paramsTrue["evolve time"]
vfs1 = VAFDyn.VFreqspace(paramsTrue["N"], 500)
@time VAFDyn.evolveVAF(vfs1, paramsTrue, timeODE1)
dfs1 = VAFDyn.makeDFSfromVFS(vfs1, paramsTrue["N"])
dfs1S = VAFDyn.sampler(dfs1, paramsTrue["sample size"])

timeODE2 = paramsTrue["evolve time"] - paramsTrue["mature time"]
vfs2 = VAFDyn.VFreqspace(paramsTrue["N"], 500)
@time VAFDyn.evolveVAF(vfs2, paramsTrue, timeODE2)
dfs2 = VAFDyn.makeDFSfromVFS(vfs2, paramsTrue["N"])
dfs2S = VAFDyn.sampler(dfs2, paramsTrue["sample size"])

# vfs = VAFDyn.VFreqspace(params["N final"], 1000)
# @time VAFDyn.evolveGrowingVAF(vfs, params, params["evolve time"])
# dfs = VAFDyn.makeDFSfromVFS(vfs, params["N final"])
# dfsS = VAFDyn.sampler(dfs, params["sample size"])

## Plot pop size over time
# p1 = plot(times_t, nLive_t,
#     label="",
#     linewidth=1.5)
# xlabel!("time")
# xlabel!("N")
# display(p1)

## ====== Get live VAF spectrum =======

# freqs_f = (0:paramsTrue["N final"]) / paramsTrue["N final"]
# liveNv_Sim_f = Vector{Int64}[]
# liveFreqs_Sim_f = Vector{Float64}[]
# freqsS_f = (0:paramsTrue["sample size"]) / paramsTrue["sample size"]
# liveNvS_Sim_f = Vector{Int64}[]
# liveFreqsS_Sim_f = Vector{Float64}[]
# for sim in 1:nSims
#     Nv_f = nV_Sim_f[sim]
#     liveNv_f = Int64[]
#     liveFreqs_f = Float64[]
#     NvS_f = nVS_Sim_f[sim]
#     liveNvS_f = Int64[]
#     liveFreqsS_f = Float64[]
#     for i in 1:length(Nv_f)
#         if Nv_f[i] != 0
#             push!(liveNv_f, Nv_f[i])
#             push!(liveFreqs_f, freqs_f[i])
#         end
#     end
#     for i in 1:length(NvS_f)
#         if NvS_f[i] != 0
#             push!(liveNvS_f, NvS_f[i])
#             push!(liveFreqsS_f, freqsS_f[i])
#         end
#     end
#     push!(liveNv_Sim_f, liveNv_f)
#     push!(liveFreqs_Sim_f, liveFreqs_f)
#     push!(liveNvS_Sim_f, liveNvS_f)
#     push!(liveFreqsS_Sim_f, liveFreqsS_f)
# end

# ##
# liveNvAv_f = Float64[]
# liveFreqsAv_f = Float64[]
# for (i,f) in enumerate(freqs_f)

#     nF_sim = [nF_f[i] for nF_f in nV_Sim_f]
#     nFR_sim = nF_sim[nF_sim .!= 0]
#     if length(nFR_sim)>0
#         push!(liveNvAv_f, sum(nFR_sim)/length(nFR_sim) )
#         push!(liveFreqsAv_f, f)
#     end

# end

# liveNvSAv_f = Float64[]
# liveFreqsSAv_f = Float64[]
# for (i,f) in enumerate(freqsS_f)

#     nF_sim = [nF_f[i] for nF_f in nVS_Sim_f]
#     nFR_sim = nF_sim[nF_sim .!= 0]
#     if length(nFR_sim)>0
#         push!(liveNvSAv_f, sum(nFR_sim)/length(nFR_sim) )
#         push!(liveFreqsSAv_f, f)
#     end

# end

## =============== Var calc ===============

# nVVar_f = var.([nV_sim_f[]])


## ======================= Plot VAF spectrum =======================
pyplot(legendfontsize=12,tickfontsize=12)
# gr()
##

freqs_f = (0:paramsTrue["N final"]) / paramsTrue["N final"]
freqsS_f = (0:paramsTrue["sample size"]) / paramsTrue["sample size"]


## ---------- complete ---------
fig1 = plot(yscale=:log10)
plot!(freqs_f[2:end], nV_sim_f[1, 2:end],
# linewidth=0,
size=(600,400),
linewidth=0.5,
# legendfontsize=12,
linealpha=0.9,
fillalpha=0.9,
fillrange=0,
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
plot!(freqs_f[2:end], nVAv_f[2:end],
color=:grey35,
linealpha=1,
label="simulations average")
plot!(freqs_f[2:end], dfsGrowth.n_f[2:end],
color=:black,
linestyle=:dash,
label="predicted average")
ylims!(10^0,10^5)
# xlims!(0,0.25)
xlims!(0,1)
display(fig1)

## ---------- sample ---------
fig2 = plot(yscale=:log10, dpi=300)
bar!(freqsS_f[2:end], nVS_sim_f[1, 2:end],
color=2,
linewidth=0,
fillalpha=0.9,
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
plot!(freqsS_f[2:end], dfsGrowthS.n_f[2:end],
color=:black,
linestyle=:dash,
label="predicted average")
ylims!(10^0,10^5)
xlims!(0,1)
display(fig2)

##




fig1 = plot(yscale=:log10, xscale=:log10)

plot!(freqs_f[2:end], nV_sim_f[simID, 2:end],
color=1,
linewidth=0,
fillrange=0,
fillcolor=:match,
fillalpha=0.5,
label="complete"
)
bar!(freqsS_f[2:end], nVS_sim_f[simID, 2:end],
color=2,
linewidth=0,
fillalpha=0.5,
label="sample"
)
ylims!(10^0,10^5)
xlims!(0.001,1)
display(fig1)

## ---------- True VAF ----------

# colorSim = "steelblue1"
# colorSim = "grey50"
cGrey = "grey55"

# single sim
# p2 = Plots.plot(freqs_f[2:end], vaf_n[2:end], yscale=:log10, dpi=500,
#     color=cGrey,
#     linecolor=cGrey,
#     fillrange=0,
#     # fillalpha=0.5,
#     # linewidth=1.3,
#     # linealpha=0.5,
#     # bar_width=1/paramsTrue["N final"],
#     # linewidth=0,
#     # linealpha=0,
#     # fillalpha=0.5,
#     label="single simulation")
p2 = Plots.bar(freqs_f[2:end], vaf_n[2:end], yscale=:log10, dpi=500,
    color=:black,
    # linecolor=1,
    # linewidth=1.3,
    fillalpha=1,
    linealpha=0,
    # bar_width=1/paramsTrue["N final"],
    # linewidth=0,
    # linealpha=0,
    # fillalpha=0.5,
    label="single simulation")



# averaged sims
plot!(freqs_f[2:end], nVAv_f[2:end],
# plot!(smoothed_f, nVAvSmoothed_f,
    color=1,
    linewidth=0.8,
    linestyle=:solid,
    # linealpha=0.7,
    label="all variants \naveraged simulations")


# let
#     xs = liveFreqsAv_f[2:end-1]
#     ys = liveNvAv_f[2:end-1]
#     model = loess(xs, ys)
#     global smoothed_f = range(extrema(xs)...; step = 0.1)
#     global nVAvSmoothed_f = predict(model, smoothed_f)
# end

# # averaged 'live' variants
# plot!(liveFreqsAv_f, liveNvAv_f,
# # plot!(smoothed_f, nVAvSmoothed_f,
#     # linewidth=0.8,
#     color=2,
#     linewidth=0.8,
#     linestyle=:dash,
#     linealpha=0.7,
#     label="'live' variants \naveraged simulations")

# growth ODE
plot!(dfsGrowth.freqs_f[2:end], dfsGrowth.n_f[2:end],
    color=:black,
    yscale=:log10,
    linewidth=1.2,
    # linestyle=:dash,
    label="model average prediction")

# fixed ODE 2
Plots.plot!(dfs2.freqs_f[2:end], dfs2.n_f[2:end], 
    color=:black,
    linewidth=1.2,
    linestyle=:dash,
    label="growth free model: t="*string(timeODE2))

xlims!(0,0.2)
# xlims!(0,1)
# xlims!(1/paramsTrue["N final"], 1)
ylims!(1*10^0, 3*10^4)
xlabel!("VAF")
ylabel!("number of variants")
title!("N = "*string(paramsTrue["N final"])*
    ", age = "*string(paramsTrue["evolve time"])*
    ", μ="*string(paramsTrue["μ"])*
    ", λ="*string(paramsTrue["λ"])*
    ", p="*string(paramsTrue["p"])
    )
display(p2)


## ---------- Sampled VAF ----------

freqsS_f = (0:paramsTrue["sample size"]) / paramsTrue["sample size"]
# single sim
p3 = bar(freqsS_f[2:end], vafS_n[2:end], yscale=:log10, dpi=500,
    bar_width=1/paramsTrue["sample size"],
    color=cGrey,
    linecolor=cGrey,
    # linewidth=1.5,
    linewidth=0,
    # linealpha=0.5,
    # fillalpha=0.5,
    fillrange=0,
    label="single simulation")

# p3 = plot(freqsS_f[2:end], vafS_n[2:end], yscale=:log10, dpi=500,
#     color=colorSim,
#     fillrange=0,
#     # linecolor=colorSim,
#     # linewidth=1.5,
#     linewidth=0,
#     # linealpha=0.5,
#     # fillalpha=0.5,
#     label="single simulation")

# averaged sims
plot!(freqsS_f[2:end], nVSAv_f[2:end],
    color=2,
    linewidth=1.2,
    linestyle=:solid,
    label="all variants \naveraged simulations")

# # averaged 'live' variants
# plot!(liveFreqsSAv_f, liveNvSAv_f,
#     color=2,
#     linewidth=1.2,
#     linestyle=:dash,
#     label="'live' variants \naveraged simulations")

# growth ODE
# plot!(dfsPDES.freqs_f[2:end], dfsPDES.n_f[2:end],
#     yscale=:log10,
#     linewidth=2,
#     label="PDE growth")
plot!(dfsGrowthS.freqs_f[2:end], dfsGrowthS.n_f[2:end],
    color=:black,
    linewidth=1.2,
    # linestyle=:dash,
    label="model average prediction growth")
# # fixed ODE 1
# plot!(dfs1S.freqs_f[2:end], dfs1S.n_f[2:end],
#     yscale=:log10,
#     linewidth=1.5,
#     label="ODE: t="*string(timeODE1))

# fixed ODE 2
plot!(dfs2S.freqs_f[2:end], dfs2S.n_f[2:end],
    color=:black,
    linewidth=1.2,
    linestyle=:dash,
    label="growth free model: t="*string(timeODE2))
xlims!(0,1)
ylims!(1*10^0, 3*10^4)
xlabel!("VAF")
ylabel!("number of variants")
title!("Sample size = "*string(paramsTrue["sample size"]))
display(p3)


## ---- Combined figure -----
# p4 = plot(p1, p2,p3, layout=3)
p4 = plot(p2,p3, layout=2, size=(1000,400))
savefig(p4, "./figures/vafODESimCompare/vafGrowthSimsLim_N"*string(paramsTrue["N final"])*"_sims"*string(nSims)*".pdf")
display(p4)

##
p5 = p4
xlims!(0,1)
# savefig(p5, "./figures/vafODESimCompare/vafGrowthSimsFull_N"*string(paramsTrue["N final"])*"_sims"*string(nSims)*".pdf")
display(p5)















# ##
# gr()
# ##
# plot6 = bar(freqs_f[2:end], vaf_n[2:end], yscale=:log10, dpi=400,
#     color=1,
#     linecolor=1,
#     linewidth=0.8,
#     linealpha=0.5,
#     fillalpha=0.5,
#     label="single simulation")
# bar!(freqsS_f[2:end], vafS_n[2:end], dpi=400,
#     color=2,
#     linecolor=2,
#     linewidth=0.8,
#     linealpha=0.5,
#     fillalpha=0.5,
#     label="single simulation")


# ylims!(10^0, 3*10^4)
# display(plot6)










# ##

# testSim = 4

# p5 = bar(liveFreqs_Sim_f[testSim], liveNv_Sim_f[testSim],
#     yscale=:log10,
#     linealpha=0,
#     fillalpha=0.5,
#     label="single sim")

# # for sim in 1:5
# #     # bar!(liveFreqs_Sim_f[sim], liveNv_Sim_f[sim],
# #     #     # bar_width=1.,
# #     #     linealpha=0,
# #     #     fillalpha=0.5,
# #     #     label="")
# #     plot!(liveFreqs_Sim_f[sim], liveNv_Sim_f[sim],
# #         # bar_width=1.,
# #         linealpha=0.8,
# #         linewidth=0.5,
# #         label="")
# # end

# plot!(dfs2.freqs_f[2:end], dfs2.n_f[2:end],
#     yscale=:log10,
#     linewidth=1.5,
#     label="ODE: t="*string(timeODE2))

# plot!(freqs_f[2:end], vafSimAv_n[2:end],
#     yscale=:log10,
#     linewidth=0.8,
#     label="sim average")
    
# plot!(liveFreqsAv_f, liveNvAv_f,
#     linewidth=0.8,
#     label="live average")

# plot!(dfsGrowth.freqs_f[2:end], dfsGrowth.n_f[2:end],
#     linewidth=1.5,
#     label="ODE: growth")


# xlims!(0,0.3)
# ylims!(0.5*10^0, 10^4)

# display(p5)

# ##

# testSim = 7

# p6 = bar(liveFreqsS_Sim_f[testSim], liveNvS_Sim_f[testSim],
#     yscale=:log10,
#     linealpha=0,
#     fillalpha=0.5,
#     label="single sim")

# # for sim in 1:5
# #     # bar!(liveFreqs_Sim_f[sim], liveNv_Sim_f[sim],
# #     #     # bar_width=1.,
# #     #     linealpha=0,
# #     #     fillalpha=0.5,
# #     #     label="")
# #     plot!(liveFreqs_Sim_f[sim], liveNv_Sim_f[sim],
# #         # bar_width=1.,
# #         linealpha=0.8,
# #         linewidth=0.5,
# #         label="")
# # end

# plot!(dfs2S.freqs_f[2:end], dfs2S.n_f[2:end],
#     yscale=:log10,
#     linewidth=1.5,
#     label="ODE: t="*string(timeODE2))

# # plot!(freqs_f[2:end], vafSimAv_n[2:end],
# #     yscale=:log10,
# #     linewidth=1.,
# #     label="sim average")

# # plot!(dfsGrowthS.freqs_f[2:end], dfsGrowthS.n_f[2:end],
# #     linewidth=1.5,
# #     label="ODE: growth")

# plot!(liveFreqsSAv_f, liveNvSAv_f,
#     # linewidth=0.8,
#     label="live average")
    
# xlims!(0,0.8)
# ylims!(0.5*10^0, 10^4)

# display(p6)













##

# simBars = Gadfly.layer(x=freqs_f[2:end], y=vaf_n[2:end], Geom.bar,
# Theme(default_color=color("grey50"))
# )
# simAvLine = Gadfly.layer(x=freqs_f[2:end],y = nVAv_f[2:end],
# # Geom.line,
# Geom.smooth(method=:loess,smoothing=0.002),
# Theme(default_color=color("steelblue1"), 
# bar_spacing=-2mm,
# line_width=1pt)
# )

# # # averaged sims
# # plot!(freqs_f[2:end], nVAv_f[2:end],
# #     color=1,
# #     linewidth=0.8,
# #     linestyle=:solid,
# #     linealpha=0.7,
# #     label="all variants \naveraged simulations")

# fig2b = Gadfly.plot(simAvLine, simBars, 
# Scale.y_log10,
# Coord.cartesian(xmin=0, xmax=1, ymin=0.1, ymax=4),
# Guide.xlabel("VAF"),
# Guide.ylabel("number of variants")
# )
# display(fig2b)
##

# set_default_plot_size(6inch, 4inch)
# fig2b |> PDF("figures/paper/gadflyTest.pdf")

##