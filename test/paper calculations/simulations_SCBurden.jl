using ProgressMeter, Statistics
using Plots
pyplot()
using StatsBase
using LaTeXStrings
using JLD2
using HypothesisTests

##
include("../../src/vafSim.jl")
using .VAFSim
include("../../src/burdenDyn.jl")
using .BurdenDyn
include("../../src/burdenSim.jl")
using .BurdenSim
include("../../src/compoundPoisson.jl")
using .CompoundPoisson

##

LOADDATA = true
SAVEFIG = false


growthRateFromNT(Nf, t) = log(Nf)/t
function extendParams!(params::Dict)
    params["ρ"] = params["λ"]*(1-params["p"])
    params["ϕ"] = params["λ"]*params["p"]
    params["N"] = params["N final"]
    γ = growthRateFromNT(params["N final"], params["mature time"])
    params["growth rate"] = γ
    return params
end


##

if LOADDATA
    @load "./data/SingleCellBurden/singleCellBurden_Nf10000_nSims1000.jld2"
else
    # ===== run limited stochastic simulation =====
    
    nSims = 500
    paramsTrue = Dict{String,Real}(
        "N initial" => 500,
        "N final" => 500,
        "μ" => 1.2,
        "λ" => 4.2,
        "p" => 0.5,
        "sample size" => 89,
        "mature time" => 15,
        "evolve time" => 80
    )
    extendParams!(paramsTrue)

    evolveTime = paramsTrue["evolve time"]
    tSaveStep = 0.5
    display(paramsTrue)
    println("runtime: ", evolveTime)
    println("timestep: ", tSaveStep)

    nCells_Sim_m = Vector{Int64}[]
    nSize_Sim_t = Vector{Int64}[]
    nCellsSample_Sim_m = Vector{Int64}[]
    scBurden_Sim_cid = Vector{Int64}[]
    scBurdenSample_Sim_cid = Vector{Int64}[]

    nCellsAvSim_m = zeros(Float64, Int(round(3*evolveTime*2paramsTrue["μ"]*paramsTrue["λ"])))
    scBurden_cid, times_t, nSize_t = BurdenSim.evolveSCBurden(paramsTrue, evolveTime, tSaveStep)
    nCellsSCB_m = BurdenSim.burdenHist(scBurden_cid)
    scBurdenSample_cid = BurdenSim.sampler(scBurden_cid, paramsTrue["sample size"])
    nCellsSample_m = BurdenSim.burdenHist(scBurdenSample_cid)
    nSizeAv_t = nSize_t/nSims
    nCellsAvSim_m[1:length(nCellsSCB_m)] += nCellsSCB_m[1:end]/nSims
    push!(nCells_Sim_m, nCellsSCB_m)
    push!(nCellsSample_Sim_m, nCellsSample_m)
    push!(nSize_Sim_t, nSize_t)
    push!(scBurden_Sim_cid, scBurden_cid)
    push!(scBurdenSample_Sim_cid, scBurdenSample_cid)

    @showprogress for sim in 2:nSims
        scBurden_cid, times_t, nSize_t = BurdenSim.evolveSCBurden(paramsTrue, evolveTime, tSaveStep)
        nCellsSCB_m = BurdenSim.burdenHist(scBurden_cid)
        scBurdenSample_cid = BurdenSim.sampler(scBurden_cid, paramsTrue["sample size"])
        nCellsSample_m = BurdenSim.burdenHist(scBurdenSample_cid)
        global nSizeAv_t += nSize_t/nSims
        nCellsAvSim_m[1:length(nCellsSCB_m)] += nCellsSCB_m[1:end]/nSims
        push!(nCells_Sim_m, nCellsSCB_m)
        push!(nCellsSample_Sim_m, nCellsSample_m)
        push!(nSize_Sim_t, nSize_t)
        push!(scBurden_Sim_cid, scBurden_cid)
        push!(scBurdenSample_Sim_cid, scBurdenSample_cid)
    end

    savename = "singleCellBurden_Nf"*string(paramsTrue["N final"])*"_nSims"*string(nSims)*".jld2"
    @save "data/SingleCellBurden/"*savename times_t nCellsSCB_m nCellsSample_m nSizeAv_t nSize_t nCellsAvSim_m nCells_Sim_m nCellsSample_Sim_m nSize_Sim_t scBurden_Sim_cid scBurdenSample_Sim_cid paramsTrue nSims

end
evolveTime=paramsTrue["evolve time"]


## ===== run ODE expected value evolve ======
# mMax = Int(round(2*paramsTrue["λ"]*evolveTime*paramsTrue["μ"]*3))
# ϵ = 0.0001
# # @time nCellsODE_m = BurdenDyn.evolveBurdenGrowth(paramsTrue, evolveTime; growthType="exponentialToLinear")
# @time nCellsODE_m = BurdenDyn.evolveBurdenGrowth(paramsTrue, evolveTime, mMax, ϵ)

## ===== Data calculations =====

# ode mutation rate
# burdenMeanODE = sum([m*nCellsODE_m[m+1] for m in 0:length(nCellsODE_m)-1])/sum(nCellsODE_m)
# burdenVarODE = sum([m^2*nCellsODE_m[m+1] for m in 0:length(nCellsODE_m)-1])/sum(nCellsODE_m) - burdenMeanODE^2
# mutRateODE = burdenVarODE/burdenMeanODE - 1

# average sim mutation rate
burdenMeanAvSim = sum([m*nCellsAvSim_m[m+1] for m in 0:length(nCellsAvSim_m)-1])/sum(nCellsAvSim_m)
burdenVarAvSim = sum([m^2*nCellsAvSim_m[m+1] for m in 0:length(nCellsAvSim_m)-1])/sum(nCellsAvSim_m) - burdenMeanAvSim^2
mutRateAvSim = burdenVarAvSim/burdenMeanAvSim - 1

# single sim mutation rate
burdenVar_sim = Vector{Float64}(undef, nSims)
burdenMean_sim = Vector{Float64}(undef, nSims)
mutRate_sim = Vector{Float64}(undef, nSims)
burdenSVar_sim = Vector{Float64}(undef, nSims)
burdenSMean_sim = Vector{Float64}(undef, nSims)
mutRateS_sim = Vector{Float64}(undef, nSims)
for sim in 1:nSims
    nCells_m = nCells_Sim_m[sim]
    nCellsS_m = nCellsSample_Sim_m[sim]
    burdenMean_sim[sim] = mean(scBurden_Sim_cid[sim])
    burdenVar_sim[sim] = var(scBurden_Sim_cid[sim])
    burdenSMean_sim[sim] = mean(scBurdenSample_Sim_cid[sim])
    burdenSVar_sim[sim] = var(scBurdenSample_Sim_cid[sim])
    mutRate_sim[sim] = burdenVar_sim[sim] / burdenMean_sim[sim] - 1
    mutRateS_sim[sim] = burdenSVar_sim[sim] / burdenSMean_sim[sim] - 1
end

mutRateSimMean = mean(mutRate_sim)
mutRateSimStd = std(mutRate_sim)
mutRateSimSampleMean = mean(mutRateS_sim)
mutRateSimSampleStd = std(mutRateS_sim)


# println("ODE estimated mutation rate: ", mutRateODE)
println("average sim estimated mutation rate: ", mutRateAvSim)
println("single sims estimated mutation rate -- mean: " , mutRateSimMean, "; std: ",  mutRateSimStd)
println("single sims sample estimated mutation rate -- mean: " , mutRateSimSampleMean, "; std: ",  mutRateSimSampleStd)



## ========== Plot pop size over time ==========
# nPlots=20

# # p1 = plot( times_t, (t->BurdenDyn.logisticGrowth(params["N initial"], params["N final"], params["growth rate"], t)).(times_t),
# plot1 = plot( times_t,
#     (t->BurdenDyn.exponentialCappedGrowth(paramsTrue["N initial"],paramsTrue["N final"], paramsTrue["growth rate"], t)).(times_t),
#     # (t->BurdenSim.exponentialToLinearGrowth(
#     #     paramsTrue["N initial"], 
#     #     paramsTrue["N final"], 
#     #     paramsTrue["exp growth rate"], 
#     #     paramsTrue["lin growth rate"], 
#     #     t
#     #     )).(times_t),
#     label="ode",
#     legend=:bottomright,
#     linewidth=1.5,
#     linestyle=:dash)
# plot!(times_t, nSizeAv_t,
#     label="sim",
#     linewidth=1.5,
#     linestyle=:dashdot
#     )
# plot!([paramsTrue["mature time"],], seriestype=:vline, label="", color="grey", linewidth=2)
# annotate!(paramsTrue["mature time"], 0, text(L"t_f", :bottom, :left))
# # xticks!([paramsTrue["mature time"],],[L"t_f",])
# # xticks!(0:10:paramsTrue["evolve time"])

# xlabel!("time")
# ylabel!("population size")
# # for i in 1:nPlots
# #     plot!(times_t, nSizeSimsAr_Sim_t[i], label="", legend=:topleft)
# # end
# display(plot1)




## ========== Plot burden histogram ==========
NMature = paramsTrue["N final"]
μTrue = paramsTrue["μ"]
# # nPlots = 50
# plot2 = plot(0:length(nCellsODE_m)-1, nCellsODE_m,
#     label="ode",
#     linewidth=2,
#     linestyle=:dash
#     )
# plot!(0:length(nCellsAvSim_m)-1, nCellsAvSim_m,
#     label="sim",
#     linewidth=2,
#     linestyle=:dashdot,
#     legend=:topleft
#     )
# annotate!(8, maximum(nCellsAvSim_m)/2, text( "var/mean - 1:\n ODE: "*string(round(mutRateODE,digits=2))*"\n Sim: "*string(round(mutRateAvSim,digits=2)), 9, :left, :top ))
# xlabel!("# mutations")
# ylabel!("# cells")
# title!("μ true = $μTrue, N mature = $NMature")
# # for sim in 10:10+nPlots
# #     plot!(0:length(nCellsSimsAr_Sim_m[sim])-1, nCellsSimsAr_Sim_m[sim],label="")
# # end
# xlims!(0,mMax/2)
# display(plot2)

# SAVEFIG ? savefig(plot2, "figures/burdenODESimCompare/scBurden_detExpCap_growthTime_"*"_nF"*string(paramsTrue["N final"])*".pdf") : nothing



## ===== Plot simulation mutation rate distribution =====
avMutTrue = round(mutRateSimMean, digits=2)
stdMutTrue = round(mutRateSimStd, digits=2)
avMutSample = round(mutRateSimSampleMean, digits=2)
stdMutSample = round(mutRateSimSampleStd, digits=2)
plot3 = stephist(mutRate_sim,
    # bins=50,
    fillalpha=0.3,
    fillrange=0,
    linestyle=:dash,
    linewidth=1.5,
    normalize=:pdf,
    label="total population \n"*
        L"\operatorname{E}(\tilde{\mu}) = %$avMutTrue"*"\n"*
        L"\sigma (\tilde{\mu}) = %$stdMutTrue")
    # *latexstring(round(mutRateSimMean, digits=2))*L"\nstd: "*latexstring(round(mutRateSimStd, digits=2)))

stephist!(mutRateS_sim,
    # bins=50,
    fillalpha=0.3,
    fillrange=0,
    linewidth=1.5,
    normalize=:pdf,
    label="sampled population \n"*
        L"\operatorname{E}(\tilde{\mu}) = %$avMutSample"*"\n"*
        L"\sigma (\tilde{\mu}) = %$stdMutSample")
    
xlabel!(L"\tilde{\mu}")
ylabel!("Density of simulations")
# title!("μ true = $μTrue, N mature = $NMature")
display(plot3)


## ===== Figure for paper =====
# pyplot()
# gr()
function getEffectiveDivisions(params)
    t = params["evolve time"]
    tM = params["mature time"]
    γ = params["growth rate"]
    ρ = params["ρ"]
    ϕ = params["ϕ"]
    Nf = params["N final"]
    Ni = params["N initial"]
    # λEff = ϕ*t + 2γ*tM + 2ρ*t - ρ*(1/Ni - 1/Nf)/γ - ρ*(t-tM)/Nf
    # λEff = (2ρ+ϕ)*t + (exp(-γ*tM)-1)*ρ/γ + 2*log(1+exp(γ*tM)) - log(4) + 2Nf*(t-tM)*γ/(Nf+1) - (t-tM)*ρ/Nf    # this only goes for Ni=1!
    λEff = (2ρ+ϕ)*t + (exp(-γ*tM)-1)*ρ/γ + 2*log(1+exp(γ*tM)) - log(4) - (t-tM)*ρ/Nf  
    return λEff
end
# getEffectiveDivisions(paramsTrue)*paramsTrue["μ"]
λEff = getEffectiveDivisions(paramsTrue)
##
cpData_id = CompoundPoisson.randComPois(λEff, paramsTrue["μ"], 800000)

##
plot4 = stephist(cpData_id,
    color=:black,
    linealpha=0.6,
    normalize=true,
    bins=250,
    linewidth=1.6,
    linestyle=:dash,
    label=L"Comp. Poisson: $\mu=1.20$")

# plot4 = plot(0:length(nCellsAvSim_m)-1, nCellsAvSim_m/sum(nCellsAvSim_m),
plot!(0:length(nCellsAvSim_m)-1, nCellsAvSim_m/sum(nCellsAvSim_m),
    color=:grey,
    linealpha=1,
    linewidth=1,
    fillalpha=0.5,
    # linestyle=:dashdot,
    # legend=:topleft,
    fillrange=0,
    label=L"sims average: $\tilde{\mu}=$"*latexstring(round(mutRateAvSim,digits=2))
    )

# plot!(0:length(nCellsODE_m)-1, nCellsODE_m / sum(nCellsODE_m),
# legend=:topright,
# color=1,
# linealpha=1,
# linewidth=3,
# grid=true,
# gridalpha=0.2,
# # size=(1200,400),
# dpi=600,
# linestyle=:dash,
# label=L"ODE: $\tilde{\mu}=$"*latexstring(round(mutRateODE, digits=2))
# )

sims = [5,]
markershapes = [:dtriangle, :diamond, :cross]
for (i, sim) in enumerate(sims)
    stephist!(scBurden_Sim_cid[sim], 
    # plot!(nCellsSample_Sim_m[i]/sum(nCellsSample_Sim_m[i]), 
    color=1,
    # bins=50, 
    normalize=true, 
    linewidth=1.1, 
    linealpha=0.8,
    # linestyle=:dash,
    # markershape=markershapes[i],
    # markershape=:dtriangle,
    # markeralpha=0.6,
    # markerstrokewidth=:0.4,
    # markerstrokealpha=:0.8,
    # markerstrokecolor=2i,
    # markersize=5,
    label=L"single sim: $\tilde{\mu}=$"*latexstring(round(mutRate_sim[sim], digits=2))
    )
    
    stephist!(scBurdenSample_Sim_cid[sim], 
    # plot!(nCellsSample_Sim_m[i]/sum(nCellsSample_Sim_m[i]), 
    color=2,
    bins=20, 
    normalize=true, 
    linewidth=1.1, 
    linealpha=0.8,
    # linestyle=:dash,
    # # markershape=markershapes[i],
    # markershape=:diamond,
    # markeralpha=0.6,
    # markerstrokewidth=:0.4,
    # markerstrokealpha=:0.8,
    # markerstrokecolor=2i+1,
    # markersize=5,
    label=L"sampled sim: $\tilde{\mu}=$"*latexstring(round(mutRateS_sim[sim], digits=2))
    )
end

xlims!(330,630)
xlabel!(L"Single cell mutational burden $m$")
ylabel!("Density of cells")
# title!("μ true = $μTrue, N mature = $NMature")
display(plot4)


##
# pyplot()
# pF = plot(p5, h3, p1, layout=3, size=(1000,800))
pF = plot(plot4, plot3, layout=2, size=(1000,400))
display(pF)

# savefig(pF, "figures/burdenODESimCompare/scBurden_detExpLin_"*string(nSims)*"sims_nI"*string(params["N initial"])*"_nF"*string(params["N final"])*".pdf")

SAVEFIG ? savefig(pF, "figures/burdenODESimCompare/scBurden_detExpCap_growthTime_"*"_nF"*string(paramsTrue["N final"])*".pdf") : nothing



##
scBurdenMultSample_cid = vcat(scBurdenSample_Sim_cid...)
scBurdenMult_cid = vcat(scBurden_Sim_cid...)
##

KSampleADTest(scBurdenMultSample_cid, cpData_id)
