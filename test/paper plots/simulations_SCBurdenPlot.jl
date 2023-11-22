using Statistics
using Plots, LaTeXStrings
# pyplot()
using StatsBase
using JLD2
using HypothesisTests

include("../../src/compoundPoisson.jl")
using .CompoundPoisson

##

SAVEFIG = true


growthRateFromNT(Nf, t) = log(Nf)/t
function extendParams!(params::Dict)
    params["ρ"] = params["λ"]*(1-params["p"])
    params["ϕ"] = params["λ"]*params["p"]
    params["N"] = params["N final"]
    γ = growthRateFromNT(params["N final"], params["mature time"])
    params["growth rate"] = γ
    return params
end

@load "./data/SingleCellBurden/singleCellBurden_Nf10000_nSims1000.jld2"

evolveTime=paramsTrue["evolve time"]

## =================== Data calculations ========================

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


## ------- Compute compound poisson distribution ---------
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
λEff = getEffectiveDivisions(paramsTrue)
cpData_id = CompoundPoisson.randComPois(λEff, paramsTrue["μ"], 800000)


## ============================= Plots ============================
using Plots, LaTeXStrings, Colors, ColorSchemes
pyplot()
theme(:default,
    minorgrid=false, 
    gridstyle=:dash, 
    fontfamily="DejaVu Sans", 
    legendfontsize=9, 
    size=(500,400),
)

## ----------- Mutational Burden distribution ------------------

figSCBurden = plot(
    legendfontsize=8,
    title = "a)",
    titleloc = :left,
    titlefont=font(20, "DejaVu Sans"),
)
stephist!(cpData_id,
    color=:black,
    linealpha=0.6,
    normalize=true,
    bins=250,
    linewidth=1.6,
    linestyle=:dash,
    label=L"Comp. Poisson: $\mu=1.20$"
)
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
simTestId = 5
stephist!(scBurden_Sim_cid[simTestId], 
    color=1,
    # bins=50, 
    normalize=true, 
    linewidth=1.1, 
    linealpha=0.8,
    label=L"single sim: $\tilde{\mu}=$"*latexstring(round(mutRate_sim[simTestId], digits=2)),
)
stephist!(scBurdenSample_Sim_cid[simTestId], 
    color=2,
    bins=20, 
    normalize=true, 
    linewidth=1.1, 
    linealpha=0.8,
    label=L"sampled sim: $\tilde{\mu}=$"*latexstring(round(mutRateS_sim[simTestId], digits=2))
)
xlims!(330,630)
xlabel!(L"Single cell mutational burden $m$")
ylabel!("Density of cells")
display(figSCBurden)
SAVEFIG && savefig(figSCBurden, "Figures/Paper/3a.pdf")


## ===== Plot simulation mutation rate distribution =====

avMutTrue = round(mutRateSimMean, digits=2)
stdMutTrue = round(mutRateSimStd, digits=2)
avMutSample = round(mutRateSimSampleMean, digits=2)
stdMutSample = round(mutRateSimSampleStd, digits=2)

figMutRateDist = stephist(mutRate_sim,
    title = "b)",
    titleloc = :left,
    titlefont=font(20, "DejaVu Sans"),
    # bins=50,
    fillalpha=0.3,
    fillrange=0,
    linestyle=:dash,
    linewidth=1.5,
    normalize=:pdf,
    label="population \n"*
        L"\operatorname{E}(\tilde{\mu}) = %$avMutTrue"*"\n"*
        L"\sigma (\tilde{\mu}) = %$stdMutTrue",
)
stephist!(mutRateS_sim,
    # bins=50,
    fillalpha=0.3,
    fillrange=0,
    linewidth=1.5,
    normalize=:pdf,
    label="sample \n"*
        L"\operatorname{E}(\tilde{\mu}) = %$avMutSample"*"\n"*
        L"\sigma (\tilde{\mu}) = %$stdMutSample"
)
    
xlabel!(L"Inferred mutation rate $\tilde{\mu}$")
ylabel!("Density of simulations")
# title!("μ true = $μTrue, N mature = $NMature")
display(figMutRateDist)
SAVEFIG ? savefig(figMutRateDist, "Figures/Paper/3b.pdf") : 0

##