using JLD2, Glob, DataFrames, FileIO
using Statistics
using StatsPlots
using DataFrames
using LaTeXStrings
##


_Nf = [1000, 2000, 4000, 7000, 10000]
nVaf1_Nf_Sim = Vector{Int}[]
nVaf1S_Nf_Sim = Vector{Int}[]
nVaf1Av_Nf = Float64[]
nVaf1Var_Nf = Float64[]
nVaf1Std_Nf = Float64[]
nVaf1SAv_Nf = Float64[]
nVaf1SVar_Nf = Float64[]
nVaf1SStd_Nf = Float64[]
vaf1Data = DataFrame(
    vaf1 = Int[],
    Nf = Int[]
)

vaf1SData = DataFrame(
    vaf1 = Int[],
    Nf = Int[]
)

for Nf in _Nf

    folder = "data/Simulations/Nf"*string(Nf)
    fnames_ = glob("singlePatientFullSim_Ni1_*.jld2", folder)

    nVar_Sim_F = Vector{Int}[]
    nVaf1_sim = Int[]
    nVaf1S_sim = Int[]
    for fname in fnames_
        paramsTrue, nVarSim_f, nVarSimS_f = load(fname, "paramsTrue", "nVarSim_f", "nVarSimS_f")
        # push!(simData, [nVarSim_f, paramsTrue["N final"]])
        push!(nVar_Sim_F, nVarSim_f)
        push!(nVaf1_sim, nVarSim_f[2])
        push!(nVaf1S_sim, nVarSimS_f[2])
    end

    push!(nVaf1_Nf_Sim, nVaf1_sim)
    push!(nVaf1S_Nf_Sim, nVaf1S_sim)

    nVaf1Av = mean(nVaf1_sim)
    nVaf1Var = var(nVaf1_sim)
    nVaf1Std = std(nVaf1_sim)

    nVaf1SAv = mean(nVaf1S_sim)
    nVaf1SVar = var(nVaf1S_sim)
    nVaf1SStd = std(nVaf1S_sim)

    push!(nVaf1Av_Nf, nVaf1Av)
    push!(nVaf1Var_Nf, nVaf1Var)
    push!(nVaf1Std_Nf, nVaf1Std)
    push!(nVaf1SAv_Nf, nVaf1SAv)
    push!(nVaf1SVar_Nf, nVaf1SVar)
    push!(nVaf1SStd_Nf, nVaf1SStd)

    vaf1DataNf = DataFrame(
        vaf1 = nVaf1_sim,
        Nf = Nf*ones(length(nVaf1_sim))
    )
    vaf1SDataNf = DataFrame(
        vaf1 = nVaf1S_sim,
        Nf = Nf*ones(length(nVaf1S_sim))
    )

    append!(vaf1Data, vaf1DataNf)
    append!(vaf1SData, vaf1SDataNf)
end

##
pyplot()
figscaler = 1
theme(:default, minorgrid=false, gridstyle=:dash, fontfamily="DejaVu Sans", size=(figscaler*600,figscaler*400))
# histogram(nVaf1_Nf_Sim[5], bins=20)

##

# fig1 = scatter(_Nf, nVaf1Av_Nf, legend=:topleft)
# xlabel!(L"population size $N$")
# ylabel!(L"number of variants at $f=1/N$")

##
fig2 = groupedboxplot(
    vaf1Data[:,:Nf], vaf1Data[:,:vaf1], 
    legend=:topleft,
    color=:grey50,
    label="Summary statistics",
    xticks=[0, 1000, 2000, 4000, 7000, 10000],
    linewidth=1,
    fillalpha=0.8,
    markersize=0,
    markeralpha=0,
)
groupeddotplot!(
    vaf1Data[:,:Nf], vaf1Data[:,:vaf1], 
    # yscale=:log10,
    label="Single simulations",
    xticks=[0, 1000, 2000, 4000, 7000, 10000],
    linewidth=0,
    fillalpha=0.7,
    # markersize=1.,
    marker=(:black,stroke(0),1.5),
)
xlabel!(L"Population size $N$")
ylabel!(L"Number of variants at $f=1/N$")

##
fig2 = groupedboxplot(
    vaf1SData[:,:Nf], vaf1SData[:,:vaf1], 
    legend=:topleft,
    color=:grey50,
    label="Summary statistics",
    xticks=[0, 1000, 2000, 4000, 7000, 10000],
    linewidth=1,
    fillalpha=0.8,
    markersize=0,
    markeralpha=0,
)
groupeddotplot!(
    vaf1SData[:,:Nf], vaf1SData[:,:vaf1], 
    # yscale=:log10,
    label="Single sampled simulations",
    xticks=[0, 1000, 2000, 4000, 7000, 10000],
    linewidth=0,
    fillalpha=0.7,
    # markersize=1.,
    marker=(:black,stroke(0),1.5),
)
xlabel!(L"Population size $N$")
ylabel!(L"Number of variants at $f=1/S$")




##
scatter(_Nf, nVaf1Std_Nf ./ nVaf1Av_Nf, legend=:topright)

##
scatter(_Nf, nVaf1SAv_Nf, legend=:topleft)

scatter(_Nf, nVaf1SStd_Nf, legend=:topleft)

##
figRelDispTotal = scatter(
    _Nf, nVaf1Std_Nf ./ nVaf1Av_Nf,
    legend=:topright,
    label="Full population",
    # color=:grey50,
    markersize=8,
    markershape=:diamond
    # size=(0.8*600,0.8*400),
)
scatter!(
    _Nf, nVaf1SStd_Nf ./ nVaf1SAv_Nf,
    legend=:topright,
    label="Sample",
    # color=:grey50,
    markersize=8,
    markershape=:utriangle,
    # size=(0.8*600,0.8*400),
)
xlabel!(L"Population size $N$")
ylabel!(L"Relative dispersion $\, \mathrm{Std}(v_{1}) / \mathrm{E}(v_{1})$")
xlims!(0, 11000)
ylims!(0,.2)

# figRelDisp = plot(figRelDispTotal, figRelDispSample, layout=(2,1), size=(figscaler*600,figscaler*500))