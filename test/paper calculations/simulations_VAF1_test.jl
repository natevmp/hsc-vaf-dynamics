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
# pyplot()
# figscaler = 1
# theme(:default, minorgrid=false, gridstyle=:dash, fontfamily="DejaVu Sans", size=(figscaler*600,figscaler*400))
pyplot()
theme(:default,
    # minorgrid=false,
    gridstyle=:dash,
    fontfamily="DejaVu Sans",
    showaxis=true,
    gridlinewidth=0.7,
    # size=(0.9*500,0.9*400),
    size=(500,400),
    legendfontsize=10,
    guidefontsize=12,
    tickfontsize=10,
)
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

##

fig1 = plot(
    range(0,1,length=lenF), nAv_f,
    yscale=:log10, xscale=:log10,
    label=L"V_f",
)
plot!(
    range(0,1,length=lenF), nAv_f .+ nStd_f,
    linestyle=:dash,
    color=2,
    alpha=0.6,
    label=L"V_f \pm \sigma_f",
)
plot!(
    range(0,1,length=lenF), nAv_f .- nStd_f,
    linestyle=:dash,
    color=2,
    alpha=0.6,
    label="",
)
xlims!(1/lenF,1E-1)
ylims!(0.1, maximum(nAv_f))
xlabel!(L"variant allele frequency  $f$")
ylabel!("number of variants")
display(fig1)

##

fig2 = plot(
    # range(0,1,length=lenF), nAv_f,
    yscale=:log10, xscale=:log10,
    # label=L"V_f",
)
plot!(
    range(0,1,length=lenF), nStd_f[nAv_f.!=0] ./ nAv_f[nAv_f.!=0],
    color=:black,
    alpha=0.8,
    label=L"\sigma_f / V_f",
)
xlims!(1/lenF,1E-1)
ylims!(1E-2, 1E2)
xlabel!(L"variant allele frequency  $f$")
# ylabel!(L"\sigma_{V_f}\, / V_f")
ylabel!(L"relative standard deviation  $\sigma_{V_f}\, / V_f$")

display(fig2)