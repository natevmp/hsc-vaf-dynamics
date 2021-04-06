##
using Plots
pyplot()
# theme(:lime)
theme(:juno)
# theme(:sand)
using JLD2

##
include("../../src/vafdyn.jl")
using .VAFDyn
include("../../src/vafSim.jl")
using .VAFSim


##
params = Dict(
    "N initial" => 500,
    "N final" => 500,
    "μ" => 1.2,
    "λ" => 4.2,
    "p" => 0.5,
    "sample size" => 89,
    "mature time" => 15,
    "evolve time" => 59
)

growthRateFromNT(Nf, t) = log(Nf)/t

function extendParams!(params::Dict)
    params["ρ"] = params["λ"]*(1-params["p"])
    params["ϕ"] = params["λ"]*params["p"]
    params["N"] = params["N final"]
    γ = growthRateFromNT(params["N final"], params["mature time"])
    params["growth rate"] = γ
    return params
end
extendParams!(params)

## simulation
nSims = 20
nVarsSCAv_Sim_f = []
nVarsSCAvS_Sim_f = []
for sim in 1:nSims
    @time times_t, nLive_t, nV_f, nVS_f, burden_m, burdenS_m, _, __, nVarsSC_cid_f, nVarsSCS_cid_f = VAFSim.birthDeathFixedGrowth(params, params["evolve time"], 0.1, singleCellSpectrum=true, showprogress=true)
    nVarsSCAv_f = vec(sum(nVarsSC_cid_f, dims=1))./(size(nVarsSC_cid_f)[1])
    nVarsSCAvS_f = vec(sum(nVarsSCS_cid_f, dims=1))./(size(nVarsSCS_cid_f)[1])
    push!(nVarsSCAv_Sim_f, nVarsSCAv_f)
    push!(nVarsSCAvS_Sim_f, nVarsSCAvS_f)
end

nVarsSCAv_sim_f = Array{Float64, 2}(undef, nSims, params["N final"])
nVarsSCAvS_sim_f = Array{Float64, 2}(undef, nSims, params["sample size"])
for i in 1:nSims
    nVarsSCAv_sim_f[i, :] .= nVarsSCAv_Sim_f[i]
    nVarsSCAvS_sim_f[i, :] .= nVarsSCAvS_Sim_f[i]
end


## PDE evolve single cell VAF spectrum

vfsSC = VAFDyn.VFreqspace(params["N final"], 100)
# VAFDyn.evolveGrowingSCVAF(vfsSC, params, params["evolve time"])
@time VAFDyn.evolveSCVAF(vfsSC, params, params["evolve time"])
dfsSC = VAFDyn.makeDFSfromVFS(vfsSC, params["N final"])
dfsSCS = VAFDyn.sampler(dfsSC, 89)

## Plot PDE evolved SC VAF spectrum
# freqs_f = (0:params["N final"])./params["N final"]
# freqsS_f = (0:params["sample size"])./params["sample size"]
fig0 = plot(dfsSC.freqs_f[2:end],dfsSC.n_f[2:end])
display(fig0)

## Plot SC VAF spectrum
freqs_f = (0:params["N final"])./params["N final"]
freqsS_f = (0:params["sample size"])./params["sample size"]

for cid in 5:10
    plot = bar(freqsS_f[2:end], nVarsSCS_cid_f[cid, :])
    # ylims!(0, 200)
    xlabel!("VAF")
    ylabel!("# variants in single cell")
    display(plot)
end

##
# nVarsSCAvS_f = vec(sum(nVarsSCS_cid_f, dims=1))./(size(nVarsSCS_cid_f)[1])
# nVarsSCAv_f = vec(sum(nVarsSC_cid_f, dims=1))./(size(nVarsSC_cid_f)[1])
nVarsSCAv_f = vec(sum(nVarsSCAv_sim_f, dims=1) ./ nSims)
nVarsSCAvS_f = vec(sum(nVarsSCAvS_sim_f, dims=1) ./ nSims)


fig1a = bar(freqsS_f[2:end], nVarsSCAvS_f)
xlabel!("VAF")
ylabel!("average variants per cell (Sample)")
display(fig1a)

fig1b = bar(freqs_f[2:end], nVarsSCAv_f)
xlabel!("VAF")
ylabel!("average variants per cell (True pop)")
display(fig1b)

##

fig2a = bar(freqsS_f[2:end], nVS_f[2:end], yscale=:log10)
xlabel!("VAF")
ylabel!("# variants in sample")
display(fig2a)

fig2b = bar(freqs_f[2:end], nV_f[2:end], yscale=:log10)
xlabel!("VAF")
ylabel!("# variants in pop")
display(fig2b)