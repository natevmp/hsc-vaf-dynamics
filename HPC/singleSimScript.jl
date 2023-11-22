
include("../src/vafSim.jl")
using .VAFSim

using JLD2

simNumber = parse(Int, ARGS[1])
N = parse(Int, ARGS[2])

<<<<<<< HEAD
=======
# simNumber = 1
# N = 150

>>>>>>> nate
## ===== User params =====
paramsTrue = Dict{String,Real}(
    "N initial" => 1,
    "N final" => N,
<<<<<<< HEAD
    # "N final" => 100,
    "μ" => 2.0,
    "λ" => 5.0,
    "p" => 0.4,
    "sample size" => 89,
    "mature time" => 15,
    "evolve time" => 60
=======
    "μ" => 2.0,
    "λ" => 5.0,
    "p" => 0.4,
    "sample size" => 100,
    "mature time" => 5,
    "evolve time" => 60,
    "pure births" => 100,
>>>>>>> nate
)
tSave = 0.1

## Parameter functions

growthRateFromNT(Nf, t) = log(Nf)/t
function extendParams!(params::Dict)
    params["ρ"] = params["λ"]*(1-params["p"])
    params["ϕ"] = params["λ"]*params["p"]
    params["N"] = params["N final"]
    γ = growthRateFromNT(params["N final"], params["mature time"])
    params["growth rate"] = γ
    return params
end
extendParams!(paramsTrue)

println(paramsTrue)

## Single Patient simulation
println("Running simulation...")

timesSim_, nCellSim_t, nVarSim_f, nVarSimS_f, nCellSim_m, nCellSimS_m, mSim_cid, mSimS_cid = VAFSim.birthDeathFixedGrowth(paramsTrue, paramsTrue["evolve time"], tSave, showprogress=true)

println("Simulation completed. Saving data...")
<<<<<<< HEAD
filename = "singlePatientFullSim_Nf"*string(paramsTrue["N final"])*"_sim"*string(simNumber, pad=2)*".jld2"
=======
filename = "singlePatientFullSim"*
    "_Ni"*string(paramsTrue["N initial"])*
    "_Nf"*string(paramsTrue["N final"])*
    "_tM"*string(paramsTrue["mature time"])*
    "_NH"*string(paramsTrue["pure births"])*
    "_sim"*string(simNumber, pad=2)*
    ".jld2"

>>>>>>> nate
@save filename paramsTrue timesSim_ nCellSim_t nVarSim_f nVarSimS_f nCellSim_m nCellSimS_m mSim_cid mSimS_cid

