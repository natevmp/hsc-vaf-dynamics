include("vafSim.jl")
using JLD2

##Script to run the simulation of neutrally evolving populations

params = Dict("N initial" => 1, "N final" => 10000, "μ" => 1, "p" => 0, "λ" => 1,"growth rate" => 1,"sample size" => 100)
tStop = 20
tSaveStep = 1



times_t, nLive_t, vaf_n_t, vafB_n_t, burden_m, burdenB_m,events = VAFSim.birthDeathGrowthLog(params, tStop, tSaveStep)



data_name = string("../data/bottleneck/dataBIGlog", "_N", params["N final"], "_mu", params["μ"],"_S", params["sample size"], "_t", tStop, "_gR", params["growth rate"], ".jld2")

@save data_name times_t nLive_t vaf_n_t vafB_n_t burden_m burdenB_m
