using Plots
using JLD2
using Distributions

include("compoundPoisson.jl")

##Script to run the simulation of neutrally evolving populations

# tmax is the amount of time that the simulation takes
tmax = 100
#r is the rate at which replication events happen per individual
r = 1
# n is the number of individuals
N = 100
# bsize is the size of the population at the bottleneck
Nbn = 100
# mu is the mutation rate
μ = 1
# p is the likelihood of self-replacement
p = 0.0

lin = 0

d = 1

data_name = string("../data/bottleneck/dataBDlog", "_N", N, "_mu", μ, "_lin", lin,"_Nbn", Nbn, "_t", tmax, "_r", r,"_d", d, ".jld2")

@load data_name vaf_n_t_run vafB_n_t_run indLive_t_run trec mFixed_t_run burden_m_t_run burdenB_m_t_run #distanceM_m distanceBM_m mLiveDist
