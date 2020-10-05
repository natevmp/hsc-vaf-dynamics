using Plots
using JLD2
using Distributions
using Statistics

include("compoundPoisson.jl")

##Script to run the simulation of neutrally evolving populations

# tmax is the amount of time that the simulation takes
trec = 1:1:10
tmax = 10
#r is the rate at which replication events happen per individual
r = 1
#do args
# d is a modifier to reduce th
#d_all = 2 .^ (1:5)
d = 2^0
# n is the number of individuals
N = 200
# bsize is the size of the population at the bottleneck
Nbn = 10
# mu is the mutation rate
#μ = 2.0
μ = 1
# p is the likelihood of self-replacement
p = 0
#C is the number of somatic cells the population produces per time step
C = 0
# lin is the linear growth rate
lin = 0

data_name = string("../data/bottleneck/dataBDlog", "_N", N, "_mu", μ, "_lin", lin,"_Nbn", Nbn, "_t", tmax, "_r", r,"_d", d, ".jld2")

@load data_name vafM_n_t vafBM_n_t indLiveM_t trec mFixedM_t burdenM_m_t burdenBM_m_t depthmaxM_t depthmaxBM_t depthmeanM_t depthmeanBM_t

plot(trec, depthmaxM_t,label="max depth",legend=:topleft)
plot!(trec, depthmeanM_t,label="mean depth")

plot!(trec, depthmaxBM_t,label="bottleneck max depth")
plot!(trec, depthmeanBM_t,label="bottleneck mean depth")



title!(string("Depth over time"))
xlabel!("time")
ylabel!("depth")


fig_name = string("../Figures/LogDepthB_N", N, "_t", tmax,"_p",p)

savefig(fig_name)
