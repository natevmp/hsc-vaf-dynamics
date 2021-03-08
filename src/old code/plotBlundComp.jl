using Plots
using JLD2
using Distributions
using Statistics

# trec is all the times the simulation records
trec = 1:1:200
tmax = trec[end]
#r is the rate at which replication events happen per individual
r = 1
#s is the percentage of surviving individuals for
d = 1
d_all = 2^0
# n is the number of individuals
N = 100
#indStart = 100
# bsize is the size of the population at the bottleneck
Nbn = 10
# mu is the mutation rate
μ = 1
#μ_all = 1 #0.2:0.2:2
# runs is the number of simulations to average over
runs = 1
# p is the likelihood of self-replacement
p = 0
# C is the number of somatic cells the population is supposed to produce per time step
#C = 0
# lin is the linear growth rate of the population
lin = 0 #5*round(10^(-2),digits=4)

data_name = string("../data/bottleneck/dataBDpres", "_N", N, "_mu", μ, "_lin", lin,"_Nbn", Nbn, "_t", tmax, "_r", r,"_d", d, ".jld2")

@load data_name vafM_n_t vafBM_n_t indLiveM_t trec mFixedM_t burdenM_m_t burdenBM_m_t #distanceM_m distanceBM_m mLiveDist

f = 1:Nbn

t = 10

vafBlund_n = 2*N*μ ./ f .* exp.(- f ./ t)

vafBlundSampled_n = 2*N*μ ./ f .* exp.(- f .* N ./ (Nbn * t))

plot(f,vafBlund_n, label = "Blund")

plot!(f,vafBlundSampled_n, label = "Blund Sampled")

plot!(f,vafM_n_t[f .+ 1,10], label = "Full Pop")

plot!(f,vafM_n_t[f .+ 1,1], label = "Full Pop t = 1")

plot!(f,vafBM_n_t[f .+ 1,10], label ="Sampled Pop")

title!(string("Comparison Blundell simulation t = 10"))
xlabel!("Prevalences")
ylabel!("# of mutations")


fig_name = string("../Figures/BlundellComp_N",N, "_t", tmax,"_mu", μ)

savefig(fig_name)
