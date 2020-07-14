using Plots
using JLD2

##Script to run the simulation of neutrally evolving populations

# tmax is the amount of time that the simulation takes
tmax = 20
#r is the rate at which replication events happen per individual
r = 1
# n is the number of individuals
N = 50
# bsize is the size of the population at the bottleneck
Nbn = 10
# mu is the mutation rate
μ = 1
# p is the likelihood of self-replacement
p = 0.0

mFixedM_all = zeros(19)

data_name = string("../data/bottleneck/dataBD", "_N", N, "_mu", μ, "_p", p,"_Nbn", Nbn, "_t", tmax, "_r", r, ".jld2")

@load data_name vafM_n vafBM_n mFixedM burdenM_m burdenBM_m

plot(1:(N*μ*2),burdenM_m[1:N*μ*2],label=string("t = ", tmax),legend=:topright,dpi=300)
k = 1

mFixedM_all[k] = mFixedM

for tmax = 40:20:80
    global k += 1
    data_name = string("../data/bottleneck/dataBD", "_N", N, "_mu", μ, "_p", p,"_Nbn", Nbn, "_t", tmax, "_r", r, ".jld2")

    @load data_name vafM_n vafBM_n mFixedM burdenM_m burdenBM_m
    mFixedM_all[k] = mFixedM
    plot!(1:(N*μ*5),burdenM_m[1:N*μ*5],label=string("t = ", tmax),legend=:topright,dpi=300)

end

fig_name = string("../Figures/burden_all_p", p, "_N", N)

savefig(fig_name)
