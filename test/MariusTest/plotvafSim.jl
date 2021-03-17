using Plots
using JLD2

##Script to run the simulation of neutrally evolving populations

# tmax is the amount of time that the simulation takes
tmax = 10
#r is the rate at which replication events happen per individual
r = 1
# n is the number of individuals
N = 200
# bsize is the size of the population at the bottleneck
Nbn = 10
# mu is the mutation rate
μ = 1
# p is the likelihood of self-replacement
p = 0.9

data_name = string("../data/bottleneck/dataBD", "_N", N, "_mu", μ, "_p", p,"_Nbn", Nbn, "_t", tmax, "_r", r, ".jld2")

@load data_name vafM_n vafBM_n mFixedM

mAll1 = 0

for k = 1:Nbn
    global mAll1 += k*vafBM_n[k+1]
end

mAll1 += mFixedM*Nbn

rEst1 = mAll1*(1+p)/2/Nbn/tmax/μ

rEstmin = mAll1/2/Nbn/tmax/μ

rEstmax = mAll1/Nbn/tmax/μ

mAll2 = 0

for k = 1:N
    global mAll2 += k*vafM_n[k+1]
end

mAll2 += mFixedM*N

rEst2 = mAll2*(1+p)/2/N/tmax/μ

println(rEstmin)

println(rEstmax)
