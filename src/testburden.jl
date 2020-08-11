using Plots
using JLD2
using Distributions

include("compoundPoisson.jl")

##Script to run the simulation of neutrally evolving populations

# tmax is the amount of time that the simulation takes
tmax = 10
#r is the rate at which replication events happen per individual
r = 1
# n is the number of individuals
#N = 600
# bsize is the size of the population at the bottleneck
Nbn = 100
# mu is the mutation rate
μ = 1
# p is the likelihood of self-replacement
p = 0.0

burdenaverage = zeros(10)
burdenvar = zeros(10)

N_all = 100:100:1000

for i = 1:10

    N = N_all[i]
    data_name = string("../data/bottleneck/dataBDtestABC", "_N", N, "_mu", μ, "_p", p,"_Nbn", Nbn, "_t", tmax, "_r", r, ".jld2")

    @load data_name vafM_n vafBM_n mFixedM burdenM_m burdenBM_m distanceM_m distanceBM_m mLiveDist

    #distanceaverage = 0

    for k = 1:Int(round((tmax*r*μ*100)))
        #global distanceaverage += distanceM_m[k]*k
        global burdenaverage[i] += burdenBM_m[k]*k
    end
    #distanceaverage = distanceaverage/(N*(N-1)/2)
    burdenaverage[i] = burdenaverage[i]/Nbn



    for k = 1:Int(round((tmax*r*μ*100)))
        global burdenvar[i] += (burdenBM_m[k]*(k)^2)
    end

    burdenvar[i] = burdenvar[i]/(Nbn) .- burdenaverage[i]^2

end

r_est = (burdenaverage.^2)./(burdenvar.-burdenaverage)./tmax
mu_est = (burdenvar.-burdenaverage)./burdenaverage

println(r_est)
println(mu_est)
