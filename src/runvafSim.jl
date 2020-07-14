include("vafSim.jl")
using JLD2

##Script to run the simulation of neutrally evolving populations

# tmax is the amount of time that the simulation takes
#tmax = 300
#r is the rate at which replication events happen per individual
r = 1
# n is the number of individuals
N = 50
# bsize is the size of the population at the bottleneck
Nbn = 10
# mu is the mutation rate
μ = 1
# runs is the number of simulations to average over
runs = 1000
# p is the likelihood of self-replacement
p = 0.0


#loop to run the simulations
for tmax = 20:20:400
    vafBM_n= zeros(Nbn+1)

    vafM_n = zeros(N+1)

    burdenBM_m= zeros(Int(round(N*μ*40*(1-p/2)/(1-p))))

    burdenM_m = zeros(Int(round(N*μ*40*(1-p/2)/(1-p))))

    mFixedM = 0


    for run = 1:runs

        vaf_n,vafB_n,mFixed,mLive,burden_m,burdenB_m = VAFSim.birthDeathAlt(N, μ*(1-p)/(1-p/2), p, Nbn, tmax, r/(1-p))

        burdenBM_m[1:mLive] .+= burdenB_m
        burdenM_m[1:mLive] .+= burden_m

        vafBM_n += vafB_n
        vafM_n += vaf_n

        mFixedM  += mFixed

    end

    #average the VAFs and save them

    burdenM_m = burdenM_m ./ runs
    burdenBM_m = burdenBM_m ./ runs

    vafM_n = vafM_n ./ runs
    vafBM_n = vafBM_n ./ runs
    mFixedM = mFixedM ./ runs

    data_name = string("../data/bottleneck/dataBD", "_N", N, "_mu", μ, "_p", p,"_Nbn", Nbn, "_t", tmax, "_r", r, ".jld2")

    @save data_name vafM_n vafBM_n mFixedM burdenM_m burdenBM_m

end
