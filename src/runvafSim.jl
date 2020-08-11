include("vafSim.jl")
using JLD2

##Script to run the simulation of neutrally evolving populations

# tmax is the amount of time that the simulation takes
tmax = 1.7
#r is the rate at which replication events happen per individual
r = 1
# n is the number of individuals
Nmax = 1400
# bsize is the size of the population at the bottleneck
Nbn = 140
# mu is the mutation rate
μ = 1
# runs is the number of simulations to average over
runs = 100
# p is the likelihood of self-replacement
p = 0.0


#loop to run the simulations
for N = 140:140:1400
    vafBM_n= zeros(Nbn+1)

    vafM_n = zeros(N+1)

    burdenBM_m= zeros(Int(round(tmax*r*μ*N*4)))

    burdenM_m = zeros(Int(round(tmax*r*μ*N*4)))

    mLiveDist = zeros(Int(round(N*μ*100)))

    distanceBM_m= zeros(Int(round(N*μ*40)))

    distanceM_m = zeros(Int(round(N*μ*40)))

    mFixedM = 0

    pc = ((Nmax/N)-1.0)/((Nmax/N)-(1/2))
    rc = r*1.0/(1.0-(pc/2.0))

    #tc = N/20*tmax
    #μc = 20/N*μ

    for run = 1:runs

        vaf_n,vafB_n,mFixed,mLive,burden_m,burdenB_m = VAFSim.birthDeathAlt(N, μ, pc, Nbn, tmax, rc)

        #burdenB_m .+= mFixedM
        #burden_m .+= mFixedM

        #distanceBM_m[1:mLive+1] .+= distanceB_m
        #distanceM_m[1:mLive+1] .+= distance_m

        burdenBM_m[(1+mFixed):(mLive+mFixed)] .+= burdenB_m
        burdenM_m[(1+mFixed):(mLive+mFixed)] .+= burden_m

        mLiveDist[mLive] += 1

        vafBM_n += vafB_n
        vafM_n += vaf_n

        mFixedM  += mFixed/runs

    end

    #average the VAFs and save them

    #distanceM_m = distanceM_m ./ runs
    #distanceBM_m = distanceBM_m ./ runs

    burdenM_m = burdenM_m ./ runs
    burdenBM_m = burdenBM_m ./ runs

    vafM_n = vafM_n ./ runs
    vafBM_n = vafBM_n ./ runs

    data_name = string("../data/bottleneck/dataBDtestabc", "_N", N, "_mu", μ, "_p", p,"_Nbn", Nbn, "_t", tmax, "_r", r, ".jld2")

    @save data_name vafM_n vafBM_n mFixedM burdenM_m burdenBM_m distanceM_m distanceBM_m mLiveDist

end
