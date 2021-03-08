include("vafSim.jl")
using JLD2

##Script to run the simulation of neutrally evolving populations

# trec is all the times the simulation records
trec = 1:1:200
tmax = trec[end]
#r is the rate at which replication events happen per individual
r = 1
# N is the number of individuals
#N = 5
N_all = 10
# bsize is the size of the population at the bottleneck
Nbn = 5
# mu is the mutation rate
μ = 1
# runs is the number of simulations to average over
runs = 100
# p is the likelihood of self-replacement
p = 0.95
#p_all = 0.1:0.1:0.8

#trec = 1:1:10

#loop to run the simulations
for N in N_all
    vafBM_n_run = zeros(Int64,Nbn+1,runs)

    vafM_n_run = zeros(Int64,N+1,runs)

    #burdenBM_m_run = zeros(Int(round(tmax*r*μ*N*4)),runs)

    #burdenM_m_run = zeros(Int(round(tmax*r*μ*N*4)),runs)

    #mLiveDist = zeros(Int(round(N*μ*100)))

    distanceBM_m = zeros(Int(round(N*μ*40)),runs)

    distanceM_m = zeros(Int(round(N*μ*40)),runs)

    mFixedM = 0

    indLiveM_t = zeros(size(trec)[1])

    #pc = ((Nmax/N)-1.0)/((Nmax/N)-(1/2))
    #rc = r*1.0/(1.0-(pc/2.0))

    #tc = N/20*tmax
    #μc = 20/N*μ

    for run = 1:runs

        vaf_n, vafB_n, mFixed, mLive, distanceB_m,distance_m = VAFSim.birthDeathAlt(N, μ, p, Nbn, tmax, r)

        distanceBM_m[1:mLive+1] .+= distanceB_m
        distanceM_m[1:mLive+1] .+= distance_m

        #burdenBM_m_run[(1+mFixed):(mLive+1+mFixed),run] .+= burdenB_m
        #burdenM_m_run[(1+mFixed):(mLive+1+mFixed),run] .+= burden_m

        #mLiveDist[mLive] += 1

        vafBM_n_run[:,run] += vafB_n
        vafM_n_run[:,run] += vaf_n

        mFixedM  += mFixed/runs
        #indLiveM_t  .+= indLive_t/runs
    end

    #average the VAFs and save them

    distanceM_m = distanceM_m ./ runs
    distanceBM_m = distanceBM_m ./ runs

    #burdenM_m = burdenM_m ./ runs
    #burdenBM_m = burdenBM_m ./ runs

    #vafM_n = vafM_n ./ runs
    #vafBM_n = vafBM_n ./ runs

    data_name = string("../data/bottleneck/dataBDdist", "_N", N, "_mu", μ, "_p", p, "_Nbn", Nbn, "_t", tmax, "_r", r, ".jld2")

    @save data_name vafM_n_run vafBM_n_run indLiveM_t trec mFixedM distanceM_m distanceBM_m

end
