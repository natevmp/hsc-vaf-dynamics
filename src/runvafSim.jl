include("vafSim.jl")
using JLD2

##Script to run the simulation of neutrally evolving populations

# trec is all the times the simulation records
trec = 1:1:10
tmax = trec[end]
#r is the rate at which replication events happen per individual
r = 1
#s is the percentage of surviving individuals for
d = 1
d_all = 2^0
# n is the number of individuals
N = 200
indStart = 10
# bsize is the size of the population at the bottleneck
Nbn = 10
# mu is the mutation rate
μ = 1
#μ_all = 1 #0.2:0.2:2
# runs is the number of simulations to average over
runs = 100
# p is the likelihood of self-replacement
p = 0
# C is the number of somatic cells the population is supposed to produce per time step
#C = 0
# lin is the linear growth rate of the population
lin = 0 #5*round(10^(-2),digits=4)

#trec = 1:1:10

#loop to run the simulations
#for μ in μ_all
    vafBM_n_t = zeros(Nbn+1,size(trec)[1])

    vafM_n_t = zeros(N+1,size(trec)[1])

    burdenBM_m_t = zeros(Int(round(tmax*r*μ*N)),size(trec)[1])

    burdenM_m_t = zeros(Int(round(tmax*r*μ*N)),size(trec)[1])

    #mLiveDist = zeros(Int(round(N*μ*100)))

    #distanceBM_m = zeros(Int(round(N*μ*40)),runs)

    #distanceM_m = zeros(Int(round(N*μ*40)),runs)

    mFixedM_t = zeros(size(trec)[1])

    indLiveM_t = zeros(size(trec)[1])

    depthmaxM_t = zeros(size(trec)[1])
	depthmaxBM_t = zeros(size(trec)[1])
	depthmeanM_t = zeros(size(trec)[1])
	depthmeanBM_t = zeros(size(trec)[1])

    #pc = ((Nmax/N)-1.0)/((Nmax/N)-(1/2))
    #rc = r*1.0/(1.0-(pc/2.0))

    #tc = N/20*tmax
    #μc = 20/N*μ

    for run = 1:runs

        vaf_n_t, vafB_n_t, mFixed_t, mLive, indLive_t , burden_m_t, burdenB_m_t, depthmax_t, depthmaxB_t, depthmean_t, depthmeanB_t = VAFSim.birthDeathLog(N, indStart, μ, p, Nbn, trec, r, d)

        #burdenB_m .+= mFixedM
        #burden_m .+= mFixedM

        #distanceBM_m[1:mLive+1] .+= distanceB_m
        #distanceM_m[1:mLive+1] .+= distance_m

        #println(typeof(depth_t))

        burdenBM_m_t[:,:] .+= burdenB_m_t ./ runs
        burdenM_m_t[:,:] .+= burden_m_t ./ runs

        #mLiveDist[mLive] += 1

        vafBM_n_t[:,:] .+= vafB_n_t ./ runs
        vafM_n_t[:,:] .+= vaf_n_t ./ runs

        mFixedM_t[:]  .+= mFixed_t ./runs
        indLiveM_t[:]  .+= indLive_t ./runs

        depthmaxM_t[:] .+= depthmax_t ./ runs
        depthmaxBM_t[:] .+= depthmaxB_t ./ runs
		depthmeanM_t[:] .+= depthmean_t ./ runs
        depthmeanBM_t[:] .+= depthmeanB_t ./ runs
        #println(typeof(depthM_t))
    end

    data_name = string("../data/bottleneck/dataBDlog", "_N", N, "_mu", μ, "_lin", lin,"_Nbn", Nbn, "_t", tmax, "_r", r,"_d", d, ".jld2")

    @save data_name vafM_n_t vafBM_n_t indLiveM_t trec mFixedM_t burdenM_m_t burdenBM_m_t depthmaxM_t depthmaxBM_t depthmeanM_t depthmeanBM_t #distanceM_m distanceBM_m mLiveDist

#end
