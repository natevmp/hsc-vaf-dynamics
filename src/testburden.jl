using Plots
using JLD2
using Distributions

include("compoundPoisson.jl")

##Script to run the simulation of neutrally evolving populations

# tmax is the amount of time that the simulation takes
tmax = 20
#r is the rate at which replication events happen per individual
r = 1
# n is the number of individuals
#N = 600
# bsize is the size of the population at the bottleneck
Nbn = 140
# mu is the mutation rate
μ = 2
# p is the likelihood of self-replacement
p = 0.0

burdenaverage_n_run = zeros(10,100)
burdenvar_n_run = zeros(10,100)

#burdenavtot = zeros(10)
#burdenerror = zeros(10)

N_all = 140:140:1400

for i = 1:10

    N = N_all[i]
    data_name = string("../data/bottleneck/dataBDerror", "_N", N, "_mu", μ, "_p", p,"_Nbn", Nbn, "_t", tmax, "_r", r, ".jld2")

    @load data_name vafM_n_run vafBM_n_run mFixedM burdenM_m_run burdenBM_m_run distanceM_m_run distanceBM_m_run mLiveDist

    #distanceaverage = 0
    for l = 1:100
        for k = 1:Int(round((tmax*r*μ*100)))
            #global distanceaverage += distanceM_m[k]*k
            global burdenaverage_n_run[i,l] += burdenBM_m_run[k,l]*k
        end
        #distanceaverage = distanceaverage/(N*(N-1)/2)
        burdenaverage_n_run[i,l] = burdenaverage_n_run[i,l]/Nbn



        for k = 1:Int(round((tmax*r*μ*100)))
            global burdenvar_n_run[i,l] += (burdenBM_m_run[k,l]*(k)^2)
        end

        burdenvar_n_run[i,l] = burdenvar_n_run[i,l]/(Nbn) .- burdenaverage_n_run[i,l]^2
    end



end

rest_n_run = (burdenaverage_n_run.^2)./(burdenvar_n_run.-burdenaverage_n_run)./tmax
muest_n_run = (burdenvar_n_run.-burdenaverage_n_run)./burdenaverage_n_run

rmean_n = zeros(10)
rsd_n = zeros(10)
mumean_n = zeros(10)
musd_n = zeros(10)

for k=1:10
    rmean_n[k]=mean(rest_n_run[k,:])
    rsd_n[k] = sqrt(var(rest_n_run[k,:]))
    mumean_n[k]=mean(muest_n_run[k,:])
    musd_n[k] = sqrt(var(muest_n_run[k,:]))
end
