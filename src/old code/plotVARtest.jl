using Plots
using JLD2
using Distributions
using Statistics

include("compoundPoisson.jl")

# trec is all the times the simulation records
trec = 1:1:1000
tmax = trec[end]
#r is the rate at which replication events happen per individual
r = 1
#s is the percentage of surviving individuals for
d = 1
d_all = 2^0
# n is the number of individuals
N = 20
indStart = 20
#N_all = 100
# bsize is the size of the population at the bottleneck
Nbn = 10
# mu is the mutation rate
μ = 1
#μ_all = 0.5:0.5:2
# runs is the number of simulations to average over
runs = 100
# p is the likelihood of self-replacement
p = 0
# C is the number of somatic cells the population is supposed to produce per time step
#C = 0
# lin is the linear growth rate of the population
lin = 1*round(10^(-3),digits=4)
#lin_all = 0.001:0.001:0.005

data_name = string("../data/bottleneck/dataBDburdenlin", "_N", N, "_mu", μ, "_lin", lin,"_Nbn", Nbn, "_t", tmax, "_r", r,"_d", d, ".jld2")

@load data_name vafM_n_t vafBM_n_t indLiveM_t trec mFixedM_t burdenM_m_t burdenBM_m_t burdenVar_t_run burdenMean_t_run #depthmaxM_t depthmaxBM_t depthmeanM_t depthmeanBM_t #distanceM_m distanceBM_m mLiveDist


bmeanlin_t = zeros(size(trec)[1])
for k = 1:Int(round(tmax*r*μ*N))
    bmeanlin_t .+= k*(burdenM_m_t[k,:])
end

bmeanlin_t = bmeanlin_t ./ indLiveM_t #.- mFixed_t_run



bvarlin_t = zeros(size(trec)[1])

for k = 1:Int(round(tmax*r*μ*N))
    bvarlin_t .+= (burdenM_m_t[k,:]*(k)^2)
end

bvarlin_t = bvarlin_t ./ indLiveM_t .- bmeanlin_t.^2

bmeanlinCor_t = bmeanlin_t .- mFixedM_t


divlin = mean(burdenVar_t_run,dims=2) ./ mean(burdenMean_t_run,dims=2)
divlin2 = mean(bvarlin_t,dims=2) ./ mean(bmeanlinCor_t,dims=2)


plot(trec,divlin ,seriestype = :scatter,label = string("corrected"),legend=:topleft)
plot!(trec,divlin2 ,seriestype = :scatter,label = string("uncorrected"),legend=:topleft)



#plot(trec,divlin ,seriestype = :scatter,label = string("var/mean"),legend=:topleft)
#plot!(trec,indLiveM_t ./ N,seriestype = :scatter,label = string("Pop Size"),legend=:topleft)


title!(string("Variance/Mean"))
xlabel!("time")
ylabel!("Var/mean")


fig_name = string("../Figures/varmeantest_N",N, "_t", tmax,"_mu", μ,"_lin",lin)

savefig(fig_name)

mu_est =


plot(1:Int(μ*tmax*5),burdenM_m_t[1:Int(μ*tmax*5),tmax]./sum(burdenM_m_t[1:Int(μ*tmax*5),tmax]),label=string("sim_lin"))
#=
poisn = 100000

cP_rand = CompoundPoisson.randComPois(r*tmax*2,μ,poisn)

cP_dist = zeros(5000)

for k=1:poisn
    cP_dist[cP_rand[k]+1] += 1
end
cP_dist = cP_dist/poisn

plot!(1:Int(μ*tmax*5),cP_dist[1:Int(μ*tmax*5)],label=string("compound Poisson approx"))
=#

#d = Normal(r*tmax*μ*2, sqrt(r*tmax*2*(μ+μ^2)))

d = Normal(mean(burdenMean_t_run),sqrt(mean(burdenVar_t_run)) )



plot!((1+Int(round(mFixedM_t[tmax]))):(Int(μ*tmax*5)+Int(round(mFixedM_t[tmax]))),pdf.(d,1:Int(μ*tmax*5)),label=string("normal approx"))

title!(string("Comparison simulation for different N"))
xlabel!("burden")
ylabel!("# of individuals")

fig_name = string("../Figures/burdentest_t", tmax,"_mu",μ,"_N", N,"_lin",lin)

savefig(fig_name)
