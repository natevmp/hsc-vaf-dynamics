using Plots
using JLD2
using Distributions

include("compoundPoisson.jl")

##Script to run the simulation of neutrally evolving populations

# tmax is the amount of time that the simulation takes
tmax = 1.5
#r is the rate at which replication events happen per individual
r = 1
# n is the number of individuals
N = 140
# bsize is the size of the population at the bottleneck
Nbn = 140
# mu is the mutation rate
μ = 1
# p is the likelihood of self-replacement
p = 0.0

#mFixedM_all = zeros(19)

#=

data_name = string("../data/bottleneck/dataBDphyldistance", "_N", N, "_mu", μ, "_p", p,"_Nbn", Nbn, "_t", tmax, "_r", r, ".jld2")

@load data_name vafM_n vafBM_n mFixedM burdenM_m burdenBM_m distanceM_m distanceBM_m mLiveDist

plot(1:(25),distanceM_m[1:25], label =string("t = ", tmax))

title!(string("Mutational distance simulation"))
xlabel!("distance")
ylabel!("# of cell comparisons")
=#
#=

for tmax =
    data_name = string("../data/bottleneck/dataBDphyldistance", "_N", N, "_mu", μ, "_p", p,"_Nbn", Nbn, "_t", tmax, "_r", r, ".jld2")

    @load data_name vafM_n vafBM_n mFixedM burdenM_m burdenBM_m distanceM_m distanceBM_m mLiveDist

    plot!(1:(25),distanceM_m[1:25], label =string("t = ", tmax))

end

=#
#=

fig_name = string("../Figures/phyldistance_N", N,"_mu", μ)

savefig(fig_name)

distanceaverage = 0
burdenaverage = 0

for k = 1:(N*μ*40)
    global distanceaverage += distanceM_m[k]*k
    global burdenaverage += burdenM_m[k]*k
end
distanceaverage = distanceaverage/(N*(N-1)/2)
burdenaverage = burdenaverage/N

burdenvar = 0

for k = 1:(N*μ*40)
    global burdenvar += (burdenM_m[k]*(k)^2)
end

burdenvar = burdenvar/(N) - burdenaverage^2

=#

#plot(1:(N*μ*20),mLiveDist[1:N*μ*20],label=string("N = ", N),legend=:topright,dpi=300)

#pc = ((Nmax/N)-1.0)/((Nmax/N)-(1/2))
#rc = r*1.0/(1.0-(pc/2.0))


N_all = 140:140:1400
N = 140

data_name = string("../data/bottleneck/dataBDtestABC", "_N", N, "_mu", μ, "_p", p,"_Nbn", Nbn, "_t", tmax, "_r", r, ".jld2")

@load data_name vafM_n vafBM_n mFixedM burdenM_m burdenBM_m

#plot(1:100,burdenM_m[1:100]/sum(burdenM_m[1:100]),label=string("simulation all"),legend=:topright,dpi=300)

plot(1:25,burdenBM_m[1:25]/sum(burdenBM_m[1:25]),label=string("simulation N =", N),legend=:topright,dpi=300)

for i = 2:10
    N = N_all[i]
data_name = string("../data/bottleneck/dataBDtestABC", "_N", N, "_mu", μ, "_p", p,"_Nbn", Nbn, "_t", tmax, "_r", r, ".jld2")

@load data_name vafM_n vafBM_n mFixedM burdenM_m burdenBM_m

#plot(1:100,burdenM_m[1:100]/sum(burdenM_m[1:100]),label=string("simulation all"),legend=:topright,dpi=300)

plot!(1:25,burdenBM_m[1:25]/sum(burdenBM_m[1:25]),label=string("simulation N =", N),legend=:topright,dpi=300)

#=
poisn = 100000

cP_rand = randComPois(r*tmax*2,μ,poisn)

cP_dist = zeros(500)

for k=1:poisn
    cP_dist[cP_rand[k]+1] += 1
end
cP_dist = cP_dist/poisn

plot!(1:100,cP_dist[1:100],label=string("compound Poisson approx"))

d = Normal(r*tmax*μ*2, sqrt(r*tmax*2*(μ+μ^2)))

plot!(1:100,pdf.(d,1:100),label=string("normal approx"))
=#

end

title!(string("Comparison simulation for different N"))
xlabel!("burden")
ylabel!("# of individuals")

fig_name = string("../Figures/burdentestall_t", tmax)

savefig(fig_name)
#=

N_all = 140:140:1400
N = 140

data_name = string("../data/bottleneck/dataBDtestABC", "_N", N, "_mu", μ, "_p", p,"_Nbn", Nbn, "_t", tmax, "_r", r, ".jld2")

@load data_name vafM_n vafBM_n mFixedM burdenM_m burdenBM_m

plot(1:10,vafM_n[2:11]/sum(vafM_n[2:11]),label=string("simulation N =", N),legend=:topright,dpi=300,yaxis=:log, ylims=(10^-5,1))


for i = 2:10
    N = N_all[i]

    data_name = string("../data/bottleneck/dataBDtestABC", "_N", N, "_mu", μ, "_p", p,"_Nbn", Nbn, "_t", tmax, "_r", r, ".jld2")

    @load data_name vafM_n vafBM_n mFixedM burdenM_m burdenBM_m

    plot!(1:10,vafBM_n[2:11]/sum(vafBM_n[2:11]),label=string("simulation N =", N),legend=:topright,dpi=300)

end

f = open("../data/Shearwater_calls_FDR0.95_all_muts.txt")

lines = readlines(f)


fline = split(lines[1],"\t")
n = size(fline)[1]-5
m = size(lines)[1]


vaf_n = zeros(Int64,n+1)
burden_ind = zeros(Int64,n)
#data = zeros(Float64,size(lines)[1],n)

for k = 2:m

	templine = split(lines[k],"\t")
	#vaf_n[parse(Int64,templine[end])+1] = vaf_n[parse(Int64,templine[end])+1] + 1
	for l = 1:n+5
		if templine[l] == "NA"
			templine[l] = "0"
		end
	end
	vaf_n[sum(parse.(Int64,templine[5:end-1]))+1] += 1
	burden_ind .+= parse.(Int64,templine[5:end-1])
end

plot!(1:10,vaf_n[2:11]/sum(vaf_n[2:11]),label="Experiment",dpi=300)

title!(string("Comparison simulation for different N"))
xlabel!("prevalence")
ylabel!("# of mutations")

fig_name = string("../Figures/vaftestall_t", tmax)

savefig(fig_name)
=#
