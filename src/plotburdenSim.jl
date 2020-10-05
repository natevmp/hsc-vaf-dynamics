using Plots
using JLD2
using Distributions

include("compoundPoisson.jl")

##Script to run the simulation of neutrally evolving populations

# tmax is the amount of time that the simulation takes
tmax = 500
#r is the rate at which replication events happen per individual
r = 1
# n is the number of individuals
N = 500
# bsize is the size of the population at the bottleneck
Nbn = 500
# mu is the mutation rate
μ = 1
# p is the likelihood of self-replacement
p = 0.0

lin = 0

d = 1
runs=100
#mFixedM_all = zeros(19)



data_name = string("../data/bottleneck/dataBDlog", "_N", N, "_mu", μ, "_lin", lin,"_Nbn", Nbn, "_t", tmax, "_r", r,"_d", d, ".jld2")

@load data_name vafM_n_t vafBM_n_t indLiveM_t trec mFixedM_t burdenM_m_t burdenBM_m_t #distanceM_m distanceBM_m mLiveDist

#mFixed_t = mean(mFixed_t_run,dims=2)

#burden_m_t = mean(burden_m_t_run,dims=3)

#for k = 1:Int(round(tmax*r*μ*N))
#    burdenMod_m_t[k] = burdenMod_m_t
#end


#burden2_m_t_run = burden_m_t_run-mFixed_t_run

bmeanlog_t = zeros(size(trec)[1],runs)
for k = 1:Int(round(tmax*r*μ*N))
    global bmeanlog_t .+= k*(burdenM_m_t[k,:])
end

bmeanlog_t = bmeanlog_t ./ indLiveM_t #.- mFixed_t_run



bvarlog_t = zeros(size(trec)[1],runs)

for k = 1:Int(round(tmax*r*μ*N))
    global bvarlog_t .+= (burdenM_m_t[k,:,:]*(k)^2)
end

bvarlog_t = bvarlog_t ./ indLiveM_t .- bmeanlog_t.^2

bmeanlogCor_t = bmeanlog_t .- mFixedM_t

#plot(1:100,mean(burden_m_t_run[1:100,tmax,:],dims = 2),seriestype = :scatter,label = "logistic")
#plot!(trec,mean(bvar_t_run,dims=2),seriestype = :scatter,label = "var log")


plot(trec,mean(bvarlog_t,dims=2) ./ mean(bmeanlogCor_t,dims=2),seriestype = :scatter,label = "const",legend=:bottomright)
#plot!(trec,mean(bvar_t_run,dims=2),seriestype = :scatter,label = "var log")


#=
Nbn = 100

data_name = string("../data/bottleneck/dataBDlog", "_N", N, "_mu", μ, "_lin", lin,"_Nbn", Nbn, "_t", tmax, "_r", r,"_d", d, ".jld2")

@load data_name vaf_n_t_run vafB_n_t_run indLive_t_run trec mFixed_t_run burden_m_t_run burdenB_m_t_run #distanceM_m distanceBM_m mLiveDist

#burden2_m_t_run[k,:,:]


bmeanconst_t_run = zeros(size(trec)[1],runs)
for k = 1:Int(round(tmax*r*μ*N))
    global bmeanconst_t_run .+= k*burden_m_t_run[k,:,:]
end

bmeanconst_t_run = bmeanconst_t_run ./ indLive_t_run


bvarconst_t_run = zeros(size(trec)[1],runs)

for k = 1:Int(round(tmax*r*μ*N))
    global bvarconst_t_run .+= (burden_m_t_run[k,:,:]*(k)^2)
end

bvarconst_t_run = bvarconst_t_run ./ indLive_t_run .- bmeanconst_t_run.^2
=#

#plot(trec,mean(bvarconst_t_run,dims=2) ./ mean(bmeanconst_t_run,dims=2),seriestype = :scatter,label = "const")
#plot!(trec,mean(bvar_t_run,dims=2),seriestype = :scatter,label = "var const")

title!(string("Variance divided by mean of burden distribution"))
xlabel!("time")
ylabel!("var/mean")


fig_name = string("../Figures/vardivmeanburdenAlt_N",N, "_t", tmax,"_mu", μ)

savefig(fig_name)





#plot(1:100,mean(burden_m_t_run[1:100,tmax,:],dims = 2),seriestype = :scatter,label = "logistic")
#plot!(trec,mean(bvar_t_run,dims=2),seriestype = :scatter,label = "var log")


plot(trec,mean(bmeanlogCor_t,dims=2),seriestype = :scatter,label = "mean const",legend=:bottomright)
plot!(trec,mean(bvarlog_t,dims=2),seriestype = :scatter,label = "var const")
#plot(trec,mean(bmeanconst_t_run,dims=2),seriestype = :scatter,label = "mean const")
#plot!(trec,mean(bvarconst_t_run,dims=2),seriestype = :scatter,label = "var const")

title!(string("Variance and mean of burden distribution"))
xlabel!("time")
ylabel!("var/mean")


fig_name = string("../Figures/vardandmeanburdenAlt_N",N, "_t", tmax,"_mu", μ)

savefig(fig_name)

plot(1:2000,burdenM_m_t[1:2000,tmax]./sum(burdenM_m_t[1:2000,tmax]))

poisn = 100000

cP_rand = CompoundPoisson.randComPois(r*tmax*2,μ,poisn)

cP_dist = zeros(2000)

for k=1:poisn
    cP_dist[cP_rand[k]+1] += 1
end
cP_dist = cP_dist/poisn

plot!(1:2000,cP_dist[1:2000],label=string("compound Poisson approx"))

d = Normal(r*tmax*μ*2, sqrt(r*tmax*2*(μ+μ^2)))




plot!(1:2000,pdf.(d,1:2000),label=string("normal approx"))

title!(string("Comparison simulation for different N"))
xlabel!("burden")
ylabel!("# of individuals")

fig_name = string("../Figures/burdentestall_t", tmax,"_mu",μ,"_N", N)

savefig(fig_name)

normalmean = 0
for k = 1:2000
    global normalmean += k*( pdf.(d,k) )
end

#normalmean = bmeanlog_t_run ./ indLive_t_run #.- mFixed_t_run



normalvar = 0

for k = 1:2000
    global normalvar += (pdf.(d,k))*(k)^2
end

normalvar = normalvar - normalmean^2

#=

N_all = 140:140:1400
N = 140

data_name = string("../data/bottleneck/dataBDtestABC", "_N", N, "_mu", μ, "_p", p,"_Nbn", Nbn, "_t", tmax, "_r", r, ".jld2")

@load data_name vafM_n vafBM_n mFixedM burdenM_m burdenBM_m

plot(1:10,vafBM_n[2:11],label=string("simulation N =", N),legend=:topright,dpi=300,yaxis=:log, ylims=(10^-2,10^4))


for i = 2:10
    N = N_all[i]

    data_name = string("../data/bottleneck/dataBDtestABC", "_N", N, "_mu", μ, "_p", p,"_Nbn", Nbn, "_t", tmax, "_r", r, ".jld2")

    @load data_name vafM_n vafBM_n mFixedM burdenM_m burdenBM_m

    plot!(1:10,vafBM_n[2:11],label=string("simulation N =", N),legend=:topright,dpi=300)

end
#=
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
=#
#plot!(1:10,vaf_n[2:11]/sum(vaf_n[2:11]),label="Experiment",dpi=300)

title!(string("Comparison simulation for different N"))
xlabel!("prevalence")
ylabel!("# of mutations")

fig_name = string("../Figures/vaftestall_t", tmax)

savefig(fig_name)
=#
