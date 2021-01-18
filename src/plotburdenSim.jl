using Plots
using JLD2
using Distributions

include("compoundPoisson.jl")

##Script to run the simulation of neutrally evolving populations

# trec is all the times the simulation records
trec = 1:1:100
tmax = trec[end]
#r is the rate at which replication events happen per individual
r = 1
#s is the percentage of surviving individuals for
d = 1
d_all = 2^0
# n is the number of individuals
N = 100
#indStart = 50
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
lin = 0

#=

data_name = string("../data/bottleneck/dataBDburdenlin", "_N", N, "_mu", μ, "_lin", lin,"_Nbn", Nbn, "_t", tmax, "_r", r,"_d", d, ".jld2")

@load data_name vafM_n_t vafBM_n_t indLiveM_t trec mFixedM_t burdenM_m_t burdenBM_m_t #distanceM_m distanceBM_m mLiveDist


bmeanconst_t = zeros(size(trec)[1])
for k = 1:Int(round(tmax*r*μ*N))
    global bmeanconst_t .+= k*(burdenM_m_t[k,:])
end

bmeanconst_t = bmeanconst_t ./ indLiveM_t #.- mFixed_t_run



bvarconst_t = zeros(size(trec)[1])

for k = 1:Int(round(tmax*r*μ*N))
    global bvarconst_t .+= (burdenM_m_t[k,:]*(k)^2)
end

bvarconst_t = bvarconst_t ./ indLiveM_t .- bmeanconst_t.^2

bmeanconstCor_t = bmeanconst_t .- mFixedM_t

divconst = mean(bvarconst_t,dims=2) ./ mean(bmeanconst_t,dims=2)

=#

#plot(trec,mean(bvarconst_t,dims=2) ./ mean(bmeanconst_t,dims=2),seriestype = :scatter,label = "const")



#=
t = 15

y1 = vafM_n_t[2:11,t] ./ (100 ./ (1:10))

y2 = vafBM_n_t[2:11,t] ./ (100 ./ (1:10))

plot(1:10,y1,seriestype = :scatter,label = "full pop",legend=:bottomleft)
plot!(1:10,y2,seriestype = :scatter,label = "sample")
#plot!(trec,mean(bvar_t_run,dims=2),seriestype = :scatter,label = "var log")


title!(string("vaf (loss) in percentage"))
xlabel!("prevalence")
ylabel!("vaf in percentage")


fig_name = string("../Figures/vafloss_N",N, "_t", t,"_mu", μ)

savefig(fig_name)
=#



lin = 5*round(10^(-3),digits=4)


plot()

for N = 40:40:160



data_name = string("../data/bottleneck/dataBDburdenlin", "_N", N, "_mu", μ, "_lin", lin,"_Nbn", Nbn, "_t", tmax, "_r", r,"_d", d, ".jld2")

@load data_name vafM_n_t vafBM_n_t indLiveM_t trec mFixedM_t burdenM_m_t burdenBM_m_t #distanceM_m distanceBM_m mLiveDist


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

divlin = mean(bvarlin_t,dims=2) ./ mean(bmeanlin_t,dims=2)

plot!(trec,divlin ,seriestype = :scatter,label = string("N = ", N),legend=:topleft)

end

#plot(trec,divlin ,seriestype = :scatter,label = string("var/mean"),legend=:topleft)
#plot!(trec,indLiveM_t ./ N,seriestype = :scatter,label = string("Pop Size"),legend=:topleft)


title!(string("Variance/Mean vs Pop Size"))
xlabel!("time")
ylabel!("Var/mean or Pop Size")


fig_name = string("../Figures/varmeanvspop_N",N, "_t", tmax,"_mu", μ,"_lin",lin)

savefig(fig_name)




#plot(1:100,mean(burden_m_t_run[1:100,tmax,:],dims = 2),seriestype = :scatter,label = "logistic")
#plot!(trec,mean(bvar_t_run,dims=2),seriestype = :scatter,label = "var log")

#=
plot(trec,mean(bmeanlogCor_t,dims=2),seriestype = :scatter,label = "mean lin",legend=:bottomright)
plot!(trec,mean(bvarlog_t,dims=2),seriestype = :scatter,label = "var lin")
#plot(trec,mean(bmeanconst_t_run,dims=2),seriestype = :scatter,label = "mean const")
#plot!(trec,mean(bvarconst_t_run,dims=2),seriestype = :scatter,label = "var const")

title!(string("Variance and mean of burden distribution"))
xlabel!("time")
ylabel!("var/mean")


fig_name = string("../Figures/burdenlin_N",N, "_t", tmax,"_mu", μ)

savefig(fig_name)
=#

#for lin = 0.001:0.001:0.005

#=

plot(1:Int(μ*tmax*5),burdenM_m_t[1:Int(μ*tmax*5),tmax]./sum(burdenM_m_t[1:Int(μ*tmax*5),tmax]),label=string("sim_const"))

lin = 5*round(10^(-3),digits=4)



data_name = string("../data/bottleneck/dataBDburdenlin", "_N", N, "_mu", μ, "_lin", lin,"_Nbn", Nbn, "_t", tmax, "_r", r,"_d", d, ".jld2")

@load data_name vafM_n_t vafBM_n_t indLiveM_t trec mFixedM_t burdenM_m_t burdenBM_m_t #distanceM_m distanceBM_m mLiveDist

plot!(1:Int(μ*tmax*5),burdenM_m_t[1:Int(μ*tmax*5),tmax]./sum(burdenM_m_t[1:Int(μ*tmax*5),tmax]),label=string("sim_lin"))

poisn = 100000

cP_rand = CompoundPoisson.randComPois(r*tmax*2,μ,poisn)

cP_dist = zeros(2000)

for k=1:poisn
    cP_dist[cP_rand[k]+1] += 1
end
cP_dist = cP_dist/poisn

plot!(1:Int(μ*tmax*5),cP_dist[1:Int(μ*tmax*5)],label=string("compound Poisson approx"))

d = Normal(r*tmax*μ*2, sqrt(r*tmax*2*(μ+μ^2)))




plot!(1:Int(μ*tmax*5),pdf.(d,1:Int(μ*tmax*5)),label=string("normal approx"))

title!(string("Comparison simulation for different N"))
xlabel!("burden")
ylabel!("# of individuals")

fig_name = string("../Figures/burdentestall_t", tmax,"_mu",μ,"_N", N,"_lin",lin)

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

=#

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
