using Plots
using JLD2
using Distributions
using Statistics

include("compoundPoisson.jl")

##Script to run the simulation of neutrally evolving populations

# tmax is the amount of time that the simulation takes
tmax = 50
#r is the rate at which replication events happen per individual
r = 1
#do args
# d is a modifier to reduce th
#d_all = 2 .^ (1:5)
d = 2^0
# n is the number of individuals
N = 200
# bsize is the size of the population at the bottleneck
Nbn = 100
# mu is the mutation rate
#μ = 2.0
μ_all = 0.2:0.2:2
# p is the likelihood of self-replacement
p = 0.0
#C is the number of somatic cells the population produces per time step
C = 0
# lin is the linear growth rate
lin = round(10^(-2),digits=4)

muta_μ = zeros(10)

popg_μ = zeros(10)
muta_t = zeros(size(trec)[1])


for i = 1:10

    μ = μ_all[i]
	data_name = string("../data/bottleneck/dataBDloglin", "_N", N, "_mu", μ, "_lin", lin,"_Nbn", Nbn, "_t", tmax, "_r", r,"_d", d, ".jld2")

	@load data_name vafM_n_t_run vafBM_n_t_run indLiveM_t trec mFixedM

	vafBM_n_t = sum(vafBM_n_t_run,dims=3) ./ size(vafBM_n_t_run)[3]
	vafM_n_t = sum(vafM_n_t_run,dims=3) ./ size(vafM_n_t_run)[3]

	fit1f_t = zeros(size(trec)[1])
	#fit1f2 = zeros(size(trec)[1])

	cor1f_t = zeros(size(trec)[1])
	#cor1f2 = zeros(size(trec)[1])

	#muta_t = zeros(size(trec)[1])

	for k = 1:size(trec)[1]


		x = N./(N-1:-1:1).-1
		y = vafM_n_t[N:-1:2,k]

		#scatter(x,y,label="Experiment",legend=:topleft,dpi=300, yaxis =:log, ylims=(10^-1,10^3))
		#plot!(1:n,dfs_fit.n_f[2:end].*mu_bn,label="ABC solution")


		cov = sum(x.*y) #- sum(x)*sum(y)

		var = sum(x.^2 ) #- sum(x)^2

		muta_t[k] = cov/var

		#b = sum(y) - a*sum(x)

		fit = muta[k].*x

		meandata = sum(y)/(N-1)
		errormean = sum( (meandata .- y).^2 )
		error1f = sum( (fit .-y).^2)
		#error1f2 = sum( ((1 ./ (1:100).^2).-(vafBM_n_t[2:101,k]./vafBM_n_t[2,k])).^2)

		fit1f_t[k] = 1 - error1f / errormean
		#fit1f2[k] = 1 -error1f2 / errormean

		cor1f_t[k] = cor( fit , y )

	end

	x = (20:50) .- mean(20:50)
	y = muta_t[20:50] .- mean(muta[20:50])
	#println(x)
	#println(y)
	cov = sum(x.*y) #- sum(x)*sum(y)

	var = sum(x.^2 ) #- sum(x)^2

	muta_μ[i] = cov/var

	x = (20:50) .- mean(20:50)
	y = indLiveM_t[20:50] .- mean(indLiveM_t[20:50])

	cov = sum(x.*y) #- sum(x)*sum(y)

	var = sum(x.^2 ) #- sum(x)^2

	popg_μ[i] = cov/var


end
plot(μ_all,muta_μ./μ_all/2, seriestype = :scatter,label = "slope of mut rate",legend=:bottomright,ylims=(0,0.012))
plot!(μ_all,popg_μ/N, seriestype = :scatter,label = "population growth")



title!(string("Relative Slope of Inferred Mutation Rate vs pop size"))
xlabel!("mu")
ylabel!("relative slope")


fig_name = string("../Figures/LogGrowth/LinR/LogLinInfMutRate_N", N, "_t", tmax)

savefig(fig_name)

b = mean(muta_t[20:50]) .- muta_μ[10].*mean(trec[20:50])
bpop = mean(indLiveM_t[20:50]) .- popg_μ[10].*mean(trec[20:50])

plot(trec,muta_t/4, seriestype = :scatter,label = "normed inferred mutations")
plot!(trec,indLiveM_t/N, seriestype = :scatter,label = "population size")
plot!(trec,(muta_μ[10].*trec .+ b)/4,label = "normed fitted mutations")
plot!(trec,(popg_μ[10].*trec .+ bpop)/N,lc = :yellow,label = "normed fitted pop",legend=:bottomright)


title!(string("Inferred Mutation Rate vs pop size over time"))
xlabel!("time")
ylabel!("mu")


fig_name = string("../Figures/LogGrowth/LinR/LogLinInfMutRate_N", N, "_t", tmax,"_mu",Int(μ))

savefig(fig_name)
