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
μ = 2.0
# p is the likelihood of self-replacement
p = 0.0
#C is the number of somatic cells the population produces per time step
C = 0
# lin is the linear growth rate
lin = round(10^(-2),digits=4)

data_name = string("../data/bottleneck/dataBDloglin", "_N", N, "_mu", μ, "_lin", lin,"_Nbn", Nbn, "_t", tmax, "_r", r,"_d", d, ".jld2")

@load data_name vafM_n_t_run vafBM_n_t_run indLiveM_t trec mFixedM

vafBM_n_t = sum(vafBM_n_t_run,dims=3) ./ size(vafBM_n_t_run)[3]
vafM_n_t = sum(vafM_n_t_run,dims=3) ./ size(vafM_n_t_run)[3]


plot(100 ./ (1:100),1 ./ (1:100), label =string(" 1 / f"), ylims=(10^-3,10^0),legend=:bottomright)


title!(string("VAF for logistic growth"))
xlabel!("inverse frequency")
ylabel!("# of mutations")

t = 50

scatter!(100 ./ (1:100),vafM_n_t[2:101,t]./vafM_n_t[2,t],color=:yellow, label = "Full Pop")



#plot!(1:20,vafBM_n_t_run[2:21,end,1]./vafBM_n_t_run[2,end,1], label =string("t = ", trec[end], ". Nc = ", indLiveM_t[end]), lc = RGB(trec[end]/trec[end],0,0))

plot!(100 ./ (1:100),1 ./ (1:100).^2, label =string(" 1 / f^2"), linewidth = 2)

fig_name = string("../Figures/LogGrowth/LinR/LogLinVAF_N", N,"_mu", μ, "_t", t,  "_d", d, "_lin", lin)

savefig(fig_name)
#=
fit1f = zeros(size(trec)[1])
#fit1f2 = zeros(size(trec)[1])

cor1f = zeros(size(trec)[1])
#cor1f2 = zeros(size(trec)[1])

muta = zeros(size(trec)[1])

for k = 1:size(trec)[1]


x = N./(N-1:-1:1).-1
y = vafM_n_t[N:-1:2,k]

#scatter(x,y,label="Experiment",legend=:topleft,dpi=300, yaxis =:log, ylims=(10^-1,10^3))
#plot!(1:n,dfs_fit.n_f[2:end].*mu_bn,label="ABC solution")


cov = sum(x.*y) #- sum(x)*sum(y)

var = sum(x.^2 ) #- sum(x)^2

muta[k] = cov/var

#b = sum(y) - a*sum(x)

fit = muta[k].*x

meandata = sum(y)/(N-1)
errormean = sum( (meandata .- y).^2 )
error1f = sum( (fit .-y).^2)
#error1f2 = sum( ((1 ./ (1:100).^2).-(vafBM_n_t[2:101,k]./vafBM_n_t[2,k])).^2)

fit1f[k] = 1 - error1f / errormean
#fit1f2[k] = 1 -error1f2 / errormean

cor1f[k] = cor( fit , y )

#=
meandata = sum((vafBM_n_t[2:101,k]./vafBM_n_t[2,k]))/size(trec)[1]
errormean = sum( ((meandata).-(vafBM_n_t[2:101,k]./vafBM_n_t[2,k])).^2)
error1f = sum( ((1 ./ (1:100)).-(vafBM_n_t[2:101,k]./vafBM_n_t[2,k])).^2)
error1f2 = sum( ((1 ./ (1:100).^2).-(vafBM_n_t[2:101,k]./vafBM_n_t[2,k])).^2)

fit1f[k] = 1 - error1f / errormean
fit1f2[k] = 1 -error1f2 / errormean

cor1f[k] = cor( (1 ./ (1:100)) , (vafM_n_t[2:101,k]./vafM_n_t[2,k]) )
cor1f2[k] = cor( (1 ./ (1:100).^2) , (vafM_n_t[2:101,k]./vafM_n_t[2,k]) )
=#
#plot!(trec[k],sum(vafBM_n_t_run[2:101,k,1]./vafBM_n_t_run[2,k,1]), seriestype = :scatter, lc = RGB(trec[k]/trec[end],0,0))
end

plot(trec,fit1f, seriestype = :scatter,label = "1/f",legend=:bottomright)
#plot!(trec,fit1f2, seriestype = :scatter,label = "1/f^2")

title!(string("R^2 of VAF-distribution"))
xlabel!("time")
ylabel!("R^2")

fig_name = string("../Figures/LogGrowth/LinR/LogLinFitplot_N", N,"_mu", μ, "_t", tmax, "_d", d, "_lin", lin)

#savefig(fig_name)

plot(trec,cor1f, seriestype = :scatter,label = "1/f",legend=:bottomright)
#plot!(trec,cor1f2, seriestype = :scatter,label = "1/f^2")

title!(string("Correlation of VAF-distribution"))
xlabel!("time")
ylabel!("correlation")

fig_name = string("../Figures/LogGrowth/LinR/LogLinCorplot_N", N,"_mu", μ, "_t", tmax, "_d", d, "_lin", lin)

#savefig(fig_name)

#plot(1:20,1 ./ (1:20), label =string(" 1 / f"))
#plot(trec[1],sum(vafBM_n_t_run[2:101,1,1]./vafBM_n_t_run[2,1,1]), seriestype = :scatter, lc = RGB(trec[1]/trec[end],0,0))

#=
width_t = zeros(size(trec)[1])
mutLive_t = zeros(size(trec)[1])
width2_t = zeros(size(trec)[1])
mutLive2_t = zeros(size(trec)[1])


hline([sum(1 ./ (1:100))],label = "1/f")
plot!(trec,width_t, seriestype = :scatter, color=:yellow,label = "Bottleneck",legend=:topleft)
plot!(trec,width2_t, seriestype = :scatter,label = "Full Pop",legend=:topleft)

title!(string("Width of VAF-distribution"))
xlabel!("time")
ylabel!("width of distribution")

#plot!(1:20,vafBM_n_t_run[2:21,end,1]./vafBM_n_t_run[2,end,1], label =string("t = ", trec[end], ". Nc = ", indLiveM_t[end]), lc = RGB(trec[end]/trec[end],0,0))

#plot!(1:20,1 ./ (1:20).^2, label =string(" 1 / f^2"), lc = :green, linewidth = 2)

hline!([sum(1 ./ (1:100).^2)],label = "1/f^2")

fig_name = string("../Figures/LogVarTimesAlt_N", N,"_mu", μ, "_t", tmax, "_d", d, "_C", C)

savefig(fig_name)
=#
mutLive_t = zeros(size(trec)[1])

mutLive2_t = zeros(size(trec)[1])

for k = 1:size(trec)[1]

mutLive_t[k] = sum(vafBM_n_t[2:101,k])

mutLive2_t[k] = sum(vafM_n_t[2:101,k])

end

plot(trec,indLiveM_t, seriestype = :scatter,label = "population size")
plot!(trec,mutLive_t, seriestype = :scatter, color=:yellow,label = "muts bottleneck",legend=:bottomright,yaxis=:log)
plot!(trec,mutLive2_t, seriestype = :scatter,label = "muts full pop")
#plot!(trec,2000/150 .*trec .+ 500,linewidth = 5, label = "linear")


title!(string("Population Size vs Number of Mutations"))
xlabel!("time")
ylabel!("#")

#plot!(1:20,vafBM_n_t_run[2:21,end,1]./vafBM_n_t_run[2,end,1], label =string("t = ", trec[end], ". Nc = ", indLiveM_t[end]), lc = RGB(trec[end]/trec[end],0,0))

#plot!(1:20,1 ./ (1:20).^2, label =string(" 1 / f^2"), lc = :green, linewidth = 2)
#hline!([sum(1 ./ (1:100))],label = "1/f")
#hline!([sum(1 ./ (1:100).^2)],label = "1/f^2")

fig_name = string("../Figures/LogGrowth/LinR/LogLinPopsize_N", N,"_mu", μ, "_t", tmax, "_d", d, "_lin", lin)

#savefig(fig_name)


#mutlive_1f = N*2*sum( 1 ./ (1:100))
plot(trec,indLiveM_t/N, seriestype = :scatter,label = "population size")
plot!(trec,muta/μ/2, seriestype = :scatter,label = "inferred mut rate",legend=:bottomright)
#plot!(trec,mutLive2_t, seriestype = :scatter,label = "muts full pop")
#plot!(trec,2000/150 .*trec .+ 500,linewidth = 5, label = "linear")


title!(string("Inferred Mutation Rate vs pop size"))
xlabel!("time")
ylabel!("mu")

#plot!(1:20,vafBM_n_t_run[2:21,end,1]./vafBM_n_t_run[2,end,1], label =string("t = ", trec[end], ". Nc = ", indLiveM_t[end]), lc = RGB(trec[end]/trec[end],0,0))

#plot!(1:20,1 ./ (1:20).^2, label =string(" 1 / f^2"), lc = :green, linewidth = 2)
#hline!([sum(1 ./ (1:100))],label = "1/f")
#hline!([sum(1 ./ (1:100).^2)],label = "1/f^2")

fig_name = string("../Figures/LogGrowth/LinR/LogLinInfMutRate_N", N,"_mu", μ, "_t", tmax, "_d", d, "_lin", lin)

#savefig(fig_name)
=#
