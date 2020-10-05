using Plots
using JLD2
using Distributions
using Statistics

include("compoundPoisson.jl")

##Script to run the simulation of neutrally evolving populations

# tmax is the amount of time that the simulation takes
tmax = 200
#r is the rate at which replication events happen per individual
r = 1
# n is the number of individuals
#N = 5
#N_all = 5:5:50
N_all = 10
# bsize is the size of the population at the bottleneck
Nbn = 5
# mu is the mutation rate
μ = 1
# p is the likelihood of self-replacement
p = 0.95

distmean = zeros(10)

for i = 1
    N = N_all[i]

    data_name = string("../data/bottleneck/dataBDdist", "_N", N, "_mu", μ,"_p", p,"_Nbn", Nbn, "_t", tmax, "_r", r,".jld2")

    @load data_name vafM_n_run vafBM_n_run indLiveM_t trec mFixedM distanceM_m distanceBM_m


    for k = 1:length(distanceM_m)
        distmean[i] += (k-1)*distanceM_m[k]
    end

    distmean[i] = distmean[i]/sum(distanceM_m)

    plot(0:50,distanceBM_m[1:51],label = "distance",legend=:bottomright)

    title!(string("Distance distribution"))
    xlabel!("distance")
    ylabel!("# of cell comparisons")


    fig_name = string("../Figures/phyldist_N",N, "_t", tmax)

    savefig(fig_name)

end

#=
x = log.(N_all) .- mean(log.(N_all))
y = distmean .- mean(distmean)

cov = sum(x.*y)

var = sum(x.^2 )

dista = cov/var

distb = mean(distmean) .- dista.*mean(log.(N_all))

distfit = dista .* log.(N_all) .+ distb

meandata = sum(distmean)/(10-1)
errormean = sum( (meandata .- distmean).^2 )
errorlog = sum( (distfit .-distmean).^2)
#error1f2 = sum( ((1 ./ (1:100).^2).-(vafBM_n_t[2:101,k]./vafBM_n_t[2,k])).^2)

fitlog = 1 - errorlog / errormean
#fit1f2[k] = 1 -error1f2 / errormean

corlog = cor( distfit , distmean )

plot(N_all,distmean,label = "mean distance",seriestype = :scatter,legend=:bottomright)
plot!(N_all,distfit,label = "fitted mean distance")

title!(string("Mean Distance for various N"))
xlabel!("N")
ylabel!("mean distance")


fig_name = string("../Figures/phyldistvarN", "_t", tmax)

savefig(fig_name)
=#
