using JLD2
include("ABC.jl")

## script to run the ABC algorithm

# trec is all the times the simulation records
trec = 1:1:10
tmax = trec[end]
#r is the rate at which replication events happen per individual
r = 1
#s is the percentage of surviving individuals for
d = 1
d_all = 2^0
# n is the number of individuals
N = 1000
indStart = 1000
# bsize is the size of the population at the bottleneck
Nbn = 100
# mu is the mutation rate
μ = 1
#μ_all = 1 #0.2:0.2:2
# runs is the number of simulations to average over
runs = 1
# p is the likelihood of self-replacement
p = 0
# C is the number of somatic cells the population is supposed to produce per time step
#C = 0
# lin is the linear growth rate of the population
lin = 0 #5*round(10^(-2),digits=4)

# n is the number of individuals
n = 100
# bsize is the size of the population at the bottleneck
#bsize = 10:5:60
# mu is the mutation rate
mu = 1

# t_bn is the fitted time based on the bottleneck pop size
#t_bn = zeros(11)
# mu_bn is the estimated mu based on the bottleneck pop size and the fitted time
#mu_bn = zeros(11)

#N_est is the estimated pop size based on the known mu
#N_est = zeros(11)
#t_est is the estimated time based on the known mu
#t_est = zeros(11)


# loop over different bottleneck sizes, run the ABC algorithm for each

#for k = 1:11
    #println(k)
    data_name = string("../data/bottleneck/dataBDpres", "_N", N, "_mu", μ, "_lin", lin,"_Nbn", Nbn, "_t", tmax, "_r", r,"_d", d, ".jld2")

    @load data_name vafM_n_t vafBM_n_t indLiveM_t trec mFixedM_t burdenM_m_t burdenBM_m_t


    t_bn,mu_bn,vafBn = ABC(vafBM_n_t[:,10],10,10,100,20,1000)

    # calculate an estimated N and time based on the known mu value
    N_est = Nbn*mu_bn/mu
    t_est = (mu_bn/mu)^2*t_bn[k]

    t = 10

    x = 1 ./ (Nbn:-1:1)
    y1 = vafBM_n_t[(Nbn+1):-1:2,t]
    y2 = vafBn[(Nbn+1):-1:2]

    plot(x,y1 ./ y1[Nbn],seriestype = :scatter,label = "sample VAF",legend=:bottomright,xlims=(0,1.05),ylims=(-0.05,1.05))
    plot!(x,y2 ./ y2[Nbn],label = "ABC derived VAF")
    plot!(x,x,label = "1/x",legend=:bottomright)

    title!(string(" t = ",t))
    xlabel!("Inverted frequency")
    ylabel!("relative M")


    fig_name = string("../Figures/Presentation/VAFPresABC_N",N, "_t", t,"_mu", μ)

    savefig(fig_name)


#end

##plots. Fitted time and mu are if we take the bottleneck population size as the
# real population size. Estimated time and N is if we know the real mu and hence
# don't need to take the bottleneck population as the assumed real population size.
# We can do the same estimate in a different direction if we know the real N or
# real time (in terms of reproduction events) instead.

# plot fitted time based on the bottleneck size

#=

h = plot(bsize,t_bn, label= "t for bottleneck size",legend=:topleft)
plot!(bsize,steps.*(bsize./n).^2, label= "x^2 estimate")
xlabel!("bottleneck size")
ylabel!("estimated time")

title!(string("Time estimates for bottleneck. True t = ", steps))
fig_name = string("../Figures/fitting/DiffBnEstT", "_n", n,"_s" ,steps, "_mu", mu)

savefig(fig_name)

#plot estimated time based on known mu
h2 = plot(bsize,t_est, label= "estimated t based on known mu",legend=:topleft)
xlabel!("bottleneck size")
ylabel!("estimated time")

title!(string("Corrected Time estimates for bottleneck. True t = ", steps))
fig_name = string("../Figures/fitting/DiffBnEstTTrue", "_n", n,"_s" ,steps, "_mu", mu)

savefig(fig_name)

#plot estimated mu based on bottleneck size and fitted time
g = plot(bsize,mu_bn, label= "t for bottleneck size",legend=:topleft)
plot!(bsize,mu.*(n./bsize), label= "1/x estimate")
xlabel!("bottleneck size")
ylabel!("estimated mu")

title!(string("Mu estimates for bottleneck. True mu = ", mu))
fig_name = string("../Figures/fitting/DiffBnEstMu", "_n", n,"_s" ,steps, "_mu", mu)

savefig(fig_name)

#plot estimated N based on known mu
j = plot(bsize,N_est, label= "N estimate based on known mu",legend=:topleft)
xlabel!("bottleneck size")
ylabel!("estimated N")

title!(string("N estimates for bottleneck. True N = ", n))
fig_name = string("../Figures/fitting/DiffBnEstNTrue", "_n", n,"_s" ,steps, "_mu", mu)

savefig(fig_name)

=#
