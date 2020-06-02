using JLD2
include("ABC.jl")

## script to run the ABC algorithm

# steps the number of reproductive events
steps = 1000
# n is the number of individuals
n = 100
# bsize is the size of the population at the bottleneck
bsize = 10:5:60
# mu is the mutation rate
mu = 1

# t_bn is the fitted time based on the bottleneck pop size
t_bn = zeros(11)
# mu_bn is the estimated mu based on the bottleneck pop size and the fitted time
mu_bn = zeros(11)

#N_est is the estimated pop size based on the known mu
N_est = zeros(11)
#t_est is the estimated time based on the known mu
t_est = zeros(11)


# loop over different bottleneck sizes, run the ABC algorithm for each

for k = 1:11
    println(k)
    data_name = string("C:/Users/ahw596/Documents/Julia/data/bottleneck/dataBD_s", steps, "_n", n,"_b", bsize[k],"_mu",mu, ".jld2")

    @load data_name VAFb_average VAF_average

    t_bn[k],mu_bn[k] = ABC(VAFb_average,10,10,100,20,bsize[k]^2/5)

    # calculate an estimated N and time based on the known mu value
    N_est[k] = bsize[k]*mu_bn[k]/mu
    t_est[k] = (N_est[k]/bsize[k])^2*t_bn[k]



end

##plots. Fitted time and mu are if we take the bottleneck population size as the
# real population size. Estimated time and N is if we know the real mu and hence
# don't need to take the bottleneck population as the assumed real population size.
# We can do the same estimate in a different direction if we know the real N or
# real time (in terms of reproduction events) instead.

# plot fitted time based on the bottleneck size
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
