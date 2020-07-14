include("VAFsimBottleneck.jl")
using JLD2

##Script to run the simulation of neutrally evolving populations

# steps is the number of reproductive events
steps=10000
#steps_all = [200,500,1000]
# n is the number of individuals
#n = 1000
n_all = 800:100:1000
# bsize is the size of the population at the bottleneck
bsize = 100
# mu is the mutation rate
mu = 0.5
# runs is the number of simulations to average over
runs = 100

p = 0
#p_all = [0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]
rho = 1

#loop to run the simulations
for n in n_all
    VAFb_average = zeros(bsize+1)

    VAF_average = zeros(n+1)

    mu_c = mu*(1000/n)

    rho_c = rho*(n/1000)

    phi_c = 2*rho*(1-n/1000)

    p = phi_c/(phi_c+rho_c)

    steps_c = Int(round((rho_c+phi_c)*steps*rho_c))

    println(string(p))
    println(string(steps_c))
    for run = 1:runs

        genes,k,c,VAFb,VAF = RunBDAlt(steps_c,n,mu,p,bsize)

        VAFb_average = VAFb_average + VAFb
        VAF_average = VAF_average + VAF

    end

    VAF_average = VAF_average./runs
    VAFb_average = VAFb_average./runs
    data_name = string("../data/bottleneck/dataBD_varp_s", steps, "_n", n, "_b", bsize,"_mu", mu, ".jld2")

    @save data_name VAFb_average VAF_average

end
