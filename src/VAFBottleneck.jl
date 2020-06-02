include("VAFsimBottleneck.jl")
using JLD2

##Script to run the simulation of neutrally evolving populations

# steps the number of reproductive events
steps = 1000
# n is the number of individuals
n = 100
# bsize is the size of the population at the bottleneck
bsize = 45
# mu is the mutation rate
mu = 1
# runs is the number of simulations to average over
runs = 100


#loop to run the simulations
for bsize = 10:5:60
    VAFb_average = zeros(bsize+1)

    VAF_average = zeros(n+1)

    println(string(bsize))
    for run = 1:runs

        genes,k,c,VAFb,VAF = RunBDShort(steps,n,mu,bsize)

        VAFb_average = VAFb_average + VAFb
        VAF_average = VAF_average + VAF

    end

    #average the VAFs and save them
    VAF_average = VAF_average./runs
    VAFb_average = VAFb_average./runs
    data_name = string("../data/bottleneck/dataBD_s", steps, "_n", n, "_b", bsize,"_mu", mu,".jld2")

    @save data_name VAFb_average VAF_average

end
