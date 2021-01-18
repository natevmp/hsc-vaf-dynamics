using Plots
using JLD2
using Distributions

include("compoundPoisson.jl")

##Script to plot the simulation of neutrally evolving populations

# trec is all the times the simulation records
trec = 1:1:200
tmax = trec[end]
#r is the rate at which replication events happen per individual
r = 1
#s is the percentage of surviving individuals for
d = 1
d_all = 2^0
# n is the number of individuals
N = 100
#indStart = 100
# bsize is the size of the population at the bottleneck
Nbn = 10
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



data_name = string("../data/bottleneck/dataBDpres", "_N", N, "_mu", μ, "_lin", lin,"_Nbn", Nbn, "_t", tmax, "_r", r,"_d", d, ".jld2")

@load data_name vafM_n_t vafBM_n_t indLiveM_t trec mFixedM_t burdenM_m_t burdenBM_m_t #distanceM_m distanceBM_m mLiveDist


x = N*μ*2 ./ (1:Nbn)


eqProg_t = zeros(Nbn)


for k = 1


for t in trec

        test = (x .- vafBM_n_t[2:Nbn+1,t]) ./ x
        #println(test)
        if test[k] < 0.1
                eqProg_t[k] = t
                k += 1
                if k == Nbn+1
                        break
                end
        end

end

end

plot(1:Nbn,eqProg_t,seriestype=:scatter, legend = false)

title!(string("near equilibrium progression"))
xlabel!("Prevalences")
ylabel!("time")


fig_name = string("../Figures/eqProg_N",N, "_t", tmax,"_mu", μ)

savefig(fig_name)


#=
y1 = vafM_n_t[(Nbn+1):-1:2,t]
y2 = vafBM_n_t[(Nbn+1):-1:2,t]



plot(x,y1 ./ y1[Nbn],seriestype = :scatter,label = "full pop VAF",legend=:bottomright,xlims=(0,1.05),ylims=(-0.05,1.05))
#plot!(x,y2 ./ y2[Nbn],seriestype = :scatter,label = "sample VAF")
plot!(x,x,label = "1/x",legend=:bottomright)

title!(string("VAF distribution t = ",t))
xlabel!("Inverted frequency")
ylabel!("relative M")


fig_name = string("../Figures/Presentation/VAFPres1_N",N, "_t", t,"_mu", μ)

savefig(fig_name)
=#
