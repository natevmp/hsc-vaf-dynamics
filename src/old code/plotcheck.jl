using Plots
using JLD2

# steps is the number of reproductive events
steps = 1000
#steps_all = [5,10,20,50,100,200,500,1000]
# n is the number of individuals
#n = 1000
n_all = 100:100:900
# bsize is the size of the population at the bottleneck
bsize = 100
# mu is the mutation rate
mu = 1
# p is the percentage of asymmetric divisions
p=0
#p_all = [0.5,0.4,0.3,0.2,0.1]
# runs is the number of simulations to average over

runs = 100

#rho = 1

diff = zeros(9)

data_name2 = string("../data/bottleneck/dataBD_varp_s", steps, "_n", 1000, "_b", bsize,"_mu", mu, ".jld2")

@load data_name2 VAFb_average VAF_average

VAF_check = VAF_average
VAFb_check = VAFb_average

for k=1:9
    n = n_all[k]
    data_name = string("../data/bottleneck/dataBD_varp_s", steps, "_n", n, "_b", bsize,"_mu", mu, ".jld2")

    @load data_name VAFb_average VAF_average

    VAF_normal = VAF_average
    VAFb_normal = VAFb_average

    diff[k] = sqrt(sum((VAFb_check[2:end] - VAFb_normal[2:end]).^2))./sum(VAFb_check[2:end])
    println(sqrt(sum((VAFb_check[2:end] - VAFb_normal[2:end]).^2))./sum(VAFb_check[2:end]))
    plot(1:bsize,VAFb_normal[2:end],label="VAF",legend=:topright,dpi=300,yaxis=:log)
    plot!(1:bsize,VAFb_check[2:end],label="VAF check")
    ylims!((maximum(VAFb_normal)*(10^(-6)),maximum(VAFb_normal)*1.5))
    xlims!((0,findall(x->x==0,VAFb_normal)[1]))
    title!(string("Comparison"))
    xlabel!("Prevalence")
    ylabel!("Number of mutations")

    fig_name = string("../Figures/Varp_n", n, "_s", steps)

    savefig(fig_name)

end

plot(n_all,diff,legend=false,dpi=300)
#plot!(1:,VAF_check[2:end],label="VAF check")

title!(string("difference for varying p and n VAFs"))
xlabel!("n")
ylabel!("Difference")

fig_name = string("../Figures/CheckbP_s", steps)

savefig(fig_name)
