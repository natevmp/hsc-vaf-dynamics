using Plots
using JLD2
using Distributions
using Statistics

params = Dict("N initial" => 1, "N final" => 10000, "μ" => 1, "p" => 0, "λ" => 1,"growth rate" => 1,"sample size" => 100)
tStop = 40
tSaveStep = 1

data_name = string("../data/bottleneck/dataBIGexp", "_N", params["N final"], "_mu", params["μ"],"_S", params["sample size"], "_t", tStop, "_gR", params["growth rate"], ".jld2")

@load data_name times_t nLive_t vaf_n_t vafB_n_t burden_m burdenB_m


x = 1 ./ (params["N final"]:-1:1)
xB = 1 ./ (params["sample size"]:-1:1)

x2 = 2 ./ (params["N final"]:-1:1) - 2 ./ (params["N final"]+1:-1:2)
#y = vaf_n

		#scatter(x,y,label="Experiment",legend=:topleft,dpi=300, yaxis =:log, ylims=(10^-1,10^3))

plot(x ,vaf_n_t[end:-1:2,tStop]./vaf_n_t[2,tStop], seriestype = :scatter,ylims=(-0.05,1.05),label = string(" Full Pop"),legend=:topleft)
plot!(xB,vafB_n_t[end:-1:2,tStop]./vafB_n_t[2,tStop], seriestype = :scatter,label = string("Bottleneck"))
plot!(x,x,label = string("1/f"))
plot!(x,x2,label = string("2/f - 2/(f+1)"))
plot!(x,x.^2,label = string("1/f^2"))
#plot!(x,,label = string("1/f"))

title!(string("VAF max N = ", params["N final"]))
xlabel!("f")
ylabel!("#")


fig_name = string("../Figures/BigSim/BigSimExpVAF_N",params["N final"], "_t", tStop,"_mu", params["μ"])

savefig(fig_name)
