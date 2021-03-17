using Plots
using JLD2
using Distributions
using Statistics

params = Dict("N initial" => 1, "N final" => 10000, "μ" => 1, "p" => 0, "λ" => 1,"growth rate" => 1,"sample size" => 100,"β" => 0.5)
tStop = 10
tSaveStep = 1
Nf = params["N final"]
S = params["sample size"]

data_name = string("../data/bottleneck/dataTrueExp", "_N", params["N final"], "_mu", params["μ"],"_S", params["sample size"], "_t", tStop, "_gR", params["growth rate"],"_beta", params["β"], ".jld2")

@load data_name times_t nLive_t vaf_n_t vafB_n_t burden_m burdenB_m

x = 1 ./ (Nf:-1:1)
xB = 1 ./ (S:-1:1)

x2 = 2 ./ (Nf:-1:1) - 2 ./ (Nf+1:-1:2)

plot(x,2 .- x,label = string("1/f"),legend=:topright,color =:darkgreen)
#plot!(x,x2,label = string("2/f - 2/(f+1)"),color =:green)
#plot!(x,x.^2,label = string("1/f^2"),color =:green1)

title!(string("VAF max N = ", Nf))
xlabel!("f")
ylabel!("#")

vafsum_n_t = zeros(Int64,Nf+1,size(nLive_t)[1])

for k = 1:size(nLive_t)[1]-1
	col = 1 - (k / (size(nLive_t)[1]-1))

	for l = 1:Nf+1
		vafsum_n_t[l,k] = sum(vaf_n_t[1:l,k])
	end
#y = vaf_n

		#scatter(x,y,label="Experiment",legend=:topleft,dpi=300, yaxis =:log, ylims=(10^-1,10^3))

		plot!(x ,vafsum_n_t[end:-1:2,k]./vafsum_n_t[2,k], seriestype = :scatter,label = string(" pop Size = ", nLive_t[k]),color = RGBA(col,col,col))
		#plot!(xB,vafB_n[end:-1:2]./vafB_n[2], seriestype = :scatter,label = string("Bottleneck"))

#plot!(x,,label = string("1/f"))
end


fig_name = string("../Figures/ExpSim/TrueExpVAF_N",params["N final"], "_t", tStop,"_mu", params["μ"],"_beta", params["β"])

savefig(fig_name)
