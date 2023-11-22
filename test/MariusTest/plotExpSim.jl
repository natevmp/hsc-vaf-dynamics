using Plots
using JLD2
using Distributions
using Statistics

params = Dict("N initial" => 1, "N final" => 10000, "μ" => 1, "p" => 0, "λ" => 1,"growth rate" => 1,"sample size" => 100,"β" => 0)
tStop = 20
tSaveStep = 1
Nf = params["N final"]
S = params["sample size"]

data_name = string("../data/bottleneck/dataTrueExp", "_N", params["N final"], "_mu", params["μ"],"_S", params["sample size"], "_t", tStop, "_gR", params["growth rate"],"_beta", params["β"], ".jld2")

@load data_name times_t nLive_t vaf_n_t vafB_n_t burden_m burdenB_m

xplot = 1:Nf
x = 1 ./ (xplot)
xB = 1 ./ (xplot[1:S])

x2 = 2 ./ (xplot .^ 2 + xplot)



vafsum_n_t = zeros(Int64,Nf+1,size(nLive_t)[1])
tMin = 11
col(k) = 1 - ((k-tMin) / (size(nLive_t)[1]-tMin))

for l = 2:Nf+1
	vafsum_n_t[l,tMin] = sum(vaf_n_t[2:l,tMin])
end


plot(xplot ,vafsum_n_t[2:end,tMin]./vafsum_n_t[2,tMin], seriestype = :scatter,label = string(" pop Size = ", nLive_t[tMin]),legend=:bottomright,xlim=(0,100),ylim=(0,3),color = RGBA(col(tMin),col(tMin),col(tMin)))

title!(string("VAF max N = ", Nf))
xlabel!("frequency")
ylabel!("# of mutations")

for k = tMin+1:size(nLive_t)[1]-1


	for l = 2:Nf+1
		vafsum_n_t[l,k] = sum(vaf_n_t[2:l,k])
	end

	plot!(xplot ,vafsum_n_t[2:end,k]./vafsum_n_t[2,k], seriestype = :scatter,label = string(" pop Size = ", nLive_t[k]),color = RGBA(col(k),col(k),col(k)))

end
plot!(xplot,2 .- x,label = string("1/f"),color =3)

fig_name = string("../Figures/ExpSim/TrueExpVAF_N",params["N final"], "_t", tStop,"_mu", params["μ"],"_beta", params["β"])

savefig(fig_name)
