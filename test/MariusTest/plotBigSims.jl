using Plots
using JLD2
using Distributions
using Statistics

params = Dict("N initial" => 1, "N final" => 10000, "μ" => 1, "p" => 0, "λ" => 1,"growth rate" => 1,"sample size" => 100)
tStop = 80
tSaveStep = 1

data_name = string("../data/bottleneck/dataBIGexp", "_N", params["N final"], "_mu", params["μ"],"_S", params["sample size"], "_t", tStop, "_gR", params["growth rate"], ".jld2")

@load data_name times_t nLive_t vaf_n_t vafB_n_t burden_m burdenB_m

tMax = 80

Nf = params["N final"]
S = params["sample size"]
xplot = (1:1:Nf) #./ Nf
x = 1 ./ xplot
xB = 1 ./ xplot[1:S]

x15 = 2 ./ (xplot .+ xplot .^ 2)
#y = vaf_n

		#scatter(x,y,label="Experiment",legend=:topleft,dpi=300, yaxis =:log, ylims=(10^-1,10^3))



tStart = 9

vafsum_n_t = zeros(Int64,Nf+1,size(nLive_t)[1])

col(k) = 1 - ((k-tStart)/(tMax-tStart)
plot(xplot ,vaf_n_t[2:end,tStart]./vaf_n_t[2,tStart], seriestype = :scatter,yaxis=:log,ylim=(10^(-5),10^(0.2)),xlim=(0,100),label = string("pop size = ",nLive_t[tStart]),legend=false,color = RGBA(col(tStart),col(tStart),col(tStart)))

for k=tStart+1:10:tMax

	for l = 2:Nf+1
		vafsum_n_t[l,k] = sum(vaf_n_t[2:l,k])
	end
	plot!(xplot ,vaf_n_t[2:end,k]./vaf_n_t[2,k], seriestype = :scatter,label = string("pop size = ",nLive_t[k]),color = RGBA(col(k),col(k),col(k)))

end

#plot(x ,vaf_n_t[end:-1:2,tStop]./vaf_n_t[2,tStop], seriestype = :scatter,ylims=(-0.05,1.05),label = string(" Full Pop"),legend=:topleft)
#plot!(xB,vafB_n_t[end:-1:2,tStop]./vafB_n_t[2,tStop], seriestype = :scatter,label = string("Bottleneck"))
plot!(xplot,x,label = string("1/f"),color =3)
plot!(xplot,x15,label = string("2/f - 2/(f+1)"), color =5)
plot!(xplot,x.^2,label = string("1/f^2"), color =4)
#plot!(x,,label = string("1/f"))

title!(string("VAF max N = ", params["N final"]))
xlabel!("frequency")
ylabel!("# of mutations")


fig_name = string("../Figures/BigSim/BigSimExpsfinVAF_N",params["N final"], "_t", tStop,"_mu", params["μ"])

savefig(fig_name)


#y = vaf_n

		#scatter(x,y,label="Experiment",legend=:topleft,dpi=300, yaxis =:log, ylims=(10^-1,10^3))

#col(k) = 1 - ((k-9)/tMax)
#=
xfit = zeros(tMax)
x15fit = zeros(tMax)
x2fit = zeros(tMax)


for k=1:tMax-1
	xfit[k] = sum((x[(Nf-nLive_t[k]+2):Nf] .- vaf_n_t[nLive_t[k]:-1:2,k] ./ vaf_n_t[2,k]).^2)#./nLive_t[k]
	x15fit[k] = sum((x15[(Nf-nLive_t[k]+2):Nf] .- vaf_n_t[nLive_t[k]:-1:2,k] ./ vaf_n_t[2,k]).^2)#./nLive_t[k]
	#x2fit[k] = sum((x.^2 .- vaf_n_t[end:-1:2,k] ./ vaf_n_t[2,k]).^2)
end

#plot(x ,vaf_n_t[end:-1:2,tStop]./vaf_n_t[2,tStop], seriestype = :scatter,ylims=(-0.05,1.05),label = string(" Full Pop"),legend=:topleft)
#plot!(xB,vafB_n_t[end:-1:2,tStop]./vafB_n_t[2,tStop], seriestype = :scatter,label = string("Bottleneck"))


plot(1:tMax-2,xfit[2:tMax-1],label = string("1/f fit"),color =:darkgreen)
plot!(1:tMax-2,x15fit[2:tMax-1],label = string("2/f - 2/(f+1) fit"), color =:cyan)
#plot!(x,x.^2,label = string("1/f^2"), color =:green1)
#plot!(x,,label = string("1/f"))

title!(string("VAF max N = ", params["N final"]))
xlabel!("time")
ylabel!("Fit")


fig_name = string("../Figures/BigSim/BigSimFitLogVAF_N",params["N final"], "_t", tStop,"_mu", params["μ"])

savefig(fig_name)
=#
