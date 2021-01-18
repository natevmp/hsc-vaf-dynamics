using JLD2

include("ABC.jl")

using Plots

f = open("../data/Shearwater_calls_FDR0.95_all_muts.txt")

lines = readlines(f)


fline = split(lines[1],"\t")
n = size(fline)[1]-5
m = size(lines)[1]


vaf_n = zeros(Int64,n+1)
burden_ind = zeros(Int64,n)
#data = zeros(Float64,size(lines)[1],n)

for k = 2:m

	templine = split(lines[k],"\t")
	#vaf_n[parse(Int64,templine[end])+1] = vaf_n[parse(Int64,templine[end])+1] + 1
	for l = 1:n+5
		if templine[l] == "NA"
			templine[l] = "0"
		end
	end
	vaf_n[sum(parse.(Int64,templine[5:end-1]))+1] += 1
	burden_ind .+= parse.(Int64,templine[5:end-1])
end

burden_prev = zeros(Float64,2000)

x1 = -100:100

#spread = Normal(0,10)

for k = 1:n
	burden_prev[burden_ind[k]+1] += 1 #pdf.(spread, x1)
end
#=
muEst = (var(burden_ind)-mean(burden_ind))/mean(burden_ind)

cellDivEst = mean(burden_ind)^2/(var(burden_ind)-mean(burden_ind))

tFEst,muFEst,vafF = ABC(vaf_n,10,10,100,20,1000)

# calculate an estimated N and time based on the known mu value
NEst = n*muFEst/muEst
tEst = (muFEst/muEst)^2*tFEst

t = 10
=#
#x = (1:n)
#y1 = vaf_n[x.+1]
x = 1 ./ (n:-1:1)
y1 = vaf_n[(n+1):-1:2]
#y2 = vafF[(n+1):-1:2]

plot(x,y1 ./ y1[n],seriestype = :scatter,label = "data VAF",legend=:topleft,xlims=(0,1.05),ylims=(-0.05,1.05))
#plot!(x,y2 ./ y2[n],label = "ABC derived VAF")
plot!(x,x,label = "1/x")

title!(string("Blood Stem Cell Data VAF"))
xlabel!("frequency")
ylabel!("relative M")


fig_name = string("../Figures/Presentation/VAFPresdata")

savefig(fig_name)

#normapprox_prev = Normal(mean(burden_ind),sqrt(mean(burden_ind)*(n-1)/n))

#normapproxemp_prev = Normal(mean(burden_ind),sqrt(var(burden_ind)))


#x = 500:1500
#=
plot(x,pdf.(normapprox_prev, x),label="model-based approximation")
plot!(x,pdf.(normapproxemp_prev, x), label="empiric approximation")
plot!(x,burden_prev[x]/sum(burden_prev),label="smoothed data")

fig_name = string("../Figures/experiment/burden")

savefig(fig_name)
=#

#t_bn,mu_bn = ABC(VAF,10,10,100,20,1000)

#N_est = n*mu_bn/1.3
#t_est = (N_est/n)^2*t_bn

#params_fit = Dict("ρ" => 1, "μ" => 1, "N" => n)

#dfs_fit = VAFDyn.DFreqspace(params_fit["N"])

#VAFDyn.evolveVAF(dfs_fit, params_fit, t_bn/n, 1/(100*n))
#=
plot(1:n,VAF[2:end],label="Experiment",legend=:topleft,dpi=300)
#plot!(1:n,dfs_fit.n_f[2:end].*mu_bn,label="ABC solution")

#title!(string("VAFdistribution, est. N =", Int(round(N_est)),". Est. t =", round(t_est,digits=-2)))
xlabel!("Prevalence")
ylabel!("Number of mutations")

fig_name = string("../Figures/experiment/VAF")

savefig(fig_name)
=#
