using JLD2

include("ABC.jl")

using Plots

using Distributions

using Statistics

files = ("PD203_55y", "PD219_58y", "PD180_65y", "PD136_73y")

ages = [55,58,65,73]

samples = [90, 28, 92, 24]

ind = 4

a = zeros(ind)
b = zeros(ind)

fit1f = zeros(ind)
#fit1f2 = zeros(ind)

cor1f = zeros(ind)
#cor1f2 = zeros(ind)


for k = 1:ind

	f = open( string("../data/HealthyTissueProject/Data/SkinData/", files[k], ".txt") )

	lines = readlines(f)


	fline = split(lines[1],"\t")
	m = size(lines)[1]

	n = 200

	vafdata_n = zeros(Int64,n)

	for l = 1:m
		templine = lines[l]

		if templine == "NA"
			templine = "0"
		end

		vafdata_n[Int(round(n*parse.(Float64,templine)))+1] += 1

	end

	x = n./(49:-1:3).-1
	y = vafdata_n[50:-1:4]/samples[k]

	scatter(x,y,label="Experiment",legend=:topleft,dpi=300, yaxis =:log, ylims=(10^-1,10^3))
	#plot!(1:n,dfs_fit.n_f[2:end].*mu_bn,label="ABC solution")
	xnorm = x .- mean(x)
	ynorm = y .- mean(y)

	cov = sum(xnorm.*ynorm) #- sum(x)*sum(y)

	var = sum(xnorm.^2 ) #- sum(x)^2

	a[k] = cov/var

	b[k] = mean(y) - a[k]*mean(x)

	fit = a[k].*x .+ b[k]


	errormean = sum( (mean(y) .- y).^2 )
	error1f = sum( (fit .-y).^2)
	#error1f2 = sum( ((1 ./ (1:100).^2).-(vafBM_n_t[2:101,k]./vafBM_n_t[2,k])).^2)

	fit1f[k] = 1 - error1f / errormean
	#fit1f2[k] = 1 -error1f2 / errormean

	cor1f[k] = cor( fit , y )
	#cor1f2[k] = cor( (1 ./ (1:100).^2) , (vafM_n_t[2:101,k]./vafM_n_t[2,k]) )

	plot!(x,fit,label="fit",legend=:topleft)

	fig_name = string("../Figures/skin/vaf", k)

	savefig(fig_name)

end

agesnorm = ages .- mean(ages)
anorm = a .- mean(a)

cov = sum(agesnorm.*anorm) #- sum(x)*sum(y)

var = sum(agesnorm.^2 ) #- sum(x)^2

a2 = cov/var

b2 = mean(a) - a2*mean(ages)

fit2 = a2.*ages .+ b2

scatter(ages,a,label="Experiment",legend=:topleft,dpi=300)
plot!(ages,fit2,label="fit")

title!(string("Mutation Rate vs Age"))
xlabel!("age")
ylabel!("inferred mut rate")

fig_name = string("../Figures/skin/mutgrowth")

savefig(fig_name)

fitnorm = fit1f .- mean(fit1f)

covfit = sum(agesnorm.*fitnorm) #- sum(ages)*sum(cor1f)

varfit = sum(agesnorm.^2 ) #- sum(ages)^2

afit = covfit/varfit

bfit = mean(fit1f) - afit*mean(ages)

fitfit = afit.*ages .+ bfit
#=
scatter(ages,fit1f,label="R^2",legend=:bottomright,dpi=300,ylims=(0,1))
plot!(ages,fitfit,label="fit")

title!(string("R^2 vs age"))
xlabel!("age")
ylabel!("R^2")

fig_name = string("../Figures/oeso/r2")

savefig(fig_name)
=#

cornorm = cor1f .- mean(cor1f)

covcor = sum(agesnorm.*cornorm) #- sum(ages)*sum(cor1f)

varcor = sum(agesnorm.^2 ) #- sum(ages)^2

acor = covcor/varcor

bcor = mean(cor1f) - acor*mean(ages)

fitcor = acor.*ages .+ bcor

scatter(ages,cor1f,label="cor",legend=:bottomright,dpi=300,ylims=(0,1))
scatter!(ages,fit1f,label="R^2",legend=:bottomright,dpi=300,ylims=(0,1))
plot!(ages,fitcor,label="cor fit")
plot!(ages,fitfit,label="R^2 fit")

title!(string("Correlation and R^2 vs Age"))
xlabel!("age")
ylabel!("cor/R^2")

fig_name = string("../Figures/skin/cor")

savefig(fig_name)
