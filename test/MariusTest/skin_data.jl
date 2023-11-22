using JLD2

include("ABC.jl")

using Plots

using Distributions

using Statistics

files = ("PD203_55y", "PD219_58y", "PD180_65y", "PD136_73y")

ages = [55,58,65,73]

ind = 4

a = zeros(ind)

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

	for k = 1:m
		templine = lines[k]

		if templine == "NA"
			templine = "0"
		end

		vafdata_n[Int(round(n*parse.(Float64,templine)))+1] += 1

	end

	x = n./(n-1:-1:2).-1
	y = vafdata_n[n:-1:3]

	#scatter(x,y,label="Experiment",legend=:topleft,dpi=300, yaxis =:log, ylims=(10^-1,10^3))


	cov = sum(x.*y) #- sum(x)*sum(y)

	var = sum(x.^2 ) #- sum(x)^2

	a[k] = cov/var

	#b = sum(y) - a*sum(x)

	fit = a[k].*x

	meandata = sum(y)/(n-1)
	errormean = sum( (meandata .- y).^2 )
	error1f = sum( (fit .-y).^2)
	#error1f2 = sum( ((1 ./ (1:100).^2).-(vafBM_n_t[2:101,k]./vafBM_n_t[2,k])).^2)

	fit1f[k] = 1 - error1f / errormean
	#fit1f2[k] = 1 -error1f2 / errormean

	cor1f[k] = cor( fit , y )
	#cor1f2[k] = cor( (1 ./ (1:100).^2) , (vafM_n_t[2:101,k]./vafM_n_t[2,k]) )

	#plot!(x,fit,label="fit",legend=:topleft)

	#fig_name = string("../Figures/skin/vaf")

	#savefig(fig_name)
end

cov = sum(ages.*a) #- sum(x)*sum(y)

var = sum(ages.^2 ) #- sum(x)^2

a2 = cov/var

#b = sum(y) - a*sum(x)

fit2 = a2.*ages

scatter(ages,a/2,label="Experiment",legend=:topleft,dpi=300)
plot!(ages,fit2/2,label="fit")

title!(string("Mutation Rate vs Age"))
xlabel!("age")
ylabel!("inferred mut rate")

fig_name = string("../Figures/skin/mutgrowth")

savefig(fig_name)

ages2 = ages .- sum(ages)/ind
fit2 = fit1f .- sum(fit1f)/ind

cov = sum(ages2.*fit2) #- sum(ages)*sum(cor1f)

var = sum(ages.^2 ) #- sum(ages)^2

afit = cov/var

bfit = sum(fit1f)/ind - afit*sum(ages)/ind

fitfit = afit.*ages .+ bfit
#=
scatter(ages,fit1f,label="R^2",legend=:bottomright,dpi=300,ylims=(0,1))
plot!(ages,fitfit,label="fit")

xlabel!("age")
ylabel!("R^2")

fig_name = string("../Figures/skin/r2")

savefig(fig_name)
=#
ages2 = ages .- sum(ages)/ind
cor2 = cor1f .- sum(cor1f)/ind

cov = sum(ages2.*cor2) #- sum(ages)*sum(cor1f)

var = sum(ages.^2 ) #- sum(ages)^2

acor = cov/var

bcor = sum(cor1f)/ind - acor*sum(ages)/ind

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
