# module driverDyn
# using DifferentialEquations
# end

using Distributions
##
using Plots
pyplot()

##


N = 50000
λ = 5
p = 0.9

ρ = λ*(1-p)
ϕ = λ*p

nBP = 6200 * 10^6
μBP = 1.2
nDrivers = 160
nBPGene = 100

μdriver = nDrivers*nBPGene*μBP / nBP

driverProbRate(t) = N * (2ρ + ϕ) * μdriver * t
driverProbAge(t) = 1 - pdf(Poisson(driverProbRate(t)), 0)

fig1 = plot(0:100, driverProbAge.(0:100), label="")
xlims!(0, 100)
ylims!(0,1)

display(fig1)