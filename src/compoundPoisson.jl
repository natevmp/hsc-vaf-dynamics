module CompoundPoisson

using Random, Distributions

"""
Draw samples from the compound Poisson distribution with Poisson distributed elements.
λ           Poisson parameter for distribution of sum length
μ           Poisson parameter for (Poisson) distribution of sum elements
nSamples    number of samples to generate
returns array of generated samples
"""
function randComPois(λ::Real, μ::Real, nSamples)
    Pλ = Poisson(λ)
    Pμ = Poisson(μ)
    n_s = rand(Pλ, nSamples)
    Y_s = sum.(rand.(Pμ, n_s))
end

function randComPois(λ::Real, p::Real, μ::Real, nSamples)
    ρ = λ*p
    ϕ = λ*(1-p)
    Pλ = Poisson(2ρ+ϕ)
    Pμ = Poisson(μ)
    n_s = rand(Pλ, nSamples)
    Y_s = sum.(rand.(Pμ, n_s))
end

end
