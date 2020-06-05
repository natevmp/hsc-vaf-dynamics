module VAFDyn

export vafMC, vafDiff, sampler, mutBurden

using StaticArrays
#using ProgressMeter
using HypergeometricFunctions
using Distributions
using PyCall
scp = pyimport("scipy.special")

struct DFreqspace
    popsize::Integer
    freqs_f::Vector{Float64}
    n_f::Vector{Float64}
    function DFreqspace(popsize::Integer)
        freqs_f = [ i/popsize for i = 0:popsize ]
        n_f = zeros(Float64, length(freqs_f))
        new(popsize, freqs_f, n_f)
    end
end

struct CFreqspace
    l::Integer
    df::Float64
    freqs_f::Vector{Float64}
    n_f::Vector{Float64}
    function CFreqspace(l::Integer)
        freqs_f = Array(range(0,1,length=l))
        n_f = zeros(Float64, l)
        df = 1.0/(l-1)
        new(l, df, freqs_f, n_f)
    end
end

include("util.jl")
include("builder.jl")
include("evolver.jl")
include("scientist.jl")


end
