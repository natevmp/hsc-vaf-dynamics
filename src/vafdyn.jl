module VAFDyn

export DFreqspace, CFreqspace, VFreqspace

using SparseArrays, LinearAlgebra, StaticArrays, HypergeometricFunctions, Distributions, Dierckx
using DifferentialEquations, DiffEqOperators, OrdinaryDiffEq
# using PyCall
# scp = pyimport("scipy.special")
# stat = pyimport("scipy.stats")

struct DFreqspace
    popsize::Int64
    freqs_f::Vector{Float64}
    n_f::Vector{Float64}
    function DFreqspace(popsize::Int64)
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

struct VFreqspace{F<:AbstractFloat}
    freqs_f::Vector{F}
    n_f::Vector{F}
end

function VFreqspace(N::Integer, l::Integer)
    # freqs_f = (i -> 1/N+a*i).()
    # freqs_f = Array{Float64, 1}(undef, l)
    freqs_f = zeros(Float64, l)
    freqs_f[1] = 0
    a = 2 * ( N-(l-1) )/( N*(l-1)*(l-2) )
    df_m = [1/N + (i-1)*a for i in 1:l-1]
    for m in 1:l-1
        freqs_f[1+m] = sum(df_m[1:m])
    end
    n_f = zeros(Float64, l)
    VFreqspace(freqs_f, n_f)
end

# function VFreqspace(N::Integer, l::Integer, ind1::Integer=1)
#     # freqs_f = (i -> 1/N+a*i).()
#     # freqs_f = Array{Float64, 1}(undef, l)
#     freqs_f = zeros(Float64, l)
#     freqs_f[1] = 0
#     a = 2 * ( N-(l-1) )/( N*(l-1)*(l-2) )
#     df_m = [1/N + (i-1)*a for i in 1:l-1]
#     for m in 1:l-1
#         freqs_f[1+m] = sum(df_m[1:m])
#     end
#     n_f = zeros(Float64, l)
#     VFreqspace(freqs_f, n_f)
# end

Base.length(vfs::VFreqspace) = length(vfs.freqs_f)

include("util.jl")
include("builder.jl")
include("evolver.jl")
include("scientist.jl")

end
