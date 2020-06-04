using FiniteDiff

"""Evolve vector for timestep using Euler method"""
function stepEuler!(n_f::AbstractVector, nC_f::AbstractVector, moveRate_f::AbstractVector, dt::Real)
    # evolve boundaries
    n_f[1] = nC_f[1] + moveRate_f[2]dt*nC_f[2]
    n_f[end] = nC_f[end] + moveRate_f[end-1]dt*nC_f[end-1]
    # evolve body
    for i in 2:length(n_f)-1
        n_f[i] = (1-2moveRate_f[i]dt)*nC_f[i] +
        moveRate_f[i-1]dt*nC_f[i-1] +
        moveRate_f[i+1]dt*nC_f[i+1]
    end
    return nothing
end

"""Evolve fixed sized population with Moran dynamics as a continuous-time Markov Chain"""
function evolveVAF(dfs::DFreqspace, params::Dict, t::Number, dt::Float64; addClones::Bool=true)
    ρ = params["ρ"]
    N = params["N"]
    μ = params["μ"]
    ϕ = params["ϕ"]
    n_f = dfs.n_f
    f_f = dfs.freqs_f

    # Calculate frequency dependent move rate
    probMove(f) = f*(1. - f)
    moveRate_f = ρ * N * probMove.(f_f)

    # Evolve state space
    nC_f = copy(n_f)
    for tt in dt:dt:t
        nC_f .= n_f
        # Evolve existing clones
        stepEuler!(n_f, nC_f, moveRate_f, dt)
        # Add new clones
        if addClones
            n_f[2] += N*μ*(ρ+ϕ/2)*dt
        end
    end
    return nothing
end

"""Evolve fixed sized population with Moran dynamics as a diffusion PDE"""
function evolveVAF(cfs::CFreqspace, params::Dict, t::Real, dt::Float64; addClones::Bool=true)
    ρ = params["ρ"]
    N = params["N"]
    μ = params["μ"]
    #event rate
    r = ρ * N
    #preallocate temp arrays
    n_f = cfs.n_f
    nC_f = zeros(length(n_f))
    f_f = cfs.freqs_f
    arg_f = zeros(length(n_f))
    g_f = map(x -> x*(1-x), f_f)
    newCloneInd = freqToInd(1/N, cfs.df)

    for tt in dt:dt:t
        nC_f .= n_f
        arg_f .= g_f .* nC_f
        # diffuse step euler
        cfs.n_f[2:end-1] = nC_f[2:end-1] .+
                            dt*r*cd2(arg_f[2:end-1], arg_f[3:end], arg_f[1:end-2], cfs.df) / N^2
        # add clones
        if addClones
            cfs.n_f[1 + newCloneInd] += r*μ*dt * 1/cfs.df
        end
    end
    return nothing
end


"""Evolve a single clone in fixed populations with Moran dynamics"""
function evolveCloneKimura(fs_f, p, t, n, cutoff::Integer)
    # native 2F1 is not very accurate; use pycall instead
    # function _1T(x, i)
    #     1/2*(i+1)*(i+2) * _₂F₁(i+3, -i, 2, (1-x)/2)
    # end
    # PyCalled hyp2f1 performs better than native 2F1
    function _1T(x, i)
        1/2*(i+1)*(i+2) * scp.hyp2f1(i+3, -i, 2, (1-x)/2)
    end
    function f(x, p, t, c::Integer)
        sum( [
                4*(2i+1)*p*(1-p)/(i*(i+1)) * _1T(1-2p, i-1) * _1T(1-2x, i-1) * exp( -i*(i+1)*t/(n) )
            for i in 1:c ] )
    end
    return map( x -> f(x, p, t, cutoff), fs_f)
end

"""apply Fokker-Planck operator"""
function fpOp(p::Array{T}, a::Array{T}, b::Array{T}, dx::Float64) where T
    return - fd1(a[2:end].*p[2:end], dx) .+ (1/2)*cd2(b.*p, dx)
end

"""Evolve growing populations with Moran and pure birth dynamics"""
function evolveGrowingVAF(cfs::CFreqspace, par::Dict, t::Real, dt::Float64; addClones::Bool=true)
    #event rates
    ρ = par["ρ"]
    γ = par["γ"]
    μ = par["μ"]
    N0 = par["N0"]

    #preallocate
    n_f = cfs.n_f
    l = length(n_f)
    nC_f = zeros(l)
    arg_f = zeros(l)
    x_f = cfs.freqs_f
    a_f = zeros(l)
    b_f = zeros(l)
    n = N0

    bFP(x, n) = n*( 2ρ/(n^2) + γ/(n+1)^2 )*x*(1-x)

    for tt in dt:dt:t
        nC_f .= n_f
        # evolve clones
        n = N0*exp(γ*tt)
        a_f .= 0
        b_f .= map(x -> bFP(x, n), x_f)
        n_f[2:end-1] = nC_f[2:end-1] .+ dt*fpOp(nC_f, a_f, b_f, cfs.df)
        # add clones
        if addClones
            cfs.n_f[1 + freqToInd(1/n, cfs.df)] += 2*(ρ + γ)*n*μ*dt * 1/cfs.df
        end
    end


end
