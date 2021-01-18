using DiffEqOperators, OrdinaryDiffEq, SparseArrays, LinearAlgebra

"""Evolve fixed sized population with Moran dynamics as a continuous-time Markov Chain"""
function evolveVAF(dfs::DFreqspace, params::Dict, t::Number; addClones::Bool=true)
    ρ = params["ρ"]
    N = params["N"]
    μ = params["μ"]
    ϕ = params["ϕ"]
    n_f = dfs.n_f
    f_f = dfs.freqs_f

    # Calculate frequency dependent move rate
    moveRate_f = ρ * N * (f->f*(1.0 - f)).(f_f)

    fluxIn = N*2μ*(ρ+ϕ/2)
    # if mutation is turned of set incoming flux to 0
    !addClones ? fluxIn = 0 : nothing

    function step!(dn_f, n_f, p, t)
        dn_f[1] = moveRate_f[2]*n_f[2]
        dn_f[end] = moveRate_f[end-1]*n_f[end-1]
        for i in 2:length(n_f)-1
            dn_f[i] = -2moveRate_f[i]*n_f[i] +
            moveRate_f[i-1]*n_f[i-1] +
            moveRate_f[i+1]*n_f[i+1]
        end
        dn_f[2] += fluxIn
    end

    n0_f = n_f
    prob = ODEProblem(step!, n0_f, (0.0, t))
    # alg = KenCarp4() #stable for stiff PDE
    alg = TRBDF2() #stablest for stiff PDE
    # alg = Tsit5()
    # alg = BS3()
    sol = solve(prob, alg, save_everystep=false)
    # alg = Euler()
    # sol = solve(prob, alg, save_everystep=false, dt=0.001)
    dfs.n_f .= sol.u[2]
end

"""Evolve fixed sized population with Moran dynamics as a diffusion PDE"""
function evolveVAF(cfs::CFreqspace, params::Dict, t::Real; addClones::Bool=true)
    N = params["N"]
    μ = params["μ"]
    ρ = params["ρ"]
    ϕ = params["ϕ"]

    inL = cfs.l-2
    dx = cfs.df
    freqInd1 = freqToInd(1/N, dx)-1 #index to add new clones
    # if the number of gridpoints is small compared to the popsize the estimated index will be 0. In this case set it to 1.
    freqInd1==0 ? freqInd1=1 : nothing

    in_x = range(dx, step=dx, length=inL)

    X = ρ/N * sparse( 1:inL,1:inL, (x->x*(1-x)).(in_x) )
    Δ = CenteredDifference(2, 2, dx, inL)
    BC = Dirichlet0BC(Float64)
    GD = Δ*BC
    L = GD*DiffEqArrayOperator(X)

    fluxIn = N*2μ*(ρ+ϕ/2)/dx
    c = spzeros(cfs.l-2)
    addClones ? c[freqInd1]=fluxIn : nothing

    function step!(du, u, p, t)
        du .= L*u .+ c
    end

    t0 = 0.0
    t1 = t
    u0 = copy(cfs.n_f[2:end-1])

    prob = ODEProblem(step!, u0, (t0, t1))
    # alg = KenCarp4()
    alg = TRBDF2() #stablest for stiff PDE
    # alg = BS3()
    sol = solve(prob, alg, save_everystep=false)

    cfs.n_f[2:end-1] .= sol.u[end]

end

"""Evolve fixed sized population with Moran dynamics as a diffusion PDE"""
function evolveVAF(vfs::VFreqspace, params::Dict, t::Real; addClones::Bool=true)
    N = params["N"]
    μ = params["μ"]
    ρ = params["ρ"]
    ϕ = params["ϕ"]

    inL = length(vfs)-2
    dx_i = vfs.freqs_f[2:end] .- vfs.freqs_f[1:end-1]
    freqInd1 = 1 #index to add new clones

    in_x = vfs.freqs_f[2:end-1]

    X = ρ/N * sparse( 1:inL,1:inL, (x->x*(1-x)).(in_x) )
    Δ = CenteredDifference(2, 2, dx_i, inL)
    BC = Dirichlet0BC(Float64)
    GD = Δ*BC
    L = GD*DiffEqArrayOperator(X)

    fluxIn = N*2μ*(ρ+ϕ/2)/dx_i[1]
    c = spzeros(inL)
    addClones ? c[freqInd1]=fluxIn : nothing

    function step!(du, u, p, t)
        du .= L*u .+ c
    end

    t0 = 0.0
    t1 = t
    u0 = copy(vfs.n_f[2:end-1])

    prob = ODEProblem(step!, u0, (t0, t1))
    # alg = KenCarp4()
    alg = TRBDF2() #stablest for stiff PDE
    # alg = BS3()
    sol = solve(prob, alg, save_everystep=false)

    vfs.n_f[2:end-1] .= sol.u[end]

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
function evolveGrowingVAF(vfs::VFreqspace, par::Dict, t::Real, dt::Float64; addClones::Bool=true)
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
    x_f = vfs.freqs_f
    a_f = zeros(l)
    b_f = zeros(l)
    # n = N0

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
            cfs.n_f[freqToInd(1/n, cfs.df)] += 2*(ρ + γ)*n*μ*dt * 1/cfs.df
        end
    end
    
end

# ====================================================
# ================ Alternate approach ================
# ====================================================

function evolveVAFev(dfs::DFreqspace, params::Dict, t::Real)
    ρ = params["ρ"]
    N = params["N"]
    μ = params["μ"]
    ϕ = params["ϕ"]
    n_f = dfs.n_f
    f_f = dfs.freqs_f

    # Calculate frequency dependent move rate
    moveRate_f = ρ * N * (f->f*(1.0 - f)).(f_f)

    fluxIn = N*2μ*(ρ+ϕ/2)

    function step!(dy_f, y_f, p, t)
        # clone size probabilities
        dy_f[1] = moveRate_f[2]*y_f[2]
        for i in 2:N
            dy_f[i] = -2moveRate_f[i]*y_f[i] +
            moveRate_f[i-1]*y_f[i-1] +
            moveRate_f[i+1]*y_f[i+1]
        end
        dy_f[N+1] = moveRate_f[N]*y_f[N]

        # vaf expected values
        dy_f[N+2:end] .= fluxIn * y_f[1:N+1]
    end

    p0_f = zeros(Float64, N+1)
    p0_f[2] = 1
    y0_f = vcat(p0_f, n_f)

    prob = ODEProblem(step!, y0_f, (0.0, t))
    alg = TRBDF2() #stablest for stiff PDE
    # alg = KenCarp4() #stable for stiff PDE
    # alg = Tsit5()
    # alg = BS3()
    sol = solve(prob, alg, save_everystep=false)
    # sol = solve(prob, alg, save_everystep=false, dt=0.001)

    dfs.n_f .= sol.u[2][N+2:end]

end

function evolveVAFvar(dfs::DFreqspace, params::Dict, t::Real)
    ρ = params["ρ"]
    N = params["N"]
    μ = params["μ"]
    ϕ = params["ϕ"]
    nVar_f = dfs.n_f
    f_f = dfs.freqs_f

    # Calculate frequency dependent move rate
    moveRate_f = ρ * N * (f->f*(1.0 - f)).(f_f)

    fluxIn = N*2μ*(ρ+ϕ/2)
    function step!(dy_f, y_f, p, t)
        # clone size probabilities
        dy_f[1] = moveRate_f[2]*y_f[2]
        for i in 2:N
            dy_f[i] = -2moveRate_f[i]*y_f[i] +
            moveRate_f[i-1]*y_f[i-1] +
            moveRate_f[i+1]*y_f[i+1]
        end
        dy_f[N+1] = moveRate_f[N]*y_f[N]

        # vaf variances
        dy_f[N+2:end] .= (fluxIn * y_f[1:N+1]) .+ (fluxIn*μ * y_f[1:N+1].^2)
    end

    p0_f = zeros(Float64, N+1)
    p0_f[2] = 1
    y0_f = vcat(p0_f, nVar_f)

    prob = ODEProblem(step!, y0_f, (0.0, t))
    alg = TRBDF2() #stablest for stiff PDE
    # alg = KenCarp4() #stable for stiff PDE
    # alg = Tsit5()
    # alg = BS3()
    sol = solve(prob, alg, save_everystep=false)
    # sol = solve(prob, alg, save_everystep=false, dt=0.001)

    dfs.n_f .= sol.u[2][N+2:end]

end

# ====================================================
# ================== Legacy methods ==================
# ====================================================

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
function evolveVAFleg(dfs::DFreqspace, params::Dict, t::Number, dt::Float64; addClones::Bool=true)
    # legacy
    ρ = params["ρ"]
    N = params["N"]
    μ = params["μ"]
    ϕ = params["ϕ"]
    n_f = dfs.n_f
    f_f = dfs.freqs_f

    # Calculate frequency dependent move rate
    probMove(f) = f*(1. - f)
    moveRate_f = ρ * N * probMove.(f_f)

    fluxIn = N*μ*(ρ+ϕ/2)*dt

    # Evolve state space
    nC_f = copy(n_f)
    for tt in dt:dt:t
        nC_f .= n_f
        # Evolve existing clones
        stepEuler!(n_f, nC_f, moveRate_f, dt)
        # Add new clones
        if addClones
            n_f[2] += fluxIn
        end
    end
    return nothing
end

"""Evolve fixed sized population with Moran dynamics as a diffusion PDE"""
function evolveVAFleg(cfs::CFreqspace, params::Dict, t::Real, dt::Float64;
    addClones::Bool=true)
    # Legacy Euler method. Don't use.
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
        @views n_f[2:end-1] = nC_f[2:end-1] .+
                            dt*r*cd2(arg_f[2:end-1],
                                    arg_f[3:end],
                                    arg_f[1:end-2], cfs.df) / N^2
        # add clones
        if addClones
            n_f[newCloneInd] += r*μ*dt * 1/cfs.df
        end
    end
    return nothing
end

"""Evolve fixed sized population with Moran dynamics as a diffusion PDE"""
function evolveVAFalt(cfs::CFreqspace, params::Dict, t::Real, dt::Float64; addClones::Bool=true)
    # uses own implementation of finite difference operator

    N = params["N"]
    μ = params["μ"]
    ρ = params["ρ"]
    ϕ = params["ϕ"]

    l = cfs.l
    dx = cfs.df
    f_x = cfs.freqs_f

    freqInd1 = freqToInd(1/N, dx) #index to add new clones
    # if the number of gridpoints is small compared to the popsize the estimated index will be 0. In this case set it to 1.
    freqInd1==0 ? freqIn=1 : nothing

    X = ρ/N * Diagonal((x->x*(1-x)).(f_x))
    Δ = Tridiagonal(ones(l-1), -2*ones(l), ones(l-1))/dx^2
    L = Δ*X

    addClones ? fluxIn=N*μ*(ρ+ϕ/2)/dx : fluxIn = 0

    function step!(du, u, p, t)
        mul!(du, L, u)
        du[freqInd1] += fluxIn
    end

    t0 = 0.0
    t1 = t
    u0 = copy(cfs.n_f)

    prob = ODEProblem(step!, u0, (t0, t1))
    alg = KenCarp4()
    # alg = BS3()
    sol = solve(prob, alg, save_everystep=false)

    cfs.n_f .= sol.u[end]

end

function evolveVAF(dfs::DFreqspace, params::Dict, t::Number, dt::Float64; addClones::Bool=true)
    # legacy func to allow dt parameter to be passed
    evolveVAF(dfs, params, t; addClones=addClones)
end

function evolveVAFfd(cfs::CFreqspace, params::Dict, t::Real, dt::Float64; addClones::Bool=true)
    # legacy func to allow dt parameter to be passed
    evolveVAFfd(cfs, params, t; addClones=addClones)
end


function evolveVAFfd(cfs::CFreqspace, params::Dict, t::Real; addClones::Bool=true)
    evolveVAF(cfs, params, t; addClones=addClones)
end

function evolveVAFfd(vfs::VFreqspace, params::Dict, t::Real; addClones::Bool=true)
    evolveVAF(vfs, params, t; addClones=addClones)
end
