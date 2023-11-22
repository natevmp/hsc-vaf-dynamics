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
    # alg = TRBDF2() #stablest for stiff PDE
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

"""Evolve fixed sized population VAF spectrum with Moran dynamics as a diffusion PDE"""
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

function evolveGrowingVAF(dfs::DFreqspace, params::Dict, t::Real; 
    t0::Real=0, addClones::Bool=true)
    Ni = params["N initial"]
    Nf = params["N final"]
    Np = params["pure births"]
    ρ(n) = n<Np ? 0 : params["ρ"]
    ϕ(n) = n<Np ? 0 : params["ϕ"]
    μ = params["μ"]
    gR = params["growth rate"]

    nT(tt) = cappedExponentialGrowth(Ni, Nf, gR, tt)
    γ(n) = cappedExpGrowthRate(n, Nf, gR)

    mrUp(m,N) = m*( 1-m/N )*ρ(N) + m*γ(N)
    mrDw(m,N) = m*( 1-m/N )*ρ(N)
    m_m = 0:Nf
    mrUp_m = zeros(Nf+1)
    mrDw_m = zeros(Nf+1)
    p = (μ, ρ, ϕ, γ, nT, m_m, mrUp_m, mrDw_m)
    function step!(dn_m, n_m, (μ, ρ, ϕ, γ, nT, m_m, mrUp_m, mrDw_m), t)
        Nind = Int(floor(nT(t))) + 1
        # Calculate frequency dependent move rate
        mrUp_m[1:Nind] .= (m->mrUp(m,nT(t))).(@view m_m[1:Nind])
        mrDw_m[1:Nind] .= (m->mrDw(m,nT(t))).(@view m_m[1:Nind])
        # dn_m[1] = mrDw_m[2]*n_m[2]
        # dn_m[Nind] = mrUp_m[Nind-1]*n_m[Nind-1]
        @views @. dn_m[2:end-1] .= -(mrUp_m[2:end-1]+mrDw_m[2:end-1])*n_m[2:end-1] +
                mrUp_m[1:end-2]*n_m[1:end-2] +
                mrDw_m[3:end]*n_m[3:end]
        dn_m[2] += nT(t)*2μ*( ρ(nT(t))+ϕ(nT(t))/2+γ(nT(t)) )
    end
    # n0_m = SVector{length(n_m)}(dfs.n_f)
    n0_m = dfs.n_f
    prob = ODEProblem(step!, n0_m, (t0, t), p)
    # alg = KenCarp4() #stable for stiff PDE
    # alg = TRBDF2(autodiff=false) #stablest for stiff PDE
    # alg = Rodas4(autodiff=false)
    # sol = solve(prob, save_everystep=false, maxiters=5e4, autodiff=false)
    sol = solve(prob, save_everystep=false)
    dfs.n_f .= sol.u[2]
end

function evolveShrinkingVAF(dfs::DFreqspace, params::Dict, t::Real; 
    t0::Real=0, addClones::Bool=true)
    Ni = params["N initial"]
    Nf = params["N final"]
    Np = params["pure births"]
    ρ(n) = n<Np ? 0 : params["ρ"]
    ϕ(n) = n<Np ? 0 : params["ϕ"]
    μ = params["μ"]
    gR = params["growth rate"]

    # nT(tt) = cappedExponentialGrowth(Ni, Nf, gR, tt)
    # γ(n) = cappedExpGrowthRate(n, Nf, gR)
    nT(tt) = round(sizeExpGrowth2Const2ExpShrink(Ni, Nmax, gR, tConst, sR, Nlo, tt))

    γ(tt) = cappedTimeExpGrowthRate(tt, Ni, Nmax, gR)
    ζ(n, t) = t>tMax+tConst ? expShrinkRateCapped(n, sR, Nmax, Nlo) : 0
    tMax = log(Nmax/Ni)/gR

    mrUp(m,N,t) = m*( 1-m/N )*ρ(N) + m*γ(t)
    mrDw(m,N,t) = m*( 1-m/N )*ρ(N) + (m/n)*ζ(n,t)
    m_m = 0:Nf
    mrUp_m = zeros(Nf+1)
    mrDw_m = zeros(Nf+1)
    p = (μ, ρ, ϕ, γ, nT, m_m, mrUp_m, mrDw_m)
    function step!(dn_m, n_m, (μ, ρ, ϕ, γ, nT, m_m, mrUp_m, mrDw_m), t)
        Nind = Int(floor(nT(t))) + 1
        # Calculate frequency dependent move rate
        mrUp_m[1:Nind] .= (m->mrUp(m,nT(t))).(@view m_m[1:Nind])
        mrDw_m[1:Nind] .= (m->mrDw(m,nT(t))).(@view m_m[1:Nind])
        # dn_m[1] = mrDw_m[2]*n_m[2]
        # dn_m[Nind] = mrUp_m[Nind-1]*n_m[Nind-1]
        @views @. dn_m[2:end-1] .= -(mrUp_m[2:end-1]+mrDw_m[2:end-1])*n_m[2:end-1] +
                mrUp_m[1:end-2]*n_m[1:end-2] +
                mrDw_m[3:end]*n_m[3:end]
        dn_m[2] += nT(t)*2μ*( ρ(nT(t))+ϕ(nT(t))/2+γ(nT(t)) )
    end
    # n0_m = SVector{length(n_m)}(dfs.n_f)
    n0_m = dfs.n_f
    prob = ODEProblem(step!, n0_m, (t0, t), p)
    # alg = KenCarp4() #stable for stiff PDE
    # alg = TRBDF2(autodiff=false) #stablest for stiff PDE
    # alg = Rodas4(autodiff=false)
    # sol = solve(prob, save_everystep=false, maxiters=5e4, autodiff=false)
    sol = solve(prob, save_everystep=false)
    dfs.n_f .= sol.u[2]
end

"""Evolve growing population VAF spectrum with Moran and pure birth dynamics"""
function evolveGrowingVAF(vfs::VFreqspace, params::Dict, t::Real; addClones::Bool=true, alg=nothing, order=2, reltol=1e-6, growthType::String="cappedExp")
    Ni = params["N initial"]
    Nf = params["N final"]
    μ = params["μ"]
    Np = params["pure births"]
    ρ(n) = n<Np ? 0 : params["ρ"]
    ϕ(n) = n<Np ? 0 : params["ϕ"]
    gR = params["growth rate"]

    # nT(tt) = cappedExponentialGrowth(Ni, Nf, gR, tt)
    nT(tt) = if growthType=="cappedExp"
            round(cappedExponentialGrowth(Ni, Nf, gR, tt))
        elseif growthType=="logistic"
            round(logisticGrowth(Nf, gR, tt))
        end
    
    γ(n) = if growthType=="cappedExp"
            cappedExpGrowthRate(n, Nf, gR)
        elseif growthType=="logistic"
            logisticGrowthRate(n, Nf, gR)
        end

    fluxInd = findall(f->f==1/Nf, vfs.freqs_f)[1] - 1
    # fluxAmp(n) = n * 2μ * ( ρ*(1-1/(2*n)) + ϕ/2 + γ(n)*(1-1/(n+1)) )
    fluxAmp(n) = n * 2μ * ( ρ(n)*(1-1/(2*n)) + ϕ(n)/2 + γ(n) )

    inL = length(vfs)-2
    dm_i = Nf * (vfs.freqs_f[2:end] .- vfs.freqs_f[1:end-1])
    in_m = Nf * vfs.freqs_f[2:end-1]
    
    dϵ = (dm_i[fluxInd] + dm_i[fluxInd+1])/2 # sets intensity of flux delta function
    BC = Dirichlet0BC(Float64) # This could be changed to a Neuman BC to obtain extinction and fixation rates

    ∇ = CenteredDifference(1, order, dm_i, inL)
    Δ = CenteredDifference(2, order, dm_i, inL)
    c = spzeros(inL)
    A = spzeros(inL, inL)
    B = spzeros(inL, inL)
    fA(m, N) = m<N ? -γ(N)*m : 0
    fB(m, N) = m<N ? ρ(N)*(1-1/(2N))*m*(1-m/N) + m*γ(N)/2 : 0

    Au = similar(vfs.n_f[2:end-1])
    Bu = similar(vfs.n_f[2:end-1])

    ∇Au = similar(vfs.n_f[2:end-1])
    ΔBu = similar(vfs.n_f[2:end-1])

    function step!(du, u, (A, B, Au, Bu, ∇Au, ΔBu, c), t)
        A[diagind(A)] .= (m -> fA(m,nT(t))).(in_m)
        B[diagind(B)] .= (m -> fB(m,nT(t))).(in_m)
        addClones ? c[fluxInd] = fluxAmp(nT(t)) / dϵ : nothing
        DiffEqOperators.mul!(Au, A, u)
        DiffEqOperators.mul!(Bu, B, u)
        # mul!(∇Au, ∇, BC*Au)
        # mul!(ΔBu, Δ, BC*Bu)
        DiffEqOperators.mul!(∇Au, ∇*BC, Au)
        DiffEqOperators.mul!(ΔBu, Δ*BC, Bu)
        du .= ∇Au .+ ΔBu .+ c
    end

    t0 = 0.0
    t1 = t
    u0 = copy(vfs.n_f[2:end-1])

    prob = ODEProblem(step!, u0, (t0, t1), (A, B, Au, Bu, ∇Au, ΔBu, c))

    if isnothing(alg)
        sol = solve(prob, save_everystep=false)
    else
        sol = solve(prob, alg, save_everystep=false, reltol=reltol)
    end

    vfs.n_f[2:end-1] .= sol.u[end] * Nf
end

"""Evolve growing population VAF spectrum with Moran and pure birth dynamics"""
function evolveShrinkingVAF(vfs::VFreqspace, params::Dict, _T::Union{Real,AbstractVector}; addClones::Bool=true, alg=nothing, order=2, reltol=1e-6)
    Ni = params["N initial"]
    Nmax = params["N max"]
    Nlo = params["N low"]
    μ = params["μ"]
    Np = params["pure births"]
    ρ(n) = n<Np ? 0 : params["ρ"]
    ϕ(n) = n<Np ? 0 : params["ϕ"]
    gR = params["growth rate"]
    sR = params["shrink rate"]

    tConst = params["t const"]


    # nT(tt) = round(sizeExpGrowth2Const2LogShrink(Ni, Nmax, gR, tConst, sR, Nlo, tt))
    nT(tt) = round(sizeExpGrowth2Const2ExpShrink(Ni, Nmax, gR, tConst, sR, Nlo, tt))

    γ(tt) = cappedTimeExpGrowthRate(tt, Ni, Nmax, gR)
    tMax = log(Nmax/Ni)/gR
    # excess death rate during shrink
    # ζ(n, t) = t>tMax+tConst ? logShrinkRate(n, sR, Nmax, Nlo) : 0
    ζ(n, t) = t>tMax+tConst ? expShrinkRateCapped(n, sR, Nmax, Nlo) : 0

    println(findall(f->f==1/Nmax, vfs.freqs_f))
    fluxInd = findall(f->f==1/Nmax, vfs.freqs_f)[1] - 1
    fluxAmp(n, t) = n * 2μ * ( ρ(n)*(1-1/(2*n)) + ϕ(n)/2 + γ(t) )

    inL = length(vfs)-2
    dm_i = Nmax * (vfs.freqs_f[2:end] .- vfs.freqs_f[1:end-1])
    in_m = Nmax * vfs.freqs_f[2:end-1]
    
    dϵ = (dm_i[fluxInd] + dm_i[fluxInd+1])/2 # sets intensity of flux delta function
    BC = Dirichlet0BC(Float64) # This could be changed to a Neuman BC to obtain extinction and fixation rates

    ∇ = CenteredDifference(1, order, dm_i, inL)
    Δ = CenteredDifference(2, order, dm_i, inL)
    c = spzeros(inL)
    A = spzeros(inL, inL)
    B = spzeros(inL, inL)
    fA(m, N, t) = m<N ? -γ(t)*m + (m/N)*ζ(Nmax-N,t) : 0
    # conditional m<N(t) set states where m>N(t) -- which occur during growth, as
    # m = 1,2,...,Nf -- to 0;
    fB(m, N, t) = m<N ? ρ(N)*(1-1/(2N))*m*(1-m/N) + m*γ(t)/2 + (m/N)*ζ(Nmax-N,t)/2 : 0

    Au = similar(vfs.n_f[2:end-1])
    Bu = similar(vfs.n_f[2:end-1])

    ∇Au = similar(vfs.n_f[2:end-1])
    ΔBu = similar(vfs.n_f[2:end-1])

    function step!(du, u, (A, B, Au, Bu, ∇Au, ΔBu, c), t)
        A[diagind(A)] .= (m -> fA(m,nT(t),t)).(in_m)
        B[diagind(B)] .= (m -> fB(m,nT(t),t)).(in_m)
        addClones ? c[fluxInd] = fluxAmp(nT(t), t) / dϵ : nothing
        mul!(Au, A, u)
        mul!(Bu, B, u)
        # mul!(∇Au, ∇, BC*Au)
        # mul!(ΔBu, Δ, BC*Bu)
        mul!(∇Au, ∇*BC, Au)
        mul!(ΔBu, Δ*BC, Bu)
        du .= ∇Au .+ ΔBu .+ c
    end

    t0 = 0.0
    t1 = _T[end]
    u0 = copy(vfs.n_f[2:end-1])

    prob = ODEProblem(step!, u0, (t0, t1), (A, B, Au, Bu, ∇Au, ΔBu, c))

    save_everystep = false
    saveat = length(_T)>1 ? _T : []
    if isnothing(alg)
        sol = solve(prob; save_everystep, saveat)
    else
        sol = solve(prob, alg; save_everystep, saveat, reltol)
    end

    n_T = (Int ∘ floor ∘ nT).(_T)
    _T_F = (n->range(0,1,length=n+1)).(n_T)

    # println("size: ", size(sol.u))
    # println("typeof: ", typeof(sol.u))
    vfs.n_f[2:end-1] .= sol.u[end] * Nmax
    return _T_F, sol.u * Nmax
end

"""Evolve growing population VAF spectrum with Moran and pure birth dynamics"""
function evolveCloneGrowingPop(vfs::VFreqspace, params::Dict, t::Real; alg=nothing, order=2, reltol=1e-6)
    Ni = params["N initial"]
    Nf = params["N final"]
    μ = params["μ"]
    ρ = params["ρ"]
    ϕ = params["ϕ"]
    gR = params["growth rate"]

    nT(tt) = cappedExponentialGrowth(Ni, Nf, gR, tt)
    γ(n) = cappedExpGrowthRate(n, Nf, gR)

    inL = length(vfs)-2
    dm_m = Nf * (vfs.freqs_f[2:end] .- vfs.freqs_f[1:end-1])
    in_m = Nf * vfs.freqs_f[2:end-1]

    # fluxInd = findall(f->f==1/Nf, vfs.freqs_f)[1] - 1
    # dϵ = (dm_m[fluxInd] + dm_m[fluxInd+1])/2 # sets intensity of flux delta function

    ∇ = CenteredDifference(1, order, dm_m, inL)
    Δ = CenteredDifference(2, order, dm_m, inL)
    A = spzeros(inL, inL)
    B = spzeros(inL, inL)
    # c = spzeros(inL)
    fA(m, N) = m<N ? - ( γ(N)*m + 2μ*(ρ+γ(N)+ϕ/2)*N*(1-m/N) ) : 0
    fB(m, N) = m<N ? ρ*m*(1-m/N) + m*γ(N)/2 + μ*(ρ*γ(N)+ϕ/2)*N*(1-m/N) : 0

    α(N) = μ*(ρ*γ(N)+ϕ/2)*N
    fL0(N) = (α(N)+N*γ(N)-2ρ)/N
    fL1(N) = (N*(γ(N)+2ρ)-α(N)*(N+1))/N
    fL2(N) = α(N)/2
    # fluxDelta(N) = 2μ*(ρ*γ(N)+ϕ/2)*N / dϵ

    Au = similar(vfs.n_f[2:end-1])
    Bu = similar(vfs.n_f[2:end-1])
    ∇Au = similar(vfs.n_f[2:end-1])
    ΔBu = similar(vfs.n_f[2:end-1])

    BC = Dirichlet0BC(Float64) # This could be changed to a Neuman BC to obtain extinction and fixation rates
    function step!(du, u, (A, B, Au, Bu, ∇Au, ΔBu), t)
        A[diagind(A)] .= (m -> fA(m,nT(t))).(in_m)
        B[diagind(B)] .= (m -> fB(m,nT(t))).(in_m) 
        BC = GeneralBC([0, fL0(nT(t)), fL1(nT(t)), fL2(nT(t))], [0.,1.,0.,0.], dm_m, 2)
        mul!(Au, A, u)
        mul!(Bu, B, u)
        mul!(∇Au, ∇*BC, Au)
        mul!(ΔBu, Δ*BC, Bu)
        du .= ∇Au .+ ΔBu
        # du .= ∇Au .+ ΔBu .+ c
    end

    u0 = copy(vfs.n_f[2:end-1])

    prob = ODEProblem(step!, u0, (0.0, t), (A, B, Au, Bu, ∇Au, ΔBu))

    if isnothing(alg)
        sol = solve(prob, save_everystep=false)
    else
        sol = solve(prob, alg, save_everystep=false, reltol=reltol)
    end

    vfs.n_f[2:end-1] .= sol.u[end]
end


function evolveCloneGrowingPop(dfs::DFreqspace, params::Dict, t::Number; alg=nothing, order=2, reltol=1e-6)
    Ni = params["N initial"]
    Nf = params["N final"]
    ρ = params["ρ"]
    ϕ = params["ϕ"]
    μ = params["μ"]
    gR = params["growth rate"]

    nT(tt) = cappedExponentialGrowth(Ni, Nf, gR, tt)
    γ(n) = cappedExpGrowthRate(n, Nf, gR)

    tUp(m,N) = m<N ? m*( 1-m/N )*ρ + m*γ(N) + 2μ*(1-m/N)*(ρ+γ(N)+ϕ/2)*N : 0
    # tUp(m,N) = m<N ? m*( 1-m/N )*ρ + m*γ(N) : 0
    tDw(m,N) = m<N ? m*( 1-m/N )*ρ : 0
    tSt(m,N) = tUp(m,N) + tDw(m,N)

    T = sparse(zeros(length(dfs.n_f), length(dfs.n_f)))
    parEvo = (nT, T)
    function step!(dp_m, p_m, (nT, T), t)
        Nt = Int(floor(nT(t)))
        T[diagind(T,-1)[1:Nt]] .= tUp.(0:Nt-1, nT(t))
        T[diagind(T)[1:Nt+1]] .= - tSt.(0:Nt, nT(t)) 
        T[diagind(T,+1)[1:Nt]] .= tDw.(1:Nt, nT(t))
        mul!(dp_m, T, p_m)
    end

    p0_m = dfs.n_f
    prob = ODEProblem(step!, p0_m, (0.0, t), parEvo)
    sol = solve(prob, save_everystep=false)
    dfs.n_f .= sol.u[2]
end


function evolveCloneGrowingPopAlt(dfs::DFreqspace, params::Dict, t::Number; alg=nothing, order=2, reltol=1e-6)
    Ni = params["N initial"]
    Nf = params["N final"]
    ρ = params["ρ"]
    ϕ = params["ϕ"]
    μ = params["μ"]
    gR = params["growth rate"]

    nT(tt) = cappedExponentialGrowth(Ni, Nf, gR, tt)
    γ(n) = cappedExpGrowthRate(n, Nf, gR)

    # tUp(m,N) = m<N ? m*( 1-m/N )*ρ + m*γ(N) : 0
    tUp(m,N) = m<N ? m*( 1-m/N )*ρ + m*γ(N) + 2μ*(1-m/N)*(ρ+γ(N)+ϕ/2)*N : 0
    tDw(m,N) = m<N ? m*( 1-m/N )*ρ : 0
    tSt(m,N) = tUp(m,N) + tDw(m,N)

    T = sparse(zeros(length(dfs.n_f)-1, length(dfs.n_f)-1))
    parEvo = (nT, T)
    function step!(dp_m, p_m, (nT, T), t)
        Nt = Int(floor(nT(t)))

        dp_m[1] = - tSt(0, nT(t)) * p_m[1]
        dp_m[2] = tUp(0, nT(t))*(p_m[end] + p_m[1]) - tSt(1, nT(t))*p_m[2] + tDw(2, nT(t))*p_m[3]
        for m in 2:Nt-1
            dp_m[1+m] = tUp(m-1, nT(t))*p_m[m] - tSt(m, nT(t))*p_m[m+1] + tDw(m+1, nT(t))*p_m[m+2]
        end
        dp_m[1+Nt] = tUp(Nt-1, nT(t))*p_m[Nt]
        dp_m[end] = - tSt(0, nT(t))*p_m[end] + tDw(1, nT(t)) * p_m[2]
    end

    p0_m = append!(deepcopy(dfs.n_f), 0)
    prob = ODEProblem(step!, p0_m, (0.0, t), parEvo)
    sol = solve(prob, save_everystep=false)
    # dfs.n_f .= sol.u[2]
    return sol.u[2]
end

# ==== Growth Functions ====
exponentialGrowth(Ni, r, t) = Ni*exp(r*t)
linearGrowth(Ni, r, t) = Ni + Ni*r*t

function cappedExponentialGrowth(Ni, K, r, t)
	return Ni*exp(r*t)<K ? Ni*exp(r*t) : K
end

function logisticGrowth(K, r, t)
    K / ( 1 + (K-1)*exp(-r*t) )
end

function logisticShrink(NMax, NMin, r, t)
    NMax - (logisticGrowth(NMax-NMin, r, t) - 1)
end

function expShrinkCapped(NMax, NMin, r, t)
    NMax - cappedExponentialGrowth(1, NMax-NMin, r, t)
end

function logShrinkRate(n, r, nMax, nMin)
    - r*(nMax-n)*(1-(nMax-n)/(nMax-nMin))
end

function expShrinkRateCapped(n, r, nMax, nMin)
    n > nMin ? r * (nMax-n) : 0
end

function sizeExpGrowth2Const2LogShrink(Ni, NMax, rG, tConst, rS, Nf, t)
    if t < log(NMax/Ni)/rG
        return exponentialGrowth(Ni, rG, t)
    elseif t < log(NMax/Ni)/rG+tConst
        return NMax
    else
        return logisticShrink(NMax, Nf, rS, t - (log(NMax/Ni)/rG+tConst))
    end
end

function sizeExpGrowth2Const2ExpShrink(Ni, NMax, rG, tConst, rS, Nf, t)
    if t < log(NMax/Ni)/rG
        return exponentialGrowth(Ni, rG, t)
    elseif t < log(NMax/Ni)/rG+tConst
        return NMax
    else
        return expShrinkCapped(NMax, Nf, rS, t - (log(NMax/Ni)/rG+tConst))
    end
end

function cappedExpGrowthRate(n, K, r)
    n < K ? r : 0
end

function cappedTimeExpGrowthRate(t, Ni, NMax, gR)
    t < log(NMax/Ni)/gR ? gR : 0
end

function logisticGrowthRate(n, NMax, gR)
    gR*(1-n/NMax)
end

function exponentialToLinearGrowth(Ni, K, rExp, rLin, t)
    tK = log(K/Ni)/rExp
    t<tK ? n=exponentialGrowth(Ni, rExp, t) : n=linearGrowth(K, rLin, t-tK)
    return n
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
