include("../vafdyn.jl")
include("../theory.jl")
using Distributions
using .VAFDyn, .Theory

##
params = Dict(
    "N initial" => 400,
    "N final" => 400,
    "μ" => 1.2 * 500 * 1e-9,
    "ρ" => 5,
    "ϕ" => 5,
    "mature time" => 5,
    "evolve time" => 20,
)
params["growth rate"] = Theory.expGrowthRateFromNT(params["N initial"], params["N final"], params["mature time"])

## ===== Discrete evolve =====

dfs = VAFDyn.DFreqspace(params["N final"])
dfs.n_f[2] = 1
@time VAFDyn.evolveCloneGrowingPop(dfs, params, params["evolve time"])
display(dfs.n_f)
println("sum: ")
println(sum(dfs.n_f[2:end-1]))
println("")

##
dfs = VAFDyn.DFreqspace(params["N final"])
dfs.n_f[2] = 1
@time vec = VAFDyn.evolveCloneGrowingPopAlt(dfs, params, params["evolve time"])
display(vec)
println("sum: ")
println(sum(vec[2:end-2]))
println("")


# ## ===== Compare with Direct Poisson prediction of occurrences
# driverProbRate(t) = params["N final"] * (2params["ρ"] + params["ϕ"]) * params["μ"] * t
# driverProbAge(t) = 1 - pdf(Poisson(driverProbRate(t)), 0)

# println("Poisson prob: ", driverProbAge(60))
# println("MC evolved prob: ", 1-vec[1])

## ===== PDE evolve =====

lVfs = 201
vfs = VAFDyn.VFreqspace(params["N final"], lVfs)
vfs.n_f[2] = 1
@time evolveCloneGrowingPop(vfs, params, params["evolve time"], 0.)
# @time evolveCloneGrowingPopReflectBC(vfs, params, params["evolve time"], 0.000000001)
# sum(vfs.n_f[2:end] .* (vfs.freqs_f[2:end] .- vfs.freqs_f[1:end-1]) * params["N final"])
# display(vfs.n_f)
#
dfs2 = VAFDyn.makeDFSfromVFS(vfs, params["N final"], normalize=false)
display(dfs2.n_f)
# sum(dfs2.n_f[2:end-1]*params["N final"])
sum(dfs2.n_f[2:end-1])
#
# # println(vfs.n_f)
# sum(vfs.n_f[1:end-1])


##
using Plots

fig = plot(yscale=:log10)
plot!(range(0, 1, length=length(vec[2:end-2])), vec[2:end-2], label="discrete MC")
plot!(vfs.freqs_f[2:end-1], vfs.n_f[2:end-1], label="PDE")
# plot!(dfs2.freqs_f[2:end-1], dfs2.n_f[2:end-1] / params["N final"], label="PDE converted", linestyle=:dash)
plot!(dfs2.freqs_f[2:end-1], dfs2.n_f[2:end-1], label="PDE converted", linestyle=:dash)

display(fig)



##
using DifferentialEquations, DiffEqOperators, SparseArrays, LinearAlgebra, Distributions, Dierckx
# using Zygote
using FiniteDifferences

##
"""Evolve growing population VAF spectrum with Moran and pure birth dynamics"""
function evolveCloneGrowingPop(vfs, params::Dict, t::Real, ϵ; alg=nothing, order=2, reltol=1e-6)
    Ni = params["N initial"]
    Nf = params["N final"]
    μ = params["μ"]
    ρ = params["ρ"]
    ϕ = params["ϕ"]
    gR = params["growth rate"]

    cappedExponentialGrowth(Ni, K, r, t) = Ni*exp(r*t)<K ? Ni*exp(r*t) : K
    cappedExpGrowthRate(n, K, r) = n < K ? r : 0
    nT(tt) = cappedExponentialGrowth(Ni, Nf, gR, tt)
    γ(n) = cappedExpGrowthRate(n, Nf, gR)

    inL = length(vfs)-2
    _m = Nf * vfs.freqs_f
    dm_m = Nf * (vfs.freqs_f[2:end] .- vfs.freqs_f[1:end-1])
    # in_m = Nf * vfs.freqs_f[2:end-1]

    # fluxInd = findall(f->f==1/Nf, vfs.freqs_f)[1] - 1
    # dϵ = (dm_m[fluxInd] + dm_m[fluxInd+1])/2 # sets intensity of flux delta function

    ∇ = CenteredDifference(1, order, dm_m, inL)
    Δ = CenteredDifference(2, order, dm_m, inL)
    A = spzeros(length(vfs), length(vfs))
    B = spzeros(length(vfs), length(vfs))
    # c = spzeros(inL)

    α(N) = 2μ*(ρ+γ(N)+ϕ/2)*N
    fA(m, N) = m<N ? - ( γ(N)*m + α(N)*(1-m/N) ) : 0
    fB(m, N) = m<N ? ρ*m*(1-m/N) + m*γ(N)/2 + α(N)*(1-m/N)/2 : 0

    # fA(m, N) = m<N ? - ( γ(N)*m ) : 0
    # fB(m, N) = m<N ? ρ*m*(1-m/N) + m*γ(N)/2 : 0

    # fL0(N) = (α(N)-N*γ(N)-2ρ)/N
    # fL1(N) = ( N*(γ(N)+2ρ) - α(N)*(N+1) )/N
    # fL2(N) = α(N)/2
    fL0(N) = (2+1/N)*α(N) - γ(N) - 2ρ
    fL1(N) = - α(N)
    # flI(t) = t==0 ? 1. : 0.
    # fluxDelta(N) = 2μ*(ρ*γ(N)+ϕ/2)*N / dϵ

    Au_m = similar(vfs.n_f)
    Bu_m = similar(vfs.n_f)
    ∇Au_m = similar(vfs.n_f)
    ΔBu_m = similar(vfs.n_f)
    ∇AuIn_m = similar(vfs.n_f[2:end-1])
    ΔBuIn_m = similar(vfs.n_f[2:end-1])


    # BC = Dirichlet0BC(Float64) # This could be changed to a Neuman BC to obtain extinction and fixation rates
    # BC = PeriodicBC(Float64)
    # BC = DirichletBC(1., 0.)
    # m0 = vfs.n_f[1]

    function step!(du_m, u_m, (A, B, Au_m, Bu_m, ∇Au_m, ΔBu_m, ∇AuIn_m, ΔBuIn_m), t)

        A[diagind(A)] .= (m -> fA(m,nT(t))).(_m)
        B[diagind(B)] .= (m -> fB(m,nT(t))).(_m)
        mul!(Au_m, A, u_m)
        mul!(Bu_m, B, u_m)
        
        # AuSpl = Spline1D(_m, Au_m)
        # BuSpl = Spline1D(_m, Bu_m)
        # ∇Au0 = forward_fdm(4, 1)(m -> AuSpl(m), 0)
        # ΔBu0 = forward_fdm(5, 2)(m -> BuSpl(m), 0)
        
        # BC = DirichletBC(u_m[1], u_m[end])
        # BC = DirichletBC(1-pdf(Poisson(α(nT(t))*t), 0), 0.)
        # BC = GeneralBC([-(∇Au0 + ΔBu0), fL0(nT(t)), fL1(nT(t)), fL2(nT(t))], [0.,1.,0.,0.], dm_m, 2)
        # BC = GeneralBC([0, fL0(nT(t)), fL1(nT(t)), fL2(nT(t))], [0.,1.,0.,0.], dm_m, 2)
        BC = GeneralBC([0, fL0(nT(t)), fL1(nT(t))], [0.,1.,0.], dm_m, 1) # zero flux boundary
        mul!(∇AuIn_m, ∇*BC, @view Au_m[2:end-1])
        mul!(ΔBuIn_m, Δ*BC, @view Bu_m[2:end-1])


        # du_m[1] = ∇Au0 + ΔBu0
        du_m[2:end-1] .= ∇AuIn_m .+ ΔBuIn_m
        # du_m[end] = 0
    end

    u0 = copy(vfs.n_f)

    prob = ODEProblem(step!, u0, (0.0, t), (A, B, Au_m, Bu_m, ∇Au_m, ΔBu_m, ∇AuIn_m, ΔBuIn_m))

    if isnothing(alg)
        sol = solve(prob, save_everystep=false)
    else
        sol = solve(prob, alg, save_everystep=false, reltol=reltol)
    end

    vfs.n_f .= sol.u[end]
end
##
function evolveCloneGrowingPopReflectBC(vfs, params::Dict, t::Real, ϵ::Real; alg=nothing, order=2, reltol=1e-6)
    Ni = params["N initial"]
    Nf = params["N final"]
    μ = params["μ"]
    ρ = params["ρ"]
    ϕ = params["ϕ"]
    gR = params["growth rate"]

    cappedExponentialGrowth(Ni, K, r, t) = Ni*exp(r*t)<K ? Ni*exp(r*t) : K
    cappedExpGrowthRate(n, K, r) = n < K ? r : 0
    nT(tt) = cappedExponentialGrowth(Ni, Nf, gR, tt)
    γ(n) = cappedExpGrowthRate(n, Nf, gR)

    inL = length(vfs)-1
    in_m = Nf * vcat(ϵ, vfs.freqs_f[2:end-1])
    dm_m = Nf.*vcat(ϵ, vfs.freqs_f[2]-ϵ, vfs.freqs_f[3:end] .- vfs.freqs_f[2:end-1])

    ∇ = CenteredDifference(1, order, dm_m, inL)
    Δ = CenteredDifference(2, order, dm_m, inL)
    A = spzeros(inL, inL)
    B = spzeros(inL, inL)

    α(N) = 2μ*(ρ+γ(N)+ϕ/2)*N
    fA(m, N) = m<N ? - ( γ(N)*m + α(N)*(1-m/N) ) : 0
    fB(m, N) = m<N ? ρ*m*(1-m/N) + m*γ(N)/2 + α(N)*(1-m/N)/2 : 0

    Au = similar(vfs.n_f[1:end-1])
    Bu = similar(vfs.n_f[1:end-1])
    ∇Au = similar(vfs.n_f[1:end-1])
    ΔBu = similar(vfs.n_f[1:end-1])

    # BC = Dirichlet0BC(Float64)
    # BC = Neumann0BC(dm_m)
    BC = RobinBC([0, 1., 0], [1., 0, 0], dm_m)

    function step!(du, u, (A, B, Au, Bu, ∇Au, ΔBu), t)
        A[diagind(A)] .= (m -> fA(m,nT(t))).(in_m)
        B[diagind(B)] .= (m -> fB(m,nT(t))).(in_m)
        mul!(Au, A, u)
        mul!(Bu, B, u)
        mul!(∇Au, ∇*BC, Au)
        mul!(ΔBu, Δ*BC, Bu)
        du .= ∇Au .+ ΔBu
    end

    u0 = zeros(Float64, length(vfs)-1)
    # u0[1] = 1/ϵ
    u0[1] = 1.

    prob = ODEProblem(step!, u0, (0.0, t), (A, B, Au, Bu, ∇Au, ΔBu))

    if isnothing(alg)
        sol = solve(prob, save_everystep=false)
    else
        sol = solve(prob, alg, save_everystep=false, reltol=reltol)
    end

    vfs.n_f[1:end-1] .= sol.u[end]
end

##
function evolveCloneGrowingPopFlux(vfs, params::Dict, t::Real, ϵ; alg=nothing, order=2, reltol=1e-6)
    Ni = params["N initial"]
    Nf = params["N final"]
    μ = params["μ"]
    ρ = params["ρ"]
    ϕ = params["ϕ"]
    gR = params["growth rate"]

    cappedExponentialGrowth(Ni, K, r, t) = Ni*exp(r*t)<K ? Ni*exp(r*t) : K
    cappedExpGrowthRate(n, K, r) = n < K ? r : 0
    nT(tt) = cappedExponentialGrowth(Ni, Nf, gR, tt)
    γ(n) = cappedExpGrowthRate(n, Nf, gR)

    inL = length(vfs)-2
    _m = Nf * vfs.freqs_f
    dm_m = Nf * (vfs.freqs_f[2:end] .- vfs.freqs_f[1:end-1])
    # in_m = Nf * vfs.freqs_f[2:end-1]

    # fluxInd = findall(f->f==1/Nf, vfs.freqs_f)[1] - 1
    # dϵ = (dm_m[fluxInd] + dm_m[fluxInd+1])/2 # sets intensity of flux delta function

    ∇ = CenteredDifference(1, order, dm_m, inL)
    Δ = CenteredDifference(2, order, dm_m, inL)
    A = spzeros(length(vfs), length(vfs))
    B = spzeros(length(vfs), length(vfs))
    # c = spzeros(inL)
    fA(m, N) = m<N ? - ( γ(N)*m + 2μ*(ρ+γ(N)+ϕ/2)*N*(1-m/N) ) : 0
    fB(m, N) = m<N ? ρ*m*(1-m/N) + m*γ(N)/2 + μ*(ρ*γ(N)+ϕ/2)*N*(1-m/N) : 0
    # fA(m, N) = m<N ? - ( γ(N)*m ) : 0
    # fB(m, N) = m<N ? ρ*m*(1-m/N) + m*γ(N)/2 : 0

    α(N) = 2μ*(ρ+γ(N)+ϕ/2)*N
    fL0(N) = (α(N)-N*γ(N)-2ρ)/N
    fL1(N) = ( N*(γ(N)+2ρ) - α(N)*(N+1) )/N
    fL2(N) = α(N)/2
    # flI(t) = t==0 ? 1. : 0.
    # fluxDelta(N) = 2μ*(ρ*γ(N)+ϕ/2)*N / dϵ

    Au_m = similar(vfs.n_f)
    Bu_m = similar(vfs.n_f)
    ∇Au_m = similar(vfs.n_f)
    ΔBu_m = similar(vfs.n_f)
    ∇AuIn_m = similar(vfs.n_f[2:end-1])
    ΔBuIn_m = similar(vfs.n_f[2:end-1])


    BC = Dirichlet0BC(Float64) # This could be changed to a Neuman BC to obtain extinction and fixation rates
    # BC = PeriodicBC(Float64)
    BC = DirichletBC(1., 0.)
    # m0 = vfs.n_f[1]

    function step!(du_m, u_m, (A, B, Au_m, Bu_m, ∇Au_m, ΔBu_m, ∇AuIn_m, ΔBuIn_m), t)

        A[diagind(A)] .= (m -> fA(m,nT(t))).(_m)
        B[diagind(B)] .= (m -> fB(m,nT(t))).(_m)
        mul!(Au_m, A, u_m)
        mul!(Bu_m, B, u_m)
        
        AuSpl = Spline1D(_m, Au_m)
        BuSpl = Spline1D(_m, Bu_m)
        # ∇Au0 = gradient(AuSpl, 0)
        # ΔBu0 = gradient(m -> gradient(BuSpl, m), 0)
        ∇Au0 = forward_fdm(4, 1)(m -> AuSpl(m), 0)
        ΔBu0 = forward_fdm(5, 2)(m -> BuSpl(m), 0)
        
        BC = DirichletBC(u_m[1], u_m[end])
        # BC = DirichletBC(1-pdf(Poisson(α(nT(t))*t), 0), 0.)
        # BC = GeneralBC([-(∇Au0 + ΔBu0), fL0(nT(t)), fL1(nT(t)), fL2(nT(t))], [0.,1.,0.,0.], dm_m, 2)
        mul!(∇AuIn_m, ∇*BC, @view Au_m[2:end-1])
        mul!(ΔBuIn_m, Δ*BC, @view Bu_m[2:end-1])


        du_m[1] = ∇Au0 + ΔBu0
        du_m[2:end-1] .= ∇AuIn_m .+ ΔBuIn_m
        du_m[end] = 0
        # du .= ∇Au .+ ΔBu .+ c
    end

    u0 = copy(vfs.n_f)

    prob = ODEProblem(step!, u0, (0.0, t), (A, B, Au_m, Bu_m, ∇Au_m, ΔBu_m, ∇AuIn_m, ΔBuIn_m))

    if isnothing(alg)
        sol = solve(prob, save_everystep=false)
    else
        sol = solve(prob, alg, save_everystep=false, reltol=reltol)
    end

    vfs.n_f .= sol.u[end]
end