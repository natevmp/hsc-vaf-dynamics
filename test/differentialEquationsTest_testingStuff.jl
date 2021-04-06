using DifferentialEquations, DiffEqOperators
using SparseArrays


# function f(u,p,t)
#     # 1.01*u
#     u = 1.
# end
# u0 = 1/2
# tspan = (0.0,1.0)
# prob = ODEProblem(f,u0,tspan)
# sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)

# using Plots
# plot(sol,linewidth=5,title="Solution to the linear ODE with a thick line",
#      xaxis="Time (t)",yaxis="u(t) (in μm)",label="My Thick Line!") # legend=false
# plot!(sol.t, t->0.5*exp(1.01t),lw=3,ls=:dash,label="True Solution!")



# Ni = params["N initial"]
# Nf = params["N final"]
# μ = params["μ"]
# ρ = params["ρ"]
# ϕ = params["ϕ"]
# gR = params["growth rate"]

N = 100

inL = N-1
dx_i = ones(N)
in_x = 1.:N-1

BC = Dirichlet0BC(Float64)
# ∇ = UpwindDifference(1, 1, dx_i, inL)
∇ = CenteredDifference(1, 2,dx_i, inL)
Δ = CenteredDifference(2, 2, dx_i, inL)
c = spzeros(inL)
fA(m, N) = m<N ? -γ(N)*m : 0
fB(m, N) = m<N ? ρ*m*(1-m/N) + m*γ(N)/2 : 0
A = spzeros(inL, inL)
B = spzeros(inL, inL)
A = - 1. * sparse( 1:inL,1:inL, in_x )
B = sparse( 1:inL,1:inL, (m -> 1*m*(1-m/N) + m*1/2 ).(in_x) )
OpB = DiffEqArrayOperator(B)

function mfA(du, u, p, t)
    
end

##
u0_m = zeros(N-1)

# (∇*BC)*DiffEqArrayOperator(A)
# (∇*BC)*DiffEqArrayOperator(A)*u0_m
# (Δ*BC)*DiffEqArrayOperator(B)*u0_m
# # (Δ*BC)*
# ( ((∇*BC)*DiffEqArrayOperator(A)) + ((Δ*BC)*DiffEqArrayOperator(B)) )*u0_m
( ((∇*BC)*DiffEqArrayOperator(A)) + ((Δ*BC)*DiffEqArrayOperator(B)) )
((∇*BC)*DiffEqArrayOperator(A))*u0_m
# ((Δ*BC)*DiffEqArrayOperator(B))*u0_m
# ∇*DiffEqArrayOperator(A)
# ∇
# (∇*DiffEqArrayOperator(A)) + (Δ*DiffEqArrayOperator(B))*BC

# GDA = ((∇*BC)*DiffEqArrayOperator(A))
# GDB = ((Δ*BC)*DiffEqArrayOperator(B))
# ( GDA + GDB )
##
function step!(du, u, (A, B, c), t)
    # nInd = Int(floor(nT(t))) + 1
    # L = GDB*DiffEqArrayOperator(B)
    # L = FPB - FPA

    # c[1] = 1*( 1+1/2+1 ) / dx_i[1]
    # du .= L*u .+ c
    du .= ( ((∇*BC)*DiffEqArrayOperator(A)) + ((Δ*BC)*DiffEqArrayOperator(B)) )*u
    # du .= ((Δ*BC)*DiffEqArrayOperator(B))*u
    # du .+= c
end

t0 = 0.0
t1 = 50
# u0 = zeros(N+1)

prob = ODEProblem(step!, u0_m, (t0, t1), (A, B, c))
# alg = KenCarp4()
alg = TRBDF2() #stablest for stiff PDE
# alg = BS3()
sol = solve(prob, alg, save_everystep=false)

##
