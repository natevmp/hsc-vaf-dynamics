

function mutBurden(dfs::DFreqspace)::Float64
    sum(dfs.n_f[2:end-1])
end

function mutBurden(cfs::CFreqspace)::Float64
    sum(cfs.n_f[2:end-1]) * cfs.df
end

function sampler(trueFs::DFreqspace, s::Integer)
    n = trueFs.popsize
    sampFs = DFreqspace(s)
    for u in 1:s
        sampFs.n_f[1+u] = 
            sum( 
                [ trueFs.n_f[1+v] * pdf(Hypergeometric(v, n-v, s), u)
                for v=1:n ]
            )
    end
    return sampFs
end

function samplerMuts(trueFs::DFreqspace, s::Integer)
    n = trueFs.popsize
    sampMuts = 
        sum(
            [ sum(
                [ trueFs.n_f[1+v] * pdf(Hypergeometric(v, n-v, s), u)
                for v=1:n ]
                )
            for u in 1:s-1 ]
            )
    return sampMuts
end

# pHypgeom(u,n,v,s) = stat.hypergeom.pmf(u, n, v, s)

# function samplerPy(trueFs::DFreqspace, s::Integer)
#     # n = true size
#     # s = sample size
#     # v = true pop variant iterator
#     # u = sample pop variant iterator
#     n = trueFs.popsize
#     sampFs = DFreqspace(s)
#     for u in 1:s
#         sampFs.n_f[1+u] = 
#             sum(
#                 [ trueFs.n_f[1+v] * pHypgeom(u, n, v, s) 
#                 for v=1:n ]
#             )
#     end
#     return sampFs
# end

# ===== Run experiments =====


function vafMC(params::Dict, time::Real, dt::Real)
    dfs = DFreqspace(params["N"])
    evolveVAF(dfs, params, time, dt)
    return dfs
end

function vafDiff(params::Dict, time::Real, dt::Real, l::Int)
    cfs = CFreqspace(l)
    println("running PDE")
    evolveVAF(cfs, params, time, dt)
    return cfs
end

