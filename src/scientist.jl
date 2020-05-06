

function mutBurden(dfs::DFreqspace)::Float64
    sum(dfs.n_f[2:end-1])
end

function mutBurden(cfs::CFreqspace)::Float64
    sum(cfs.n_f[2:end-1]) * cfs.df
end

function sampler(trueFs::DFreqspace, n::Integer)
    N = trueFs.popsize
    sampFs = DFreqspace(n)
    for m in 1:n
        sampFs.n_f[m] = 
            sum( 
                [ trueFs.n_f[v] * pdf(Hypergeometric(v, N-v, n), m)
                for v=1:N ]
            )
    end
    trueMuts = mutBurden(trueFs)
    sampMuts = 
        trueMuts - 
        sum( 
            [ trueFs.n_f[v]*
            ( pdf(Hypergeometric(v, N-v, n), n) 
            + pdf(Hypergeometric(v, N-v, n), 0) ) 
            for v=1:N ]
        )

    return sampFs, sampMuts
end


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

