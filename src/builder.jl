
function makeCFreqspace(l::Integer)

    # infinitesimal step in frequency space
    df = 1/(l-1)
    freqs_f = Array(range(0, 1; length=l))
    nVAF_f = zeros(Float64, l)
    return (nVAF_f, freqs_f, df)
end

function addClones(dfs::DFreqspace, size::Integer, n::Real=1.)
    dfs.n_f[size+1] += n
end

function addClones(cfs::CFreqspace, freq::Real, n::Real=1.)
    ind = Integer(round(freq / cfs.df ))
    cfs.n_f[1+ind] += n * 1/cfs.df
end