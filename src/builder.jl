
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

"""
Create discrete probability space (DFreqspace object) from irregularly spaced probability density object (VFreqspace) through spline interpolation and scaling
"""
function makeDFSfromVFS(vfs::VFreqspace, N::Integer; normalize=false)
    if !normalize
        cnSpline = Spline1D(vfs.freqs_f[2:end-1], vfs.n_f[2:end-1])
        dfs = DFreqspace(N)
        dfs.n_f[2:end-1] .= cnSpline(dfs.freqs_f[2:end-1])/N
    else
        Δf_f = ( (vfs.freqs_f[3:end]-vfs.freqs_f[2:end-1]) .+ (vfs.freqs_f[2:end-1]-vfs.freqs_f[1:end-2]) )/2
        cnSpline = Spline1D(vfs.freqs_f[2:end-1], N * vfs.n_f[2:end-1] .* Δf_f)
        dfs = DFreqspace(N)
        dfs.n_f[2:end-1] .= cnSpline(dfs.freqs_f[2:end-1])
    end
    return dfs
end

