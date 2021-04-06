include("../src/vafdyn.jl")
using .VAFDyn
using Plots
plotly()

vfs = VAFDyn.VFreqspace(100000, 500)

println(length(vfs.freqs_f))
println(length(vfs.n_f))
println(length(vfs))
println(vfs.freqs_f[2])

h1 = scatter( vfs.freqs_f, ones(length(vfs)), ylims=(0,2))
# xlims!((0, 1))
# ylims!((0, 2))
display(h1)