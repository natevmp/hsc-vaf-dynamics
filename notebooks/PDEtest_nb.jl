# %% codecell
include("src/vafdyn.jl")
using .VAFDyn
# %% codecell
using Plots
plotlyjs()
# gr()
# %% codecell
evolveTime = 59
params = Dict(
    "ρ"=>1.0,
    "μ"=>1.2,
    "ϕ"=>9.0,
    "N"=>2000
)
# %% codecell
dfs = VAFDyn.DFreqspace(params["N"])
@time VAFDyn.evolveVAF(dfs, params, evolveTime, 0.001);
# %% codecell
# l0 = 4001
# l1 = 3001
# l2 = 2501
l3 = 2001
l4 = 1001
l5 = 769
l6 = 500

# vfs0 = VAFDyn.VFreqspace(params["N"],l0)
# @time VAFDyn.evolveVAFfd(vfs0, params, evolveTime);
#
# vfs1 = VAFDyn.VFreqspace(params["N"],l1)
# @time VAFDyn.evolveVAFfd(vfs1, params, evolveTime);
#
# vfs2 = VAFDyn.VFreqspace(params["N"],l2)
# @time VAFDyn.evolveVAFfd(vfs2, params, evolveTime);

vfs3 = VAFDyn.VFreqspace(params["N"],l3)
@time VAFDyn.evolveVAFfd(vfs3, params, evolveTime);

vfs4 = VAFDyn.VFreqspace(params["N"],l4)
@time VAFDyn.evolveVAFfd(vfs4, params, evolveTime);

vfs5 = VAFDyn.VFreqspace(params["N"],l5)
@time VAFDyn.evolveVAFfd(vfs5, params, evolveTime);

vfs6 = VAFDyn.VFreqspace(params["N"],l6)
@time VAFDyn.evolveVAFfd(vfs6, params, evolveTime);
# %% codecell
# using Interpolations

# # linear interpolation
# vfs6Interpol = LinearInterpolation(vfs6.freqs_f, vfs6.n_f);

# difference_f = dfs.n_f[3:end-1] .- vfs6Interpol.(dfs.freqs_f[3:end-1]);
# %% codecell
# plotting

h = plot(dfs.freqs_f[2:end-1], dfs.n_f[2:end-1], label="MC clone size: "*string(params["N"]), yaxis=:log10, linewidth=3)
# plot!(vfs0.freqs_f[2:end-1], vfs0.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l0),
#     linestyle=:dash, linewidth=2)
# plot!(vfs1.freqs_f[2:end-1], vfs1.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l1),
#     linestyle=:dash, linewidth=2)
# plot!(vfs2.freqs_f[2:end-1], vfs2.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l2),
#     linestyle=:dash, linewidth=2)
plot!(vfs3.freqs_f[2:end-1], vfs3.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l3),
    linestyle=:solid, linewidth=2)
plot!(vfs4.freqs_f[2:end-1], vfs4.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l4),
    linestyle=:dashdot, linewidth=2)
plot!(vfs5.freqs_f[2:end-1], vfs5.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l5),
    linestyle=:dash, linewidth=2)
plot!(vfs6.freqs_f[2:end-1], vfs6.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l6),
    linestyle=:dot, linewidth=2)
# xlims!((0, 0.1))
# ylims!((10^-1, 10^4))
xlabel!("variant frequency")
ylabel!("number of variants")
title!("Discrete Freqspace")
display(h)

h2 = plot(dfs.freqs_f[2:end-1], dfs.n_f[2:end-1], label="MC clone size: "*string(params["N"]), yaxis=:log10, linewidth=3)
plot!(vfs3.freqs_f[2:end-1], vfs3.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l3),
    linestyle=:solid, linewidth=2)
plot!(vfs4.freqs_f[2:end-1], vfs4.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l4),
    linestyle=:dashdot, linewidth=2)
plot!(vfs5.freqs_f[2:end-1], vfs5.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l5),
    linestyle=:dash, linewidth=2)
plot!(vfs6.freqs_f[2:end-1], vfs6.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l6),
    linestyle=:dot, linewidth=2)
xlims!((0, 0.2))
ylims!((10^-1, 10^4))
xlabel!("variant frequency")
ylabel!("number of variants")
title!("Discrete Freqspace")
display(h2)

# h3 = plot(dfs.freqs_f[3:end-1], difference_f)
# display(h3)
# %% markdown
# # Large Population
# %% codecell
evolveTimeB = 59
paramsB = Dict(
    "ρ"=>1.0,
    "ϕ"=>9.0,
    "μ"=>1.2,
    "N"=>10000
)
# %% codecell
dfsB = VAFDyn.DFreqspace(paramsB["N"])
@time VAFDyn.evolveVAF(dfsB, paramsB, evolveTimeB, 0.001);
# %% codecell
lB1 = 2001
lB2 = 1001
lB3 = 769
lB4 = 500

vfsB1 = VAFDyn.VFreqspace(paramsB["N"],lB1)
@time VAFDyn.evolveVAFfd(vfsB1, paramsB, evolveTimeB);

vfsB2 = VAFDyn.VFreqspace(paramsB["N"],lB2)
@time VAFDyn.evolveVAFfd(vfsB2, paramsB, evolveTimeB);

vfsB3 = VAFDyn.VFreqspace(paramsB["N"],lB3)
@time VAFDyn.evolveVAFfd(vfsB3, paramsB, evolveTimeB);

vfsB4 = VAFDyn.VFreqspace(paramsB["N"],lB4)
@time VAFDyn.evolveVAFfd(vfsB4, paramsB, evolveTimeB);
# %% codecell
# plotting

hB = plot(dfsB.freqs_f[1:end], dfsB.n_f[1:end], label="MC clone size: "*string(paramsB["N"]), yaxis=:log10, linewidth=3)
plot!(vfsB1.freqs_f[2:end-1], vfsB1.n_f[2:end-1]/paramsB["N"], label="PDE space length: "*string(lB1),
    linestyle=:solid, linewidth=2)
plot!(vfsB2.freqs_f[2:end-1], vfsB2.n_f[2:end-1]/paramsB["N"], label="PDE space length: "*string(lB2),
    linestyle=:dashdot, linewidth=2)
plot!(vfsB3.freqs_f[2:end-1], vfsB3.n_f[2:end-1]/paramsB["N"], label="PDE space length: "*string(lB3),
    linestyle=:dash, linewidth=2)
plot!(vfsB4.freqs_f[2:end-1], vfsB4.n_f[2:end-1]/paramsB["N"], label="PDE space length: "*string(lB4),
    linestyle=:dot, linewidth=2)
xlabel!("variant frequency")
ylabel!("number of variants")
title!("Discrete Freqspace")
display(hB)

hB2 = hB
xlims!((0, 0.2))
ylims!((10^-5, 10^4))
xlabel!("variant frequency")
ylabel!("number of variants")
title!("Discrete Freqspace")
display(hB2)
# %% codecell
# Small population size
evolveTime = 59
params = Dict(
    "μ"=>1.2,
    "ρ"=>1.0,
    "ϕ"=>9.0,
    "N"=>1000
)

# dfs = VAFDyn.DFreqspace(params["N"])
# @time VAFDyn.evolveVAF(dfs, params, evolveTime);

l0 = 151
l1 = 101
l2 = 151
l3 = 201
l4 = 301
l5 = 401
l6 = 500
l7 = 751
l8 = 1001

vfs0 = VAFDyn.VFreqspace(params["N"],l0)
println(vfs0.freqs_f[1:5])
println(vfs0.freqs_f[end-3:end])
@time VAFDyn.evolveVAFfd(vfs0, params, evolveTime);

# vfs1 = VAFDyn.VFreqspace(params["N"],l1)
# println(vfs1.freqs_f[1:5])
# println(vfs1.freqs_f[end-3:end])
# @time VAFDyn.evolveVAFfd(vfs1, params, evolveTime);
#
# vfs2 = VAFDyn.VFreqspace(params["N"],l2)
# println(vfs2.freqs_f[1:5])
# println(vfs2.freqs_f[end-3:end])
# @time VAFDyn.evolveVAFfd(vfs2, params, evolveTime);
#
# vfs3 = VAFDyn.VFreqspace(params["N"],l3)
# println(vfs3.freqs_f[1:5])
# println(vfs3.freqs_f[end-3:end])
# @time VAFDyn.evolveVAFfd(vfs3, params, evolveTime);
#
# vfs4 = VAFDyn.VFreqspace(params["N"],l4)
# println(vfs4.freqs_f[1:5])
# println(vfs4.freqs_f[end-3:end])
# @time VAFDyn.evolveVAFfd(vfs4, params, evolveTime);
#
# vfs5 = VAFDyn.VFreqspace(params["N"],l5)
# println(vfs5.freqs_f[1:5])
# println(vfs5.freqs_f[end-3:end])
# @time VAFDyn.evolveVAFfd(vfs5, params, evolveTime);
#
# vfs6 = VAFDyn.VFreqspace(params["N"],l6)
# println(vfs6.freqs_f[1:5])
# println(vfs6.freqs_f[end-3:end])
# @time VAFDyn.evolveVAFfd(vfs6, params, evolveTime);
#
# vfs7 = VAFDyn.VFreqspace(params["N"],l7)
# println(vfs7.freqs_f[1:5])
# println(vfs7.freqs_f[end-3:end])
# @time VAFDyn.evolveVAFfd(vfs7, params, evolveTime);
#
# vfs8 = VAFDyn.VFreqspace(params["N"],l8)
# println(vfs8.freqs_f[1:5])
# println(vfs8.freqs_f[end-3:end])
# @time VAFDyn.evolveVAFfd(vfs8, params, evolveTime);

# %% codecell

plot(vfs0.freqs_f[2:end]-vfs0.freqs_f[1:end-1])


# %% codecell

# plotting

h = plot(dfs.freqs_f[2:end-1], dfs.n_f[2:end-1], label="MC clone size: "*string(params["N"]), yaxis=:log10, linewidth=4)
plot!(vfs0.freqs_f[2:end-1], vfs0.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l0),
    linestyle=:dash, linewidth=2)
plot!(vfs1.freqs_f[2:end-1], vfs1.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l1),
    linestyle=:dash, linewidth=2)
plot!(vfs2.freqs_f[2:end-1], vfs2.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l2),
    linestyle=:dash, linewidth=2)
plot!(vfs3.freqs_f[2:end-1], vfs3.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l3),
    linestyle=:dash, linewidth=2)
plot!(vfs4.freqs_f[2:end-1], vfs4.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l4),
    linestyle=:dash, linewidth=2)
plot!(vfs5.freqs_f[2:end-1], vfs5.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l5),
    linestyle=:dash, linewidth=2)
plot!(vfs6.freqs_f[2:end-1], vfs6.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l6),
    linestyle=:dash, linewidth=2)
# plot!(vfs7.freqs_f[2:end-1], vfs7.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l7),
#     linestyle=:dash, linewidth=2)
# plot!(vfs8.freqs_f[2:end-1], vfs8.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l8),
#     linestyle=:dash, linewidth=2)
xlabel!("variant frequency")
ylabel!("number of variants")
title!("Discrete Freqspace")
display(h)

h2 = plot(dfs.freqs_f[2:end-1], dfs.n_f[2:end-1], label="MC clone size: "*string(params["N"]), yaxis=:log10, linewidth=3)
plot!(vfs0.freqs_f[2:end-1], vfs0.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l0),
    linestyle=:dash, linewidth=2)
# plot!(vfs1.freqs_f[2:end-1], vfs1.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l1),
#     linestyle=:dash, linewidth=2)
# plot!(vfs2.freqs_f[2:end-1], vfs2.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l2),
#     linestyle=:dash, linewidth=2)
# plot!(vfs3.freqs_f[2:end-1], vfs3.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l3),
#     linestyle=:dash, linewidth=2)
plot!(vfs4.freqs_f[2:end-1], vfs4.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l4),
    linestyle=:dash, linewidth=2)
plot!(vfs5.freqs_f[2:end-1], vfs5.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l5),
    linestyle=:dash, linewidth=2)
plot!(vfs6.freqs_f[2:end-1], vfs6.n_f[2:end-1]/params["N"], label="PDE space length: "*string(l6),
    linestyle=:dash, linewidth=2)
xlims!((0, 0.2))
ylims!((10^-1, 10^4))
xlabel!("variant frequency")
ylabel!("number of variants")
title!("Discrete Freqspace")
display(h2)
