include("../src/vafSim.jl")

using .VAFSim
using Statistics, Plots

gr()
theme(:juno)


params = Dict(
    "N"=>500,
    "t"=>10,
    "ρ"=>2.0,
    "ϕ"=>8.0,
    "μ"=>1.2,
    "Nb"=>80
)

# single simulation
function singleSim(params)
    N = params["N"]
    t = params["t"]
    ρ = params["ρ"]
    ϕ = params["ϕ"]
    μ = params["μ"]
    NSample = params["Nb"]
    tsteps = Integer(round(N*t*(ρ+ϕ)))
    # println("number of evolveSteps: "*string(steps))
    vaf_n, vafB_n, nMuts_cell, mFixed, mLive = VAFSim.birthDeathShort(N, μ, NSample, tsteps)
    return vaf_n, vafB_n
end

# @time vaf_n, vafB_n, mFixed, mLive = singleSim(params)
# @time vaf_n, vafB_n, mFixed, mLive = singleSim(params)
# @time vaf_n, vafB_n, mFixed, mLive = singleSim(params)

function averageSims(params, nSims)
    N = params["N"]
    t = params["t"]
    ρ = params["ρ"]
    ϕ = params["ϕ"]
    μ = params["μ"]
    Nb = params["Nb"]

    vafs_sim_n = Array{Float64}(undef, nSims, N+1)
    vafBs_sim_n = Array{Float64}(undef, nSims, Nb+1)

    println("running sims:")
    for i in 1:nSims
        vaf_n, vafB_n = singleSim(params)
        vafs_sim_n[i, :] = vaf_n
        vafBs_sim_n[i, :] = vafB_n
    end

    vafAv_n = vec(mean(vafs_sim_n, dims=1))
    vafBAv_n = vec(mean(vafBs_sim_n, dims=1))

    vafStd_n = vec(std(vafs_sim_n, mean=vafAv_n, dims=1))
    vafBStd_n = vec(std(vafBs_sim_n, mean=vafBAv_n, dims=1))

    return vafAv_n, vafBAv_n, vafStd_n, vafBStd_n
end

# println("number of mutations fixated: ", mFixed)
# println("number of mutations in pop: ", mLive)
# println(mLive == sum(vaf_n))
# println(vaf_n)

@time vafAv_n, vafBAv_n, vafStd_n, vafBStd_n = averageSims(params, 100)

# println(size(vafStd_n))
# println(vafStd_n)
freqsN_f = 0:1/params["N"]:1
freqsNB_f = 0:1/params["Nb"]:1

# h1 = plot(freqsN_f[2:end], vafAv_n[2:end], yaxis=:log10, ylims=(10^-3,10^3), xlims=(0,0.2), label="true", color=1)
# plot!(freqsN_f[2:end], vafAv_n[2:end] + vafStd_n[2:end], color=1, linestyle=:dash, label="")
# plot!(freqsN_f[2:end], vafAv_n[2:end] - vafStd_n[2:end], color=1, linestyle=:dash, label="")
# plot!(freqsNB_f[2:end], vafBAv_n[2:end], label="sampled", color=2)
# plot!(freqsNB_f[2:end], vafBAv_n[2:end] + vafBStd_n[2:end], color=2, linestyle=:dash, label="")
# plot!(freqsNB_f[2:end], vafBAv_n[2:end] - vafBStd_n[2:end], color=2, linestyle=:dash, label="")
# xlabel!("variant frequency")
# ylabel!("number of variants")
# title!("t = "*string(params["t"]))
# display(h1)
#
# h2 = plot(freqsN_f[2:end], vafAv_n[2:end], ylims=(10^-3,10^3), xlims=(0,0.2), label="true", color=1)
# plot!(freqsN_f[2:end], vafAv_n[2:end] + vafStd_n[2:end], color=1, linestyle=:dash, label="")
# plot!(freqsN_f[2:end], vafAv_n[2:end] - vafStd_n[2:end], color=1, linestyle=:dash, label="")
# plot!(freqsNB_f[2:end], vafBAv_n[2:end], label="sampled", color=2)
# plot!(freqsNB_f[2:end], vafBAv_n[2:end] + vafBStd_n[2:end], color=2, linestyle=:dash, label="")
# plot!(freqsNB_f[2:end], vafBAv_n[2:end] - vafBStd_n[2:end], color=2, linestyle=:dash, label="")
# xlabel!("variant frequency")
# ylabel!("number of variants")
# title!("t = "*string(params["t"]))
# display(h2)
#
# h3 = plot(freqsN_f[2:end], vafStd_n[2:end], label="true", xlims=(0,0.2))
# plot!(freqsNB_f[2:end], vafBStd_n[2:end], label="sampled")
# xlabel!("variant frequency")
# ylabel!("VAF standard deviation")
# title!("t = "*string(params["t"]))
# display(h3)
