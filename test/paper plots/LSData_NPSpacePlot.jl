using JLD2, LaTeXStrings

##
NOptInterpol_tM_P = Vector{Float64}[]
lVfs_ = Int[]
tM_ = Int[]

let
    @load "data/LSData_NPSpace/LSData_NPSpaceInference_tM3_lVFS500.jld2" _p NOptInterpol_p
    push!(NOptInterpol_tM_P, NOptInterpol_p)
    push!(tM_, 3)
    global _p = _p
end
let
    @load "data/LSData_NPSpace/LSData_NPSpaceInference_tM5_lVFS500.jld2" NOptInterpol_p
    push!(NOptInterpol_tM_P, NOptInterpol_p)
    push!(tM_, 5)
end
let
    @load "data/LSData_NPSpace/LSData_NPSpaceInference_tM8_lVFS500.jld2" NOptInterpol_p
    push!(NOptInterpol_tM_P, NOptInterpol_p)
    push!(tM_, 8)
end
let
    @load "data/LSData_NPSpace/LSData_NPSpaceInference_tM10_lVFS500.jld2" NOptInterpol_p
    push!(NOptInterpol_tM_P, NOptInterpol_p)
    push!(tM_, 10)
end


##
using Plots
pyplot()
theme(:juno, size=(0.9*600,0.9*400), minorgrid=false, gridstyle=:dash, fontfamily="DejaVu Sans",
showaxis=true, gridlinewidth=0.7)
##
palette = cgrad(:tableau_red_blue, 4, categorical=true, rev=true)

fig1 = plot()
for i in 1:length(NOptInterpol_tM_P)-1
    plot!(
        _p, NOptInterpol_tM_P[i], 
        label="tM: "*string(tM_[i]), 
        # marker=:circle,
        color=palette[i],
    )
end
xlabel!(L"p")
ylabel!(L"N")
xlims!(0, 1)
# ylims!(2E4, 4E5)
ylims!(2E4, 3E5)

display(fig1)

# savefig(fig1, "Figures/LSData/NPSpace_mu1.2.pdf")