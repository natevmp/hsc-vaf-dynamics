using CSV, DataFrames
using Gadfly, Cairo, Fontconfig
using Statistics
##

dataIn1 = CSV.File(
    "./data/blood769869-sup-table1.csv", 
    normalizenames=true,
    header=2,
    datarow=3
) |> DataFrame


##
dataIn2 = CSV.File(
    "./data/blood769869-sup-table2.csv", 
    normalizenames=true,
    header=2,
    datarow=3
) |> DataFrame
dataZink = DataFrame(gene=dataIn2.Gene_Symbol, vaf=dataIn2.VAF, age=dataIn2.Age_at_blood_draw)

## ===== plot all genes ===== 
fig1 = plot(
    dataZink, 
    x=:age, y=:vaf,
    color=:gene,
    Geom.point
)
display(fig1)
# ===== plot selected genes ===== 
genetag_r = falses(size(dataZink)[1])
for gene in ["TP53", "DNMT3A", "TET2", "JAK2"]
    genetag_r = genetag_r .| (dataZink[!, :gene].==gene)
end

fig2 = plot(
    dataZink[genetag_r, :], 
    x=:age, y=:vaf,
    color=:gene,
    Geom.point
)
display(fig2)

## ====================== Functions for getting averages ======================

function binAgeData!(data::DataFrame; ageBins=10)
    lData = size(data, 1)

    edges_ = range(0, 100; length=ageBins+1)
    centers_ = edges_[1:end-1] .+ (edges_[2]-edges_[1])/2

    data.ageBin = zeros(lData)
    for i in 1:lData
        indAge = findmin(abs.( data.age[i] .- centers_ ))[2]
        ageBin = centers_[indAge]
        data.ageBin[i] = ageBin
    end
    return data
end

function averageVafAge(data::DataFrame; ageBins=10)
    
    edges_ = range(0, 100; length=ageBins+1)
    _age = edges_[1:end-1] .+ (edges_[2]-edges_[1])/2
    
    vafAv_age = Vector{Union{Float64, Missing}}(undef, length(_age))
    vafVar_age = Vector{Union{Float64, Missing}}(undef, length(_age))
    vafAv_age[:] .= missing
    vafVar_age[:] .= missing
    for (i,age) in enumerate(_age)
        inds_ = (data[:, :ageBin] .== age)
        if isempty(data.vaf[inds_])
            vafAv_age[i] = missing
            vafVar_age[i] = missing
        else
            vafAv_age[i] = mean(data.vaf[inds_])
            if length(data.vaf[inds_])<2
                vafVar_age[i] = missing
            else
                vafVar_age[i] = var(data.vaf[inds_])
            end
        end
    end
    return _age, vafAv_age, vafVar_age
end



## ======= Plotting =======
ageBins = 25


binAgeData!(dataZink, ageBins=ageBins)
_age, vafAv_age, vafVar_age = averageVafAge(dataZink; ageBins=ageBins)
dataPoolAv = DataFrame(age=_age, all_genes=vafAv_age)
dataPoolVar = DataFrame(age=_age, all_genes=vafVar_age)


genes_ = ["TP53", "DNMT3A", "TET2", "JAK2", "ASXL1", "PPM1D"]
genetagAll_r = falses(size(dataZink,1))
for gene in genes_
    genetag_r = (dataZink[!, :gene].==gene)
    genetagAll_r .= genetagAll_r .| (dataZink[!, :gene].==gene)
    _age, vafAv_age, vafVar_age = averageVafAge(dataZink[genetag_r, :], ageBins=ageBins)
    dataPoolAv[:, Symbol(gene)] = vafAv_age
end
##
subfig1 = layer(
    dataPoolAv,
    x=:age, y=:all_genes,
    Geom.line
)

subfig2 = layer(
    dataZink, 
    x=:age, y=:vaf,
    color=:gene,
    Geom.point
)


fig3 = plot(
    subfig1, subfig2,
    Guide.xlabel("age"),
    Guide.ylabel("VAF"),
    style(point_size=2pt, highlight_width=0.5pt)
)
display(fig3)

fig3 |> PDF("Figures/ZinkData/vafAges_allGenes.pdf")

##

genes_ = ["TP53", "DNMT3A", "TET2", "JAK2"]
genetag_r = falses(size(dataZink)[1])
for gene in genes_
    genetag_r = genetag_r .| (dataZink[!, :gene].==gene)
end

plotgrid = Array{Plot}(undef, (2,2))
for (i, gene) in enumerate(genes_)
    fig = plot(
        dataZink[dataZink[!, :gene].==gene, :],
        x=:age,
        y=:vaf,
        size=[2pt],
        color=:gene
    )
    plotgrid[i] = fig
end
fig4 = gridstack(plotgrid, )
display(fig4)

fig4 |> PDF("Figures/ZinkData/vafAges_selectedGenes.pdf")

##





