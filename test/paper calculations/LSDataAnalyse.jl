using DelimitedFiles, JLD2, Statistics

## ===== Load Data =====

untypedM = readdlm("data/Shearwater_calls_FDR0.95_all_muts.txt", '\t', Any; skipstart=1)
untypedM = untypedM[:, 5: end-1]
replace!(untypedM, "NA"=>0)

variants_var_col = Array{Int}(untypedM)

## ===== Order data =====
HSCMask_col = fill(false, size(variants_var_col, 2))
HSCMask_col[1:73] .= true
# HSCMask_col[125:end] .= true
HPCMask_col = .!HSCMask_col

sampleSize = sum(HSCMask_col)


## ===== variant allele frequencies =====
prevHSC_var = sum(variants_var_col[:, HSCMask_col], dims=2)

function makeHist(varPrev_var, N::Integer)
    n_f = zeros(Int64, N+1)
    for m in 1:N
        n_f[1+m] = length(findall(x -> x==m, varPrev_var))
    end
    return n_f
end

nVHSC_f = makeHist(prevHSC_var, sampleSize)
freqs_f = range(0, 1; length=sampleSize+1)


##

SCBurdenHSC_CID = vec(sum(variants_var_col[:, HSCMask_col], dims=1))
SCBurdenHSCMean = mean(SCBurdenHSC_CID)
SCBurdenHSCVar = var(SCBurdenHSC_CID)

## save data

# :SCBurdenHSC_CID
# :SCBurdenHSCMean
# :SCBurdenHSCVar
# :freqs_f
# :nVHSC_f
# :sampleSize
filename = "./data/LSDataStats.jld2"
filename = "./data/LSDataStatsBM.jld2"
# filename = "./data/LSDataStatsTot.jld2"
@save filename SCBurdenHSC_CID SCBurdenHSCMean SCBurdenHSCVar freqs_f nVHSC_f sampleSize


