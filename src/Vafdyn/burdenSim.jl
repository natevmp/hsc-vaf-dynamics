module BurdenSim

using Random
using Distributions
# using SparseArrays

function cappedExponentialGrowth(Ni, K, r, t)
	return Ni*exp(r*t)<K ? Ni*exp(r*t) : K
end

function evolveSCBurden(params::Dict, evolveTime::Number, tSaveStep::Number)
    Ni = params["N initial"]
	Nf = params["N final"]
	μ = params["μ"]
	p = params["p"]
	λ = params["λ"]
	gR = params["growth rate"]
	S = params["sample size"]
    ρ = λ*(1-p)
    ϕ = λ*p

    mutDist = Poisson(μ)
    dt(n) = rand( Exponential(1/(λ*n)) )
    nTime(tt) = cappedExponentialGrowth(Ni, Nf, gR, tt)
    
    scBurden_cid = zeros(Int64, Ni)
    times_t = Float64[]
    tSaves_t = tSaveStep:tSaveStep:evolveTime
    
    t = 0
    nLive = Ni
    push!(times_t, t)
	nLive_t = Int64[]
	push!(nLive_t, Ni)
    tCounter = 1

    while t<evolveTime
        
        while nLive < round(nTime(t))
            # grow pop
            divCID = rand(1:nLive)
            push!(scBurden_cid, scBurden_cid[divCID])
            nLive += 1
            scBurden_cid[divCID] += rand(mutDist)
            scBurden_cid[end] += rand(mutDist)
        end
        
        t += dt(nLive)
        event = rand()
        if event > p
            # symmetric divisions
            srCID = rand(1:nLive)
            diffCID = rand(1:nLive)
            scBurden_cid[diffCID] = scBurden_cid[srCID]
            scBurden_cid[srCID] += rand(mutDist)
            if diffCID != srCID
                scBurden_cid[diffCID] += rand(mutDist)
            end
        else
            # asymmetric division
            divCID = rand(1:nLive)
            scBurden_cid[divCID] += rand(mutDist)
        end

        # save time-resolved data
		if t>tSaves_t[tCounter]
			push!(nLive_t, nLive)
			push!(times_t, t)
			tCounter += 1
		end
    end

    return scBurden_cid, times_t, nLive_t
end

function burdenHist(scBurden_cid::AbstractVector{Int})
    mMax = maximum(scBurden_cid)
	nCells_m = zeros(Int64, mMax+1)
	# sort the genes into their respective frequencies
	for m in scBurden_cid
		nCells_m[1 + m] += 1
	end
	return nCells_m
end

function sampler(scBurden_cid::AbstractVector{Int}, S::Int)
    # bottleneck_inds = randperm(nLive)[1:S]
    if S > length(scBurden_cid)
        println("warning: sampe size must be smaller than population size")
    end
    return shuffle(scBurden_cid)[1:S]
end


end