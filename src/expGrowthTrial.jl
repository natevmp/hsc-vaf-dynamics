module expTrial

include("vafSim.jl")
using Random
using Distributions
using SparseArrays
using CSV

function birthDeathGrowthExp(params, tStop, tSaveStep)

	Nin = params["N initial"]
	Nf = params["N final"]
	μ = params["μ"]
	#p = params["p"]
	β = params["β"]
	λ = params["λ"]
	gR = params["growth rate"]
	S = params["sample size"]

	maxMuts = Integer(round(200*μ*Nf))
	muts_loc_cell = falses(maxMuts, Nf)
	# muts_loc_cell = sparse(falses(maxMuts, N))
	mutPrevs_loc = zeros(Int16, maxMuts)

	mLive = 0
	mFixed = 0

	#γ(n) = gR #logisticGrowthRate(n, gR, Nf)
	dt() = rand(  Exponential( 1/( (gR) ) )  )

	tSaves_t = tSaveStep:tSaveStep:tStop
	t_size = size(tSaves_t)[1]

	vaf_n_t = zeros(Int64,Nf+1,t_size)
	vafB_n_t = zeros(Int64,S+1,t_size)
	#burden_m_t = zeros(Int(round(tmax*r*μ*N)),size(trec)[1])

	times_t = Float64[]
	push!(times_t, 0.)
	nLive_t = Int[]
	push!(nLive_t, Nin)

	nLive = Nin
	t = 0
	t += dt()
	# nSymDivs = 0
	# nAsymDivs = 0

	tCounter = 1
	events = 0
	while (t < tStop)&(nLive < Nf)
		events += 1
		#pGrowth = gR/(λ+gR)

		#pSym = λ*(1-p)/(λ+gR)

		#if nLive < Nf
		#pGrowth = 1
		#pSym = 0
		#else
			#pGrowth = 0
			#pSym = 1-p
		#end
		# pASym = λ*p/(λ+γ(nLive))

		#event = rand()
		#if event < pGrowth	# ===== growth =====
		nNow = nLive
		for k = 1:nNow
			# choose individual cell for symmetric division
			divCID = k #rand(1:nNow)
			#cdeath = rand()
			if rand()> β
				nLive += 1
				newCID = nLive
			# add copy of self-renewing cell ID to mutation matrix
				muts_loc_cell[:, newCID] = muts_loc_cell[:, divCID]
				mutPrevs_loc += muts_loc_cell[:, newCID]
			# randomly mutate daughters with on average μ mutations
				mLive += VAFSim.mutateCell!(muts_loc_cell, mutPrevs_loc, mLive, newCID, μ)
			end
			mLive += VAFSim.mutateCell!(muts_loc_cell, mutPrevs_loc, mLive, divCID, μ)

			if nLive >= Nf
				break
			end

		end
		#end

		# clean up the gene by removing all mutations that can't change anymore
		# mLive, mFixed = cleanGenes!(muts_loc_cell, mutPrevs_loc, nLive, mLive, mFixed)

		t += dt()

		if (t>tSaves_t[tCounter])|(nLive >= Nf)
			push!(nLive_t, nLive)
			push!(times_t, t)
			bottleneck_inds = randperm(Nf)[1:S]
			if nLive > S
				vafB_n_t[:,tCounter] = VAFSim.VAFcalc(muts_loc_cell[:, bottleneck_inds], S, mLive)
			end
			vaf_n_t[:,tCounter] = VAFSim.VAFcalc(muts_loc_cell, Nf, mLive)

			tCounter += 1
		end
	end

	#return

	bottleneck_inds = randperm(Nf)[1:S]

	burden_m = VAFSim.burdencalc(muts_loc_cell, nLive, mLive)
	if nLive > S
		burdenB_m = VAFSim.burdencalc(muts_loc_cell[:, bottleneck_inds], S, mLive)
	else
		burdenB_m = 0
	end

	return times_t, nLive_t, vaf_n_t, vafB_n_t, burden_m, burdenB_m, events
end

end
