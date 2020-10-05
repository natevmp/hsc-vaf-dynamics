module VAFSim

using Random
using Distributions
using SparseArrays

export birthDeathShort

export birthDeathAlt


"""
Simulate the population.
N:				size of the population
Nbn:			size of the population at the bottleneck
tmax:			the amount of time the simulation runs
r:				cell division rate
muts_loc_cell:	matrix for which mutations all individuals hold
mLive:			number of currently non-fixated mutations
mFixed:			number of currently fixated mutations
μ:				mutation rate
p:				probability of asymmetric (vs symmetric) division
vafB_n:			VAF at the bottleneck
vaf_n:			VAF for the full population
"""
function birthDeathAlt(N, μ, p, Nbn, tmax, r)
	"""
	maxMuts is the maximum number of mutations to keep track of.
	number is based on experience, but it throws an error if the simulation
	hits the limit just once, so it can't bias the result
	it might be better to use variable length, not quite sure though
	==> Nate: there's a package called ElasticArrays which implements variable length multidimensional arrays.
	"""
	maxMuts = Integer(round(200*μ*N*(1-p/2)/(1-p)))
	muts_loc_cell = falses(maxMuts, N)
	# muts_loc_cell = sparse(falses(maxMuts, N))
	mutPrevs_loc = zeros(Int16, maxMuts)
	mLive = 0
	mFixed = 0

	divExp = Exponential(1/(r*N))
	t = 0
	# dt = randexp()/(r*N)
	dt = rand(divExp)
	t += dt

	while t < tmax

		# choose individual cells for death/birth
		deathCID = rand(1:N)
		# symmetric replication event inside if-loop
		if rand()>p
			birthCID = rand(1:N)
			# kill dying cell and replace with copy of dividing cell
			mutPrevs_loc -= muts_loc_cell[:, deathCID]
			mutPrevs_loc += muts_loc_cell[:, birthCID]
			muts_loc_cell[:, deathCID] = muts_loc_cell[:, birthCID]

			# randomly mutate new individual with on average μ mutations
			if birthCID != deathCID
				nMuts = rand(Poisson(μ))
				if nMuts > 0
					for i = 1:nMuts
						mLive += 1
						muts_loc_cell[mLive, birthCID] = true
						mutPrevs_loc[mLive] += 1
					end
				end
			end
		end
		# asymmetric replication events are equivalent to only mutations happening
		# in a single cell
		nMuts = rand(Poisson(μ))
		if nMuts > 0
			for i = 1:nMuts
				mLive += 1
				muts_loc_cell[mLive, deathCID] = true
				mutPrevs_loc[mLive] += 1
			end
		end


		# clean up the gene by removing all mutations that can't change anymore

		mLive, mFixed = cleanGenes!(muts_loc_cell, mutPrevs_loc, N, mLive, mFixed)

		# dt = randexp()/(r*N)
		dt = rand(divExp)
		t += dt
	end

	# after simulation, record VAF before and after bottleneck
	bottleneck_inds = randperm(N)[1:Nbn]
	#burdenB_m = burdencalc(muts_loc_cell[:, bottleneck_inds], Nbn, mLive)
	#burden_m = burdencalc(muts_loc_cell, N, mLive)

	distanceB_m = createtree(muts_loc_cell[:, bottleneck_inds], Nbn, mLive)
	distance_m = createtree(muts_loc_cell, N, mLive)

	#distanceB_m = distancecalc(muts_loc_cell[:, bottleneck_inds], Nbn, mLive)
	#distance_m = distancecalc(muts_loc_cell, N, mLive)


	vafB_n = VAFcalc(muts_loc_cell[:, bottleneck_inds], Nbn, mLive)
	vaf_n = VAFcalc(muts_loc_cell, N, mLive)

	return vaf_n, vafB_n, mFixed, mLive, distanceB_m,distance_m
end

function birthDeathShort(N, μ, Nbn, steps)
	#=
	maxMuts is the maximum number of mutations to keep track of.
	number is based on experience, but it throws an error if the simulation
	hits the limit just once, so it can't bias the result
	it might be better to use variable length, not quite sure though
	==> Nate: there's a package called ElasticArrays which implements variable length multidimensional arrays.
	=#
	maxMuts = Integer(round(10*μ*N))
	muts_loc_cell = falses(maxMuts, N)
	# muts_loc_cell = sparse(falses(maxMuts, N))
	mutPrevs_loc = zeros(Int16, maxMuts)
	mLive = 0
	mFixed = 0

	for step = 1:steps
		# choose individual cells for death/birth
		birthCID = rand(1:N)
		deathCID = rand(1:N)
		# kill dying cell and replace with copy of dividing cell
		mutPrevs_loc -= muts_loc_cell[:, deathCID]
		mutPrevs_loc += muts_loc_cell[:, birthCID]
		muts_loc_cell[:, deathCID] = muts_loc_cell[:, birthCID]

		# randomly mutate new individual with on average μ mutations
		nMuts = rand(Poisson(μ))
		if nMuts > 0
			for i = 1:nMuts
				mLive += 1
				muts_loc_cell[mLive, deathCID] = true
				mutPrevs_loc[mLive] += 1
			end
		end
		# clean up the gene by removing all mutations that can't change anymore
		mLive, mFixed = cleanGenes!(muts_loc_cell, mutPrevs_loc, N, mLive, mFixed)
	end

	# record VAF before and after bottleneck
	bottleneck_inds = randperm(N)[1:Nbn]

	vafB_n = VAFcalc(muts_loc_cell[:, bottleneck_inds], Nbn, mLive)
	vaf_n = VAFcalc(muts_loc_cell, N, mLive)

	# single cell mutational burden
	nMuts_cell = mutBurdenSingleCell(muts_loc_cell)

	return vaf_n, vafB_n, nMuts_cell, mFixed, mLive
end

"""
Simulate a logistically growing population.
N:				maximum size of the population
Nbn:			size of the population at the bottleneck
tmax:			the amount of time the simulation runs
r:				cell division rate
d:				death rate
muts_loc_cell:	matrix for which mutations all individuals hold
mLive:			number of currently non-fixated mutations
mFixed:			number of currently fixated mutations
indLive:		individuals currently alive
μ:				mutation rate
p:				probability of asymmetric (vs symmetric) division
vafB_n:			VAF at the bottleneck
vaf_n:			VAF for the full population
"""
function birthDeathLog(N, indStart, μ, p, Nbn, trec, r, d)
	tmax = trec[end]
	maxMuts = Integer(round(200*μ*N*(1-p/2)/(1-p)))
	muts_loc_cell = falses(maxMuts, N)
	# muts_loc_cell = sparse(falses(maxMuts, N))
	mutPrevs_loc = zeros(Int16, maxMuts)
	mLive = 0
	mFixed = 0


	indLive = indStart
	divExp = Exponential(1/(r*indLive))
	t = 0
	# dt = randexp()/(r*N)
	dt = rand(divExp)
	t += dt
	k = 1

	vafB_n_t = zeros(Int64,Nbn+1,size(trec)[1])
	vaf_n_t = zeros(Int64,N+1,size(trec)[1])

	burdenB_m_t = zeros(Int(round(tmax*r*μ*N)),size(trec)[1])

	burden_m_t = zeros(Int(round(tmax*r*μ*N)),size(trec)[1])

	indLive_t = zeros(size(trec)[1])

	mFixed_t = zeros(size(trec)[1])

	depthmax_t = zeros(size(trec)[1])
	depthmaxB_t = zeros(size(trec)[1])
	depthmean_t = zeros(size(trec)[1])
	depthmeanB_t = zeros(size(trec)[1])

	while t < tmax
		#println(t)
		#println(indLive)
		q = (1-indLive/N)./d
		h = 0
		# choose individual cells for death/birth
		deathCID = rand(1:indLive)
		# symmetric replication event inside if-loop
		if rand()<q
			birthCID = rand(1:indLive)
			indLive += 1
			deathCID = indLive
			h = 1
		elseif rand()>p
			birthCID = deathCID
			while birthCID==deathCID
				birthCID = rand(1:indLive)
			end
			h = 1
		end

		if h==1
			# kill dying cell and replace with copy of dividing cell
			mutPrevs_loc -= muts_loc_cell[:, deathCID]
			mutPrevs_loc += muts_loc_cell[:, birthCID]
			muts_loc_cell[:, deathCID] = muts_loc_cell[:, birthCID]

			# randomly mutate new individual with on average μ mutations

			nMuts = rand(Poisson(μ))
			if nMuts > 0
				for i = 1:nMuts
					mLive += 1
					muts_loc_cell[mLive, birthCID] = true
					mutPrevs_loc[mLive] += 1
				end
			end

		end
		# asymmetric replication events are equivalent to only mutations happening
		# in a single cell
		nMuts = rand(Poisson(μ))
		if nMuts > 0
			for i = 1:nMuts
				mLive += 1
				muts_loc_cell[mLive, deathCID] = true
				mutPrevs_loc[mLive] += 1
			end
		end


		# clean up the gene by removing all mutations that can't change anymore

		mLive, mFixed = cleanGenes!(muts_loc_cell, mutPrevs_loc, N, mLive, mFixed)

		# dt = randexp()/(r*N)
		divExp = Exponential(1/(r*indLive))
		dt = rand(divExp)
		t += dt

		if t > trec[k]
			bottleneck_inds = randperm(indLive)[1:Nbn]
			#println(VAFcalc(muts_loc_cell[:, bottleneck_inds], Nbn, mLive))
			vafB_n_t[:,k] = VAFcalc(muts_loc_cell[:, bottleneck_inds], Nbn, mLive)
			vaf_n_t[1:indLive+1,k] = VAFcalc(muts_loc_cell, indLive, mLive)
			burdenB_m_t[1+mFixed:mLive+1+mFixed,k] = burdencalc(muts_loc_cell[:, bottleneck_inds], Nbn, mLive)
			burden_m_t[1+mFixed:mLive+1+mFixed,k] = burdencalc(muts_loc_cell, indLive, mLive)
			indLive_t[k] = indLive
			mFixed_t[k] = mFixed

			inc_parent_off = VAFSim.inclusion(muts_loc_cell,N,mLive)
			incB_parent_off = VAFSim.inclusion(muts_loc_cell[:, bottleneck_inds],Nbn,mLive)

			depthmax_t[k] = maximum(sum(inc_parent_off,dims=1))
			depthmaxB_t[k] = maximum(sum(incB_parent_off,dims=1))
			depthmean_t[k] = mean(sum(inc_parent_off,dims=1))
			depthmeanB_t[k] = mean(sum(incB_parent_off,dims=1))

			k += 1
		end
	end

	return vaf_n_t, vafB_n_t, mFixed_t, mLive, indLive_t , burden_m_t, burdenB_m_t, depthmax_t, depthmaxB_t, depthmean_t, depthmeanB_t #, muts_loc_cell #, distanceB_m,distance_m
end

function birthDeathLogC(N, indStart, μ, C, Nbn, trec, r, d)
	p = C/(N*r)
	tmax = trec[end]
	maxMuts = Integer(round(200*μ*N*(1-p/2)/(1-p)))
	muts_loc_cell = falses(maxMuts, N)
	# muts_loc_cell = sparse(falses(maxMuts, N))
	mutPrevs_loc = zeros(Int16, maxMuts)
	mLive = 0
	mFixed = 0


	indLive = indStart
	divExp = Exponential(1/(r*indLive))
	t = 0
	# dt = randexp()/(r*N)
	dt = rand(divExp)
	t += dt
	k = 1

	vafB_n_t = zeros(Int64,Nbn+1,size(trec)[1])
	vaf_n_t = zeros(Int64,N+1,size(trec)[1])
	indLive_t = zeros(size(trec)[1])

	while t < tmax
		#println(t)
		#println(indLive)
		q = (1-indLive/N)./d
		h = 0
		p = C/(indLive*r)
		# choose individual cells for death/birth
		deathCID = rand(1:indLive)
		# symmetric replication event inside if-loop

		if rand()<p

		elseif rand()<q
			birthCID = rand(1:indLive)
			indLive += 1
			deathCID = indLive
			h = 1
		elseif rand()>p
			birthCID = deathCID
			while birthCID==deathCID
				birthCID = rand(1:indLive)
			end
			h = 1

		end

		if h==1
			# kill dying cell and replace with copy of dividing cell
			mutPrevs_loc -= muts_loc_cell[:, deathCID]
			mutPrevs_loc += muts_loc_cell[:, birthCID]
			muts_loc_cell[:, deathCID] = muts_loc_cell[:, birthCID]

			# randomly mutate new individual with on average μ mutations

			nMuts = rand(Poisson(μ))
			if nMuts > 0
				for i = 1:nMuts
					mLive += 1
					muts_loc_cell[mLive, birthCID] = true
					mutPrevs_loc[mLive] += 1
				end
			end

		end
		# asymmetric replication events are equivalent to only mutations happening
		# in a single cell
		nMuts = rand(Poisson(μ))
		if nMuts > 0
			for i = 1:nMuts
				mLive += 1
				muts_loc_cell[mLive, deathCID] = true
				mutPrevs_loc[mLive] += 1
			end
		end


		# clean up the gene by removing all mutations that can't change anymore

		mLive, mFixed = cleanGenes!(muts_loc_cell, mutPrevs_loc, N, mLive, mFixed)

		# dt = randexp()/(r*N)
		divExp = Exponential(1/(r*indLive))
		dt = rand(divExp)
		t += dt

		if t > trec[k]
			bottleneck_inds = randperm(indLive)[1:Nbn]
			#println(VAFcalc(muts_loc_cell[:, bottleneck_inds], Nbn, mLive))
			vafB_n_t[:,k] = VAFcalc(muts_loc_cell[:, bottleneck_inds], Nbn, mLive)
			vaf_n_t[1:indLive+1,k] = VAFcalc(muts_loc_cell, indLive, mLive)
			indLive_t[k] = indLive
			k += 1
		end
	end



	# after simulation, record VAF before and after bottleneck
	#bottleneck_inds = randperm(indLive)[1:Nbn]
	#burdenB_m = burdencalc(muts_loc_cell[:, bottleneck_inds], Nbn, mLive)
	#burden_m = burdencalc(muts_loc_cell, indLive, mLive)

	#distanceB_m = createtree(muts_loc_cell[:, bottleneck_inds], Nbn, mLive)
	#distance_m = createtree(muts_loc_cell, N, mLive)

	#distanceB_m = distancecalc(muts_loc_cell[:, bottleneck_inds], Nbn, mLive)
	#distance_m = distancecalc(muts_loc_cell, N, mLive)


	#vafB_n = VAFcalc(muts_loc_cell[:, bottleneck_inds], Nbn, mLive)
	#vaf_n = VAFcalc(muts_loc_cell, indLive, mLive)

	return vaf_n_t, vafB_n_t, mFixed, mLive, indLive_t #, burden_m, burdenB_m,muts_loc_cell #, distanceB_m,distance_m
end

function birthDeathLogLin(N, indStart, μ, p, lin, Nbn, trec, r, d)

	tmax = trec[end]
	maxMuts = Integer(round(200*μ*N*(1-p/2)/(1-p)))
	muts_loc_cell = falses(maxMuts, N*5)
	# muts_loc_cell = sparse(falses(maxMuts, N))
	mutPrevs_loc = zeros(Int16, maxMuts)
	mLive = 0
	mFixed = 0


	indLive = indStart
	divExp = Exponential(1/(r*indLive))
	t = 0
	# dt = randexp()/(r*N)
	dt = rand(divExp)
	t += dt
	k = 1

	vafB_n_t = zeros(Int64,Nbn+1,size(trec)[1])
	vaf_n_t = zeros(Int64,N*5,size(trec)[1])

	burdenB_m_t = zeros(Int(round(tmax*r*μ*N)),size(trec)[1])

    burden_m_t = zeros(Int(round(tmax*r*μ*N)),size(trec)[1])

	indLive_t = zeros(size(trec)[1])

	mFixed_t = zeros(size(trec)[1])


	while t < tmax
		#println(t)
		#println(indLive)
		q = (1-indLive/N)./d
		if q < lin
			q = lin
		end

		h = 0
		# choose individual cells for death/birth
		deathCID = rand(1:indLive)
		# symmetric replication event inside if-loop

		if rand()<p

		elseif rand()<q
			birthCID = rand(1:indLive)
			indLive += 1
			deathCID = indLive
			h = 1
		elseif rand()>p
			birthCID = deathCID
			while birthCID==deathCID
				birthCID = rand(1:indLive)
			end
			h = 1

		end

		if h==1
			# kill dying cell and replace with copy of dividing cell
			mutPrevs_loc -= muts_loc_cell[:, deathCID]
			mutPrevs_loc += muts_loc_cell[:, birthCID]
			muts_loc_cell[:, deathCID] = muts_loc_cell[:, birthCID]

			# randomly mutate new individual with on average μ mutations

			nMuts = rand(Poisson(μ))
			if nMuts > 0
				for i = 1:nMuts
					mLive += 1
					muts_loc_cell[mLive, birthCID] = true
					mutPrevs_loc[mLive] += 1
				end
			end

		end
		# asymmetric replication events are equivalent to only mutations happening
		# in a single cell
		nMuts = rand(Poisson(μ))
		if nMuts > 0
			for i = 1:nMuts
				mLive += 1
				muts_loc_cell[mLive, deathCID] = true
				mutPrevs_loc[mLive] += 1
			end
		end


		# clean up the gene by removing all mutations that can't change anymore

		mLive, mFixed = cleanGenes!(muts_loc_cell, mutPrevs_loc, N, mLive, mFixed)

		# dt = randexp()/(r*N)
		rc = N/indLive
		divExp = Exponential(1/(rc*indLive))
		dt = rand(divExp)
		t += dt

		if t > trec[k]
			bottleneck_inds = randperm(indLive)[1:Nbn]
			#println(VAFcalc(muts_loc_cell[:, bottleneck_inds], Nbn, mLive))
			vafB_n_t[:,k] = VAFcalc(muts_loc_cell[:, bottleneck_inds], Nbn, mLive)
			vaf_n_t[1:indLive+1,k] = VAFcalc(muts_loc_cell, indLive, mLive)
			burdenB_m_t[1+mFixed:mLive+1+mFixed,k] = burdencalc(muts_loc_cell[:, bottleneck_inds], Nbn, mLive)
			burden_m_t[1+mFixed:mLive+1+mFixed,k] = burdencalc(muts_loc_cell, indLive, mLive)
			indLive_t[k] = indLive
			mFixed_t[k] = mFixed

			k += 1
		end
	end



	# after simulation, record VAF before and after bottleneck
	#bottleneck_inds = randperm(indLive)[1:Nbn]
	#burdenB_m = burdencalc(muts_loc_cell[:, bottleneck_inds], Nbn, mLive)
	#burden_m = burdencalc(muts_loc_cell, indLive, mLive)

	#distanceB_m = createtree(muts_loc_cell[:, bottleneck_inds], Nbn, mLive)
	#distance_m = createtree(muts_loc_cell, N, mLive)

	#vafB_n = VAFcalc(muts_loc_cell[:, bottleneck_inds], Nbn, mLive)
	#vaf_n = VAFcalc(muts_loc_cell, indLive, mLive)

	return vaf_n_t, vafB_n_t, mFixed, mLive, indLive_t , burden_m_t, burdenB_m_t #, muts_loc_cell #, distanceB_m,distance_m
end



"""
Remove mutations that already died out or fixated.
N:				the number of individuals
muts_loc_cell:	the matrix for which mutations all individuals hold
mLive:			the number of currently live (non-fixated) mutations
mFixed:			the number of currently fixated mutations
"""
function cleanGenes!(muts_loc_cell::AbstractArray{Bool, 2}, mutPrevs_loc, N, mLive, mFixed)
	locOut = mLive
	# for i = 0:locOut-1
	for loc in locOut:-1:1
		# find and remove fixated genes
		if mutPrevs_loc[loc] == N
			# swap fixed and live var mutations
			muts_loc_cell[loc, :] = muts_loc_cell[mLive, :]
			muts_loc_cell[mLive, :] .= false
			mutPrevs_loc[loc] = mutPrevs_loc[mLive]
			mutPrevs_loc[mLive] = 0
			# update final live mutation index
			mLive -= 1
			mFixed += 1
		end
		# find and remove genes that died out

		if (mutPrevs_loc[loc] == 0)&&(mLive>0)
			# swap dead and live var mutations
			muts_loc_cell[loc, :] = muts_loc_cell[mLive, :]
			muts_loc_cell[mLive, :] .= false
			mutPrevs_loc[loc] = mutPrevs_loc[mLive]
			mutPrevs_loc[mLive] = 0
			# update final live mutation index
			mLive -= 1
		end
	end

	return mLive, mFixed
end

"""
Calculate the VAF distribution.
N:				number of individuals
muts_loc_cell:	matrix for which mutations all individuals hold
mLive:			number of currently non-fixated mutations
vaf_n:			VAF distribution
"""
function VAFcalc(muts_loc_cell::AbstractArray{Bool, 2}, N, mLive)
	# mutPrev_loc = sum(muts_loc_cell, dims=2)
	vaf_n = zeros(Int64, N+1)
	# sort the genes into their respective frequencies
	for i = 1:mLive
		# vaf_n[1+mutPrev_loc[i]] += 1
		vaf_n[1 + sum(muts_loc_cell[i, :])] += 1
	end
	return vaf_n
end


function burdencalc(muts_loc_cell::AbstractArray{Bool, 2}, N, mLive)
	# mutPrev_loc = sum(muts_loc_cell, dims=2)
	burden_m = zeros(Int64, mLive+1)
	# sort the genes into their respective frequencies
	for i = 1:N
		# vaf_n[1+mutPrev_loc[i]] += 1
		burden_m[1 + sum(muts_loc_cell[:, i])] += 1
	end
	return burden_m
end

function distancecalc(muts_loc_cell::AbstractArray{Bool, 2}, N, mLive)
	# mutPrev_loc = sum(muts_loc_cell, dims=2)
	distance_m = zeros(Int64, mLive+1)
	# sort the genes into their respective frequencies
	for i = 1:N
		for j = (i+1):N
			# vaf_n[1+mutPrev_loc[i]] += 1
				distance_m[1 + sum(abs.(muts_loc_cell[:, i]-muts_loc_cell[:,j]))] += 1

		end
	end
	return distance_m
end

function inclusion(muts_loc_cell::AbstractArray{Bool, 2}, N, mLive)

	inc_parent_off = zeros(mLive,mLive)

	for k = 1:mLive
		for l = 1:mLive
			if k != l
				if sum(muts_loc_cell[k,:])>0
					if minimum(muts_loc_cell[k,:]-muts_loc_cell[l,:]) == 0
						inc_parent_off[k,l] = 1
					end
				end
			end
		end
	end
	return inc_parent_off
end

function createtree(muts_loc_cell::AbstractArray{Bool, 2}, N, mLive)

	distance_m = zeros(Int64, mLive+1)

	mutscurrent_loc_cell = deepcopy(muts_loc_cell)

	h = 1
	cdis = 0
	cind = N

	maxloops = 100000
	loop = 1
	#println(muts_loc_cell)
	while cind > 1
		loop += 1
		nocdis = 1
		#println(string("loop = ", loop, "; cdis = ", cdis, "; cind = ", cind))
		#println(muts_loc_cell)
		for k = cind:-1:1
			#println(string("k =", k))
			for l = (k-1):-1:1
				#println(string("k =", k, "; l =", l, "; distance = ", sum(abs.(mutsc_loc_cell[:, k]-mutsc_loc_cell[:,l])) ))
				if k != l
					#println(string(" l =", l, "; distance = ", sum(abs.(mutscurrent_loc_cell[:, k]-mutscurrent_loc_cell[:,l])) ))
					if sum(abs.(mutscurrent_loc_cell[:, k]-mutscurrent_loc_cell[:,l])) == cdis
						#println(muts_loc_cell)
						newmuttemp = ((mutscurrent_loc_cell[:, k] + mutscurrent_loc_cell[:, l]).>1)
						distance_m[1 + sum(abs.(mutscurrent_loc_cell[:, k]-newmuttemp))] += 1
						distance_m[1 + sum(abs.(mutscurrent_loc_cell[:, l]-newmuttemp))] += 1
						mutscurrent_loc_cell[:,k] =  mutscurrent_loc_cell[:,cind]
						mutscurrent_loc_cell[:,l] = newmuttemp
						mutscurrent_loc_cell[:,cind] .= 0
						cind -= 1
						nocdis = 0
						cdis = 0
						#println(muts_loc_cell)
						break
					end
				end
			end
		end
		#println(muts_loc_cell)
		if nocdis == 1
			cdis += 1
		end
		if loop > maxloops
			println("createtree stopped early")
			break
		end
	end
	return distance_m
end

"""
Get mutational burden of each single cell in the population.
"""
function mutBurdenSingleCell(muts_loc_cell::AbstractArray{Bool, 2})
	nMuts_cell = vec(sum(muts_loc_cell, dims=1))
end




end
