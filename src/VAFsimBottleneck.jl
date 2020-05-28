using Random
using Distributions

## function to remove genes that already died out or fixated and keep track of
##how many mutations are still 'in transit' and how many fixated
# n is the number of individuals
# genes is the matrix for which mutations all individuals hold
# k is the number of currently non-fixated mutations
# c is the number of currently fixated mutations
function clean_genes(n,genes::Array{Int8,2},k,c)
	genes_sum = sum(genes,dims=2)
	jump_size = 0
	i = 1
	k_o = k
	for i = 0:k_o-1
		# find fixated genes
		if genes_sum[k_o-i] == n
			c = c + 1
			genes[k_o-i,1:n] = genes[k,1:n]
			genes[k,1:n] = zeros(Int8,1,n)
			k = k - 1
			jump_size = jump_size + 1
		end
		# find genes that died out
		if genes_sum[k_o-i] == 0
			genes[k_o-i,1:n] = genes[k,1:n]
			genes[k,1:n] = zeros(Int8,1,n)
			k = k - 1
		end
		i=i+1
	end

	return k,c
end

## function to calculate the VAF distribution
# n is the number of individuals
# genes is the matrix for which mutations all individuals hold
# k is the number of currently non-fixated mutations
# VAF is the VAF distribution
function VAFcalc(genes::Array{Int8,2},n,k)
	gene_sum = sum(genes,dims=2)
	VAF = zeros(Int64,n+1)
	# sort the genes into their respective frequencies
	for i = 1:k
		VAF[gene_sum[i]+1] = VAF[gene_sum[i]+1]+1
	end
	return VAF
end


## function for the actual simulation of the population
# n is the number of individuals
# steps is the number of reproductive events
# genes is the matrix for which mutations all individuals hold
# k is the number of currently non-fixated mutations
# c is the number of currently fixated mutations
# bsize is the size of the population at the bottleneck
# mu is the mutation rate
# VAFb is the VAF at the bottleneck
# VAF is the VAF for the full population
function BirthdeathShort(n,steps,genes::Array{Int8,2},k,c,bsize,mu,VAFb::Array{Float64,1},VAF::Array{Float64,1})
	# loop over all steps
	for step = 1:steps

		# clean up the gene by removing all mutations that can't change anymore
		c_last = c
		k,c = clean_genes(n,genes,k,c)

		# choose individuals for death/birth
		birth = rand(1:n)
		death = rand(1:n)
		genes[:,death] = genes[:,birth]

		# randomly mutate new individual with on average mu mutations
		mutstep = rand(Poisson(mu))

		if mutstep > 0
			for i = 1:mutstep
				k=k+1

				genes[k,death] = 1
			end
		end


	end

	# after simulation, record VAF before and after bottleneck

	bottleneck = randperm(n)[1:bsize]
	VAFb= VAFcalc(genes[:,bottleneck],bsize,k)

	VAF = VAFcalc(genes,n,k)


	return k,c, VAFb, VAF
end


## Initializer for the main function, supposedly makes everything run slightly faster
## n is the number of individuals
## steps is the number of reproductive events
## genes is the matrix for which mutations all individuals hold
## k is the number of currently non-fixated mutations
## c is the number of currently fixated mutations
## bsize is the size of the population at the bottleneck
## mu is the mutation rate
## VAFb is the VAF at the bottleneck
##VAF is the VAF for the full population

function RunBDShort(steps,n,mu,bsize)

	# initialize all values

	# max genes is the maximum number of mutations to keep track of.
	# number is based on experience, but it throws an error if the simulation
	# hits the limit just once, so it can't bias the result
	# it might be better to use variable length, not quite sure though
	max_genes = 15*mu*n
	genes = zeros(Int8,max_genes,n)
	k = 0
	c = 0
	VAFb = zeros(Float64,bsize+1)
	VAF = zeros(Float64,n+1)
	# run actual simulation
	k,c,VAFb,VAF = BirthdeathShort(n,steps,genes,k,c,bsize,mu,VAFb,VAF)

	return genes,k,c,VAFb,VAF
end
