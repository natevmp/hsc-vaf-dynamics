using PyCall
stats = pyimport("scipy.stats")

##function to find the VAF after a bottleneck based on hypergeometric drawing
# VAF is the VAF distribution before the bottleneck
# n is the population before the bottleneck
# bsize is the population size at the bottleneck

function sample(VAF,n,bsize)
    VAF_hyper = zeros(bsize)

    #loop over all original frequencies
    for l = 1:n
        h = stats.hypergeom(n,bsize,l)
        #loop over the bottleneck frequencies and calculate how much of the original
        # distribution converts to there
        for k = 1:bsize

            VAF_hyper[k] = VAF_hyper[k] + (VAF[l+1]*h.pmf(k))
        end
    end
    return VAF_hyper
end
