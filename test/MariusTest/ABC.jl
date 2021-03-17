include("vafdyn.jl")


using StatsBase
using .VAFDyn
using Plots
using JLD2
using PyCall
stats = pyimport("scipy.stats")

##function to find the best fitting time for a given VAF distribution
# VAF is the given VAF distribution that we want to fit

# esize is the number of 'employed bees' of the algorithm

# osize is the number of 'onlooker bees'

# cnumber is the cycle number, which is the maximum number of cycles we search

# limit is how often a 'bee' searches around the same spot without improving
# until it gives up

# maxsteps tells us the maximum number of steps we want to suppose, so that we
# can search in [0,maxsteps]

function ABC(VAF,esize,osize,cnumber,limit,maxsteps)

    n = length(VAF)-1
    # time solution is the array of solutions that each 'bee' is searching at
    time_solution = zeros(esize)
    # fitness is the array of fitness values of each solution
    fitness = zeros(esize)
    # repeats is the array of times each bee has searched at its current solution
    # without improvement
    repeats = zeros(esize)

    #randomly initialize all solutions
    time_solution = rand(esize).*maxsteps
    counter = 0

    h = 0

    VAF_normalized = VAF ./ VAF[2]

    while h == 0

        ## fitness calculation
        # we loop over the population of bees and everyone looks around its chosen
        # spot to find an improvement
        for k = 1:esize

            tempsteps = time_solution[k] - maxsteps/20 + rand()*maxsteps/10

            params_fit = Dict("ρ" => 1, "μ" => 1, "N" => n, "ϕ" => 0)

            dfs_fit = VAFDyn.DFreqspace(params_fit["N"])

            VAFDyn.evolveVAF(dfs_fit, params_fit, tempsteps/n, 1/(100*n))

            fit_normalized = dfs_fit.n_f ./ dfs_fit.n_f[2]

            fitnesstemp = 1/(1+sum((fit_normalized[2:end] - VAF_normalized[2:end]) .^ 2))


            if fitnesstemp > fitness[k]

                time_solution[k] = tempsteps
                fitness[k] = fitnesstemp
                repeats[k] = 0

            else
                repeats[k] += 1
            end

        end

        ##onlooker phase
        # onlooker bees choose good food sources and try to improve the local solution further
        for l=1:osize

            t_new = StatsBase.sample(1:esize,Weights(fitness))
            tempsteps = time_solution[t_new] - maxsteps/20 + rand()*maxsteps/10

            params_fit = Dict("ρ" => 1, "μ" => 1, "N" => n, "ϕ" => 0)

            dfs_fit = VAFDyn.DFreqspace(params_fit["N"])

            VAFDyn.evolveVAF(dfs_fit, params_fit, tempsteps/n, 1/(100*n))

            fit_normalized = dfs_fit.n_f ./ dfs_fit.n_f[2]

            fitnesstemp = 1/(1+sum((fit_normalized[2:end] - VAF_normalized[2:end]) .^ 2))


            if fitnesstemp > fitness[t_new]

                time_solution[t_new] = tempsteps
                fitness[t_new] = fitnesstemp
            end

        end

        ## scout phase
        # we again loop over all 'bees', and, if it has repeatedly searched at its
        # spot without finding an improvement, it instead changes location
        for l=1:esize

            if repeats[l] >= limit
                if fitness[l] < maximum(fitness)
                    time_solution[l] = rand().*maxsteps


                end
            end
        end
        if counter >= cnumber
            h = 1
        end

        counter += 1
    end

    tempsteps = time_solution[argmax(fitness)]

    params_fit = Dict("ρ" => 1, "μ" => 1, "N" => n, "ϕ" => 0)

    dfs_fit = VAFDyn.DFreqspace(params_fit["N"])

    VAFDyn.evolveVAF(dfs_fit, params_fit, tempsteps/n, 1/(10*n))

    #estimate mu based on the number of mutations at frequency 1
    mu_est = VAF[2] / dfs_fit.n_f[2]

    C_solution = time_solution[argmax(fitness)]  / n

    vafSol = dfs_fit.n_f


    return time_solution[argmax(fitness)],mu_est,vafSol
end
