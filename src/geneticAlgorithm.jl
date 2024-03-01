include("helperFunctions.jl")

function deriveNextGeneration(Xk::Matrix, f::Function, pDict::Dict;
    verbose::Bool=false)

    n, p = size(Xk)
    Xkp1 = myfill(Xk, -75.0)
    Fkp1 = myzeros(Fk)
    survived = zeros(p)
    popSize = 0
    keepAddingToPopulation = true
    fevals = 0
    actions = Dict(
        :bothParentsSurvived => 0, :onlyOneParentSurvived => 0, :noParentSurvived => 0, :crossover => 0, :genFitnessImproved => 0,
        :genFitnessNotImproved => 0, :mutation => 0,
        :mutationFailure => 0, :mutationSuccess => 0, :OneChild => 0, :OneChildOneParent => 0,
        :OneChildBothParents => 0, :parentsSelected => 0)


    # add the top survivor, it is assumed that Xk, Fk have the best values at the front
    popSize += 1
    Xkp1[:, popSize] = Xk[:, 1]
    Fkp1[popSize] = Fk[1]

    # TODO increment suitable actions for the top survivor

    while keepAddingToPopulation

        if popSize >= p
            keepAddingToPopulation = false
            myprintln(verbose, "Next Generation Complete.")
        end

        p1, p2 = selectParents(Xk, Fk)
        actions[:parentsSelected] += 1

        o = crossover(p1, p2) # offspring
        actions[:crossover] += 1

        om, mutations_1mut = mutate(o, f, pDict)
        fevals += 1
        actions[:mutations] += mutations_1mut

        Xkp1, Fkp1, popSize, survived, actions_1dAA = decideAndAdd(p1, p2, om, Xkp1, Fkp1, popSize, pDict) # select the child, and if choosing to let the parent survive by default, the best min(2, p-popSize) parents
        actions = merge(+, actions, actions_1dAA)

    end

    Xkp1, Fkp1 = fittestFirst(Xkp1, Fkp1)

    if Fkp1[1] > Fk[1]
        actions[:genFitnessImproved] += 1
        myprintln(verbose, "Fittest individual now even fitter.")
    else
        actions[:genFitnessNotImproved] += 1
        myprintln(verbose, "Fittest individual has same fitness as previous generation.")
    end

    return Xkp1, Fkp1, fevals, actions
end

"""
    createInitialPopulation(x0, popSize, f, pDict; deviationFactor = 0.1, verbose::Bool = false)

Generates an initial population for a genetic algorithm, with individuals arrayed around a central point and scaled by a deviation factor. The population is then evaluated for fitness, and the fittest individual is positioned first in the population array.

# Arguments
- `x0`: A vector representing the central point in the parameter space around which the initial population is generated.
- `popSize`: The size of the population, i.e., the number of individuals to generate.
- `f`: The fitness function to be minimized. It must accept an individual's parameter vector and the parameter dictionary `pDict`, and it should return a scalar fitness value.
- `pDict`: A dictionary containing additional parameters required by the fitness function.
- `deviationFactor`: A factor that controls the spread of the initial population around the central point. Defaults to 0.1.
- `verbose`: If `true`, enables the printing of function-related messages.

The function uses a Halton sequence to create diverse starting points for the population, ensuring a broad and non-collapsing initial search space. The deviation factor is used to fine-tune the spread of the population around the initial central point `x0`.

The fitness of each individual in the population is evaluated, and these values are used to determine the fittest individual. The population and fitness arrays are then adjusted to position the fittest individual at the start.

# Notes
- The Halton sequence is a low-discrepancy, quasi-random sequence used to generate sample points that are more uniformly distributed than with simple random sampling.
- The fitness function evaluations will increment the total function evaluations count by the population size (`popSize`) after the call to this function.
- The function assumes the fitness function is to be minimized and inverts the fitness values so that higher values represent fitter individuals.

# Returns
- `X0`: An `n x p` matrix representing the initial population, where `n` is the number of parameters for each individual, and `p` is the population size.
- `F0`: A vector of inverted fitness values corresponding to `X0`, with the fittest individual's fitness at the first index.

# Examples
```julia
# Generate the initial population
X0, F0 = createInitialPopulation(x0, popSize, your_fitness_function, pDict, deviationFactor=0.1)
```
"""
function createInitialPopulation(x0, popSize, f, pDict; # increment fevals by p after the call
    deviationFactor = 0.1,
    verbose::Bool = false)

    n = length(x0)
    p = popSize
    # Generate n x p matrix using Halton sequence
    haltonMatrix = sampleSpaceHalton(n, p) 

    # Adjust the matrix to fit the deviation factor and center around x0
    X0 = zeros(n, p)
    for i in 1:n
        for j = 1:p
            # Scale and shift the Halton sequence points
            X0[i, j] = x0[i] * (1 - deviationFactor) + (x0[i] * 2 * deviationFactor) * haltonMatrix[i, j]
        end
    end

    X0[:, end] = x0 
    fvals = zeros(p)

    for i ∈ 1:p
        fvals[i] = f(X0[:, i], pDict, getGradientToo=false)
    end

    F0 = ones(p) .+ maximum(fvals) .- fvals # F0, X0 are still corresponding

    # Find the fittest individual (max in F0)
    fittest_index = argmax(F0)

    # Swap the fittest individual to the front if it's not already there
    if fittest_index != 1
        X0[:, [1, fittest_index]] = X0[:, [fittest_index, 1]]
        F0[1], F0[fittest_index] = F0[fittest_index], F0[1]
    end

    return X0, F0

end

# x0 = pr.x0
# f = pr.objective
# n = length(x0)
# p = n+1
# X0, F0 = createInitialPopulation(x0, p, f, pr.p)
# function selectParents(Xk, Fk,;
#     verbose::Bool = false)

#     pdf = Fk./sum(Fk)
#     cdf = cumsum(pdf)
# end