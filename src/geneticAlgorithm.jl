include("helperFunctions.jl")

function deriveNextGeneration(Xk,
    Fk,
    fk, 
    f::Function, 
    pDict::Dict;
    delta = 0.1,
    dftol = 1e-15,
    deviation = 0.1,
    Dist = randn,
    parentsSurvive = true,
    verbose::Bool=false)

    n, p = size(Xk)
    # @show p
    Xkp1 = myfill(Xk, -75.0)
    Fkp1 = myzeros(Fk)
    survived = zeros(p)
    popAdded = 0
    keepAddingToPopulation = true
    fevals = 0
    actions = Dict(
        :bothParentsSurvived => 0, 
        :crossover => 0, 
        :genFitnessImproved => 0,
        :genFitnessNotImproved => 0,
        :mutation => 0,
        :noParentSurvived => 0,  
        :onlyOneParentSurvived => 0,
        :parentsSelected => 0)


    # add the top survivor, it is assumed that Xk, Fk have the best values at the front
    popAdded += 1
    Xkp1[:, popAdded] = Xk[:, 1]
    Fkp1[popAdded] = Fk[1]
    survived[popAdded] = 1

    # TODO increment suitable actions for the top survivor

    while keepAddingToPopulation

        if popAdded >= p
            keepAddingToPopulation = false
            myprintln(verbose, "Next Generation Complete.")
            break
        end

        p1, p2 = selectParents(Xk, Fk, verbose=verbose)
        fevals += p
        actions[:parentsSelected] += 1

        o = crossover(p1, p2) # offspring
        actions[:crossover] += 1

        om, mutations_1mut = mutation(o, f, pDict, 
        delta=delta, deviation=deviation, Dist=Dist, verbose=verbose)
        fevals += 1
        actions[:mutation] += mutations_1mut

        # myprintln(verbose, "popAdded = $popAdded before decideAndAdd()")

        # select the child, and if choosing to let the parent survive by default, the best min(2, p-popAdded) parents
        Xkp1, Fkp1, popAdded, survived, actions_1dAA = 
            decideAndAdd(p1, p2, om, Xkp1, Fkp1, popAdded, survived, parentsSurvive=parentsSurvive, verbose=false)

        actions = merge(+, actions, actions_1dAA)

        # myprintln(verbose, "popAdded = $popAdded after decideAndAdd()")

    end

    Xkp1, Fkp1, fkp1 = reevaluateFitness(Xkp1, Fkp1, f, pDict, verbose=verbose)
    fevals += p

    Xkp1, Fkp1 = fittestFirst(Xkp1, Fkp1, verbose=verbose)

    if fkp1 < fk
        actions[:genFitnessImproved] += 1
        myprintln(verbose, "Fittest individual now even fitter.")
    elseif fkp1 == fk
        actions[:genFitnessNotImproved] += 1
        myprintln(verbose, "Fittest individual has same  fitness as previous generation.")
    elseif fkp1 - fk < dftol 
        actions[:genFitnessNotImproved] += 1
        myprintln(verbose, "Fittest individual has pretty much same fitness as previous generation.")
    else
        @error "How come fitness has decreased? Need to investigate."
    end

    return Xkp1, Fkp1, fkp1, fevals, actions
end

function reevaluateFitness(Xkp1, Fkp1, f::Function, pDict::Dict;
    verbose::Bool = false)

    n, p = size(Xkp1)
    fvals = zeros(p)

    for i ∈ 1:p
        fvals[i] = f(Xkp1[:, i], pDict, getGradientToo=false)
    end

    fkp1 = minimum(fvals)
    Fkp1 = ones(p) .+ maximum(fvals) .- fvals # F0, X0 are still corresponding
    Xkp1, Fkp1 = fittestFirst(Xkp1, Fkp1, verbose=verbose)
    myprintln(verbose, "New Generation created with fittest point $(Xkp1[:, 1]) having fitness of $(Fkp1[1]).")

    return Xkp1, Fkp1, fkp1
end

"""
    createInitialPopulation(x0, popSize, f, pDict; deviationFactor=0.1, verbose=false)

Generates an initial population for a genetic algorithm, positioning the fittest individual at the forefront of the population array. This setup provides a starting point for the evolutionary process, with diversity introduced via the Halton sequence and controlled spread around a central point `x0`.

# Arguments
- `x0`: Central point in the parameter space, serving as a reference for generating the initial population.
- `popSize`: Size of the population to be generated.
- `f`: Fitness function used to evaluate each individual. It should accept a parameter vector and a dictionary of additional parameters `pDict`, returning a scalar fitness value.
- `pDict`: Dictionary containing additional parameters required by the fitness function.
- `deviationFactor`: Controls the variance of the initial population from the central point `x0`. Defaults to 0.1.
- `verbose`: Enables detailed logging of the function's operations if set to `true`.

# Returns
- `X0`: A matrix of the initial population's parameter vectors, with dimensions `n x p`, where `n` is the number of parameters, and `p` is the population size. The fittest individual is placed in the first column.
- `F0`: A vector of the inverted fitness values for the initial population, aligned with `X0`, positioning the highest fitness value first.
- `f0`: The minimum fitness value among the initial population, indicating the fittest individual's performance.

# Notes
- Utilizes the Halton sequence for generating a diverse initial population around `x0`, ensuring a broad exploration of the parameter space.
- Fitness values are calculated inversely to the objective function values to align with a maximization strategy: higher fitness values indicate better individuals.
- The population and fitness arrays are adjusted so that the fittest individual based on `F0` is positioned first, facilitating direct access for further genetic operations.

# Examples
```julia
# Define the central point and population size
x0 = [1.0, 2.0]
popSize = 50

# Define the fitness function and additional parameters
your_fitness_function(x, pDict) = sum(x.^2)  # Example fitness function
pDict = Dict()

# Generate the initial population
X0, F0, f0 = createInitialPopulation(x0, popSize, your_fitness_function, pDict, deviationFactor=0.1, verbose=true)
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

    f0 = minimum(fvals)
    F0 = ones(p) .+ maximum(fvals) .- fvals # F0, X0 are still corresponding
    X0, F0 = fittestFirst(X0, F0, verbose=verbose)
    myprintln(verbose, "Initial Generation created with fittest point $(X0[:, 1]) having fitness of $(F0[1]).")


    return X0, F0, f0

end

"""
    fittestFirst(Xk::Matrix, Fk::Vector; verbose=false)

Reorders the population matrix `Xk` and the corresponding fitness vector `Fk` to position the fittest individual first. This adjustment is crucial for genetic algorithms and similar evolutionary methods, providing immediate access to the current optimal solution.

# Arguments
- `Xk`: A matrix where each column represents the parameter vector of an individual in the population.
- `Fk`: A vector containing the fitness values for each individual in `Xk`, with higher values indicating better fitness in the context of maximization.
- `verbose`: If set to `true`, enables the output of detailed logging information during the operation.

# Returns
- `Xk`: The updated population matrix with the fittest individual's parameter vector in the first column.
- `Fk`: The updated fitness vector with the highest fitness value moved to the first position.

# Method
1. Determine the index of the individual with the highest fitness value in `Fk`.
2. If this individual is not already in the first position, perform a swap operation to move their data (both parameter vector and fitness value) to the front.
3. Log detailed information if `verbose` is enabled, providing insights into the reordering process and the fitness of the relocated individual.

# Usage Example
```julia
# Given a population matrix `Xk` and corresponding fitness values `Fk`
Xk = [1.0 2.0; 3.0 4.0]  # Example population
Fk = [0.5, 0.7]          # Fitness values

# Reorder the population to highlight the fittest individual
Xk, Fk = fittestFirst(Xk, Fk, verbose=true)

# The output will show the fittest individual now at the first index, with logs if verbose is true
```
"""
function fittestFirst(Xk::Matrix, Fk::Vector;
    verbose::Bool = false)
    # Identify the index of the fittest individual
    fittest_index = argmax(Fk)
    F_fittest = Fk[fittest_index]
    myprintln(verbose, "Fittest Individual is at Index $(fittest_index) with Fitness of $(F_fittest).")
    # Swap to bring the fittest individual to the front if necessary
    if fittest_index != 1
        myprintln(verbose, "Since the index is not 1, we swap it to become at index 1.")
        Xk[:, [1, fittest_index]] = Xk[:, [fittest_index, 1]]
        Fk[1], Fk[fittest_index] = Fk[fittest_index], Fk[1]
    end

    return Xk, Fk
end

"""
    selectParents(Xk, Fk; verbose=false)

Selects two unique parents from a population for crossover in a genetic algorithm, ensuring that individuals with higher fitness have a greater chance of being selected. The selection process is probability-based, with each individual's chance of selection proportional to its fitness.

# Arguments
- `Xk`: A matrix representing the current population, where each column corresponds to an individual's parameter vector.
- `Fk`: A vector containing the fitness values for each individual in `Xk`. Higher values indicate better fitness, assuming a maximization problem.
- `verbose`: If set to `true`, prints detailed information about the selection process.

# Returns
- `p1`, `p2`: Tuples representing the selected parents. Each tuple contains:
    - `idx`: The index of the individual in the population.
    - `x`: The parameter vector of the selected individual.
    - `F`: The fitness value of the selected individual.

# Method
1. Normalize the fitness values in `Fk` to create a probability density function (pdf), ensuring the probability of selecting each individual is proportional to its fitness.
2. Utilize the `selectTwoUniqueIndices` function to determine two unique indices based on the pdf, promoting diversity in the selection of parents.
3. Retrieve the parameter vectors (`x`) and fitness values (`F`) for the selected indices from `Xk` and `Fk`, respectively.
4. Organize the selected individuals' information into tuples (`p1`, `p2`) for subsequent use in crossover and mutation operations.

# Examples
```julia
# Given a population matrix `Xk` and corresponding fitness values `Fk`
Xk = [1.0 2.0; 3.0 4.0]  # Example population matrix
Fk = [0.5, 0.7]          # Corresponding fitness values

# Execute parent selection from the population
parent1, parent2 = selectParents(Xk, Fk)

# Output the details of the selected parents
println("Parent 1: ", parent1)
println("Parent 2: ", parent2)
```
""" 
function selectParents(Xk, Fk;
    verbose::Bool = false)

    pdf = Fk./sum(Fk)

    idx1, idx2 = selectTwoUniqueIndices(pdf)

    x1, x2 = Xk[:, idx1], Xk[:, idx2]
    F1, F2 = Fk[idx1], Fk[idx2]
    
    p1 = (idx=idx1, x=x1, F=F1)
    p2 = (idx=idx2, x=x2, F=F2)

    myprintln(verbose, "Choosing parents at indices $idx1 and $idx2 with Fitness values of $F1 and $F2.")
    return p1, p2

end

"""
    selectTwoUniqueIndices(pdf::Vector; verbose::Bool=false)

Selects two unique indices from a population based on the probability density function (pdf) representing selection probabilities for each individual.

# Arguments
- `pdf`: A vector representing the probability density function for the selection probabilities of individuals within a population.
- `verbose`: If `true`, prints detailed information about the selection process.

# Returns
- A tuple containing two unique indices selected based on the pdf. These indices correspond to individuals within the population that have been chosen for further operations (e.g., crossover).

# Method
1. Convert the input pdf to a cumulative distribution function (cdf) for the purpose of selection.
2. Randomly select the first index by generating a random value and identifying its position within the cdf.
3. Create a modified pdf (`pdfm1`) by setting the probability of the already selected index to zero, effectively removing it from consideration for the second selection. Renormalize `pdfm1` to ensure it sums to 1.
4. Convert the updated `pdfm1` to a cdf and select the second index, ensuring it is different from the first.
5. If `verbose` is `true`, print the indices of the selected parents.

# Examples
```julia
pdf = [0.1, 0.2, 0.3, 0.4] # Example pdf representing selection probabilities
idx1, idx2 = selectTwoUniqueIndices(pdf, verbose=true)

# Output may show something like:
# "Choosing parent 1 at index 3."
# "Choosing parent 2 at index 1."
```
"""
function selectTwoUniqueIndices(pdf::Vector;
    verbose::Bool = false)
    # Convert pdf to cdf for the first selection
    # Select the first index
    cdf = cumsum(pdf)
    cdf[end] = 1.0
    rand_val1 = rand()

    idx1 = findfirst(≥(rand_val1), cdf)
    myprintln(verbose, "Choosing parent 1 at index $(idx1).")

    # Set the probability of the selected index to zero and renormalize pdf
    pdfm1 = deepcopy(pdf)
    pdfm1[idx1] = 0
    pdfm1 /= sum(pdfm1)

    # Convert the updated pdf to cdf for the second selection
    cdf = cumsum(pdfm1)
    cdf[end] = 1.0

    # Select the second index
    rand_val2 = rand()
    idx2 = findfirst(≥(rand_val2), cdf)
    myprintln(verbose, "Choosing parent 2 at index $(idx2) from the remaining parents.")
    return idx1, idx2

end

"""
    crossover(p1, p2; verbose=false)

Performs crossover between two parent individuals to produce an offspring. The offspring inherits attributes through a weighted average of the parents' attributes, with the weighting factor determined by the relative fitness of the parents.

# Arguments
- `p1`: A tuple representing the first parent, containing the parent's index (`idx`), parameter vector (`x`), and fitness value (`F`).
- `p2`: A tuple representing the second parent, structured similarly to `p1`.
- `verbose`: If set to `true`, prints detailed information about the crossover process.

# Returns
- `o`: A tuple representing the offspring, which includes an index (`idx`), parameter vector (`x`), and fitness value (`F`). The index is set to `-1` to denote that the offspring is not part of the current generation, and the fitness value is also `-1`, indicating that it has not yet been evaluated.

# Method
1. Extracts the parameter vectors (`x`) and fitness values (`F`) from both parents.
2. Calculates the weighting factor `theta` as the ratio of the first parent's fitness to the sum of both parents' fitness values. This determines the influence of each parent on the offspring's attributes.
3. Computes the offspring's parameter vector (`x0`) as a weighted average of the parents' parameter vectors, applying `theta`.
4. Initializes the offspring's `idx` and `F` to `-1`, marking it as a new, unevaluated individual distinct from its parents.

# Usage Example
```julia
# Assume two parent tuples, each with an index, parameter vector, and fitness value
p1 = (idx=1, x=[1.0, 2.0], F=5.0)
p2 = (idx=2, x=[2.0, 3.0], F=3.0)

# Generate an offspring through crossover
offspring = crossover(p1, p2, verbose=true)

# The offspring's parameter vector is a weighted blend of the parents' vectors
```
"""
function crossover(p1, p2;
    verbose::Bool = false)
    
    @unpack x, F = p1
    x1, F1 = x, F
    @unpack x, F = p2
    x2, F2 = x, F

    theta = F1/(F1+F2)

    x0 = theta*x1 + (1-theta)*x2

    o = (idx=-1, x=x0, F=-1)

    return o

end

"""
    mutation(o, f, pDict; delta=0.1, deviation=0.1, Dist=randn, verbose=false)

Applies mutations to the parameter vector of an offspring `o` with a certain probability and evaluates the mutated offspring's fitness.

# Arguments
- `o`: The offspring to be mutated, represented as a tuple containing its index (`idx`), parameter vector (`x`), and fitness value (`F`).
- `f`: The fitness function used to evaluate the mutated offspring. It must accept an individual's parameter vector and a dictionary of additional parameters, returning a scalar fitness value.
- `pDict`: A dictionary of parameters required by the fitness function.
- `delta`: The probability of mutation for each element in the offspring's parameter vector.
- `deviation`: The magnitude of mutation, influencing how much an individual parameter can change during mutation.
- `Dist`: A distribution function used to generate random values for the mutation. The default is the standard normal distribution (`randn`).
- `verbose`: If set to `true`, enables printing of function-related messages.

# Returns
- `om`: The mutated offspring, structured similarly to the input `o` but with the potentially altered parameter vector and newly evaluated fitness value.
- `mutations`: The total number of mutations applied to the offspring.

# Method
1. Iterates over each element in the offspring's parameter vector.
2. With probability `delta`, mutates the element by multiplying it by `(1 + deviation * Dist())`, where `Dist()` generates a random value from the specified distribution.
3. After potentially mutating each element, evaluates the fitness of the mutated parameter vector using the provided fitness function `f` and parameters `pDict`.
4. Constructs a new tuple for the mutated offspring with its index set to `-1` (indicating it is a new generation) and its fitness value updated based on the mutation.

# Usage Example
```julia
# Define an offspring and the fitness function
offspring = (idx=1, x=[1.0, 2.0], F=5.0)
function your_fitness_function(x, pDict)
    # Fitness calculation...
end

# Perform mutation on the offspring
mutated_offspring, mutation_count = mutation(offspring, your_fitness_function, Dict(), delta=0.2)

# mutated_offspring contains the potentially mutated parameter vector and its new fitness value
```
"""
function mutation(o, f, pDict;
    delta = 0.1,
    deviation = 0.1,
    Dist = randn,
    verbose::Bool = false)

    mutations = 0
    xo = o.x
    xom = deepcopy(xo)

    for i ∈ eachindex(xom)
        if rand() < delta
            xom[i] *= (1 + deviation*Dist())
            mutations += 1
        end
    end

    Fom = f(xom, pDict, getGradientToo=false)

    om = (idx=-1, x=xom, F=Fom)

    return om, mutations
end

function decideAndAdd(p1, p2, om, Xkp1, Fkp1, popAdded, survived;
    parentsSurvive::Bool = true,
    verbose::Bool = false)

    actions = Dict(
        :bothParentsSurvived => 0,
        :crossover => 0,
        :genFitnessImproved => 0,
        :genFitnessNotImproved => 0,
        :mutation => 0,
        :noParentSurvived => 0,
        :onlyOneParentSurvived => 0,
        :parentsSelected => 0)

    n, p = size(Xkp1)

    childrenAdded = 0
    parentsAdded = 0

    if popAdded < p
        popAdded += 1 # add the child
        Xkp1[:, popAdded] = om.x
        Fkp1[popAdded] = om.F
        childrenAdded += 1
        myprintln(verbose, "Adding child")
    end

    if popAdded < p && parentsSurvive
        if survived[p1.idx] == 0
            popAdded += 1 # add the first parent
            Xkp1[:, popAdded] = p1.x
            Fkp1[popAdded] = p1.F
            survived[p1.idx] = 1
            myprintln(verbose, "Adding Parent 1")
        else
            myprintln(verbose, "Parent 1 already added")
        end
        parentsAdded += 1
    end

    if popAdded < p && parentsSurvive
        if survived[p2.idx] == 0
            popAdded += 1 # add the second parent
            Xkp1[:, popAdded] = p2.x
            Fkp1[popAdded] = p2.F
            survived[p1.idx] = 1
            myprintln(verbose, "Adding Parent 2")
        else
            myprintln(verbose, "Parent 2 already added")
        end
        parentsAdded += 1
    end

    if parentsAdded == 2
        actions[:bothParentsSurvived] = 1
    elseif parentsAdded == 1
        actions[:onlyOneParentSurvived] = 1
    elseif childrenAdded == 1
        actions[:noParentSurvived] = 1
    elseif childrenAdded == 0
        @error "floc"
    end

    return Xkp1, Fkp1, popAdded, survived, actions
    
end

# pdf = [0.22915395984352244, 0.03417449314058, 0.13089575337747889, 0.22024884567883649, 0.02160210334521976, 0.14400617789385242, 0.21991866672050991]

# cdf = [0.22915395984352244
# 0.26332845298410246
# 0.3942242063615813
# 0.6144730520404178
# 0.6360751553856376
# 0.78008133327949
# 0.9999999999999998];

# cdf[end] = 1.0

# selectTwoUniqueIndices(pdf)
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