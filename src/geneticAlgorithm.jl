include("helperFunctions.jl")

function deriveNextGeneration(Xk,
    Fk, 
    f::Function, 
    pDict::Dict;
    delta = 0.1,
    deviation = 0.1,
    Dist = randn,
    parentsSurvive = true,
    verbose::Bool=false)

    n, p = size(Xk)
    @show p
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

        p1, p2 = selectParents(Xk, Fk)
        fevals += p
        actions[:parentsSelected] += 1

        o = crossover(p1, p2) # offspring
        actions[:crossover] += 1

        om, mutations_1mut = mutation(o, f, pDict, 
        delta=delta, deviation=deviation, Dist=Dist)
        fevals += 1
        actions[:mutation] += mutations_1mut

        println("popAdded = $popAdded before decideAndAdd()")

        # select the child, and if choosing to let the parent survive by default, the best min(2, p-popAdded) parents
        Xkp1, Fkp1, popAdded, survived, actions_1dAA = 
        decideAndAdd(p1, p2, om, Xkp1, Fkp1, popAdded, survived, verbose=true) 
        actions = merge(+, actions, actions_1dAA)

        println("popAdded = $popAdded after decideAndAdd()")

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

    X0, F0 = fittestFirst(X0, F0)

    return X0, F0

end

"""
    fittestFirst(Xk::Matrix, Fk::Vector)

Reorders a population matrix `Xk` and its corresponding fitness vector `Fk` to ensure that the fittest individual is positioned first. This operation is useful in genetic algorithms and other evolutionary strategies to easily track and access the current best solution.

# Arguments
- `Xk`: A matrix representing the current population, where each column is an individual's parameter vector.
- `Fk`: A vector of fitness values corresponding to each individual in `Xk`. Higher values indicate better fitness, assuming a maximization problem.

# Returns
- `Xk`: The reordered population matrix with the fittest individual in the first column.
- `Fk`: The reordered fitness vector with the fitness value of the fittest individual in the first index.

# Method
1. Identify the index of the fittest individual based on the maximum fitness value in `Fk`.
2. If the fittest individual is not already in the first position, swap the first column of `Xk` with the column corresponding to the fittest individual. Perform a similar swap for the first element of `Fk` with the element at the fittest individual's index.

# Usage Example
```julia
# Assume a population matrix `Xk` and their fitness values `Fk`
Xk = [1.0 2.0; 3.0 4.0]  # Example population matrix
Fk = [0.5, 0.7]          # Corresponding fitness values

# Apply the function to reorder the population
Xk, Fk = fittestFirst(Xk, Fk)

# `Xk` and `Fk` now have the fittest individual in the first column/index
```
"""
function fittestFirst(Xk::Matrix, Fk::Vector)
    # Identify the index of the fittest individual
    fittest_index = argmax(Fk)
    # Swap to bring the fittest individual to the front if necessary
    if fittest_index != 1
        Xk[:, [1, fittest_index]] = Xk[:, [fittest_index, 1]]
        Fk[1], Fk[fittest_index] = Fk[fittest_index], Fk[1]
    end

    return Xk, Fk
end

"""
    selectParents(Xk, Fk)

Selects two unique parents from a population for crossover in a genetic algorithm. The selection is based on the fitness values of the population members, ensuring that individuals with higher fitness have a greater chance of being selected.

# Arguments
- `Xk`: A matrix representing the current population, where each column is an individual in the population space.
- `Fk`: A vector of fitness values corresponding to each individual in `Xk`. Higher values indicate better fitness.

# Returns
- `p1`, `p2`: Two tuples representing the selected parents. Each tuple contains:    
    - `idx`: The index of the individual in the original population.
    - `x`: The parameter vector of the selected individual.
    - `F`: The fitness value of the selected individual.

# Method
1. Normalize the fitness values `Fk` to create a probability density function (pdf), where the probability of selecting each individual is proportional to its fitness.
2. Use the `selectTwoUniqueIndices` function to choose two unique indices based on the pdf, ensuring diversity in parent selection.
3. Extract the parameter vectors (`x`) and fitness values (`F`) for the selected indices from `Xk` and `Fk`, respectively.
4. Package the selected individuals' information into tuples (`p1`, `p2`) for easy use in crossover and mutation operations.

# Examples
```julia
# Assume a population matrix `Xk` and their fitness values `Fk`
Xk = [1.0 2.0; 3.0 4.0]  # Example population
Fk = [0.5, 0.7]          # Example fitness values

# Select two parents from the population
parent1, parent2 = selectParents(Xk, Fk)

# Display the selected parents' details
println("Parent 1: ", parent1)
println("Parent 2: ", parent2)
```
"""
function selectParents(Xk, Fk)

    pdf = Fk./sum(Fk)

    idx1, idx2 = selectTwoUniqueIndices(pdf)

    x1, x2 = Xk[:, idx1], Xk[:, idx2]
    F1, F2 = Fk[idx1], Fk[idx2]
    
    p1 = (idx=idx1, x=x1, F=F1)
    p2 = (idx=idx2, x=x2, F=F2)

    return p1, p2

end

# p1, p2 = selectParents(Xk, Fk);
# p1.x
# p2.x

"""
    selectTwoUniqueIndices(pdf::Vector)

Selects two unique indices based on the probability density function (pdf) of selection probabilities.

# Arguments
- `pdf`: A vector representing the probability density function for parent selection probabilities.

# Returns
- A tuple containing two unique indices selected based on the pdf.

# Method
1. Select the first index using the pdf and convert the pdf to a cdf for selection.
2. Set the selected index's probability in the pdf to zero and renormalize the pdf.
3. Convert the updated pdf to a cdf and select the second index.

# Examples
```julia
pdf = [0.1, 0.2, 0.3, 0.4] # Example pdf for a population
idx1, idx2 = selectTwoUniqueIndices(pdf)
```
"""
function selectTwoUniqueIndices(pdf::Vector)
    # Convert pdf to cdf for the first selection
    # Select the first index
    cdf = cumsum(pdf)
    cdf[end] = 1.0
    # @show cdf
    rand_val1 = rand()
    # @show rand_val1
    # @show idx1 = findfirst(≥(rand_val1), cdf)
    idx1 = findfirst(≥(rand_val1), cdf)


    # Set the probability of the selected index to zero and renormalize pdf
    pdfm1 = deepcopy(pdf)
    pdfm1[idx1] = 0
    pdfm1 /= sum(pdfm1)

    # Convert the updated pdf to cdf for the second selection
    cdf = cumsum(pdfm1)
    # @show cdf
    cdf[end] = 1.0

    # Select the second index
    rand_val2 = rand()
    # @show rand_val2
    idx2 = findfirst(≥(rand_val2), cdf)

    return idx1, idx2

end

"""
    crossover(p1, p2)

Performs crossover between two parent individuals to produce a single offspring. The offspring's attributes are a weighted average of the parents' attributes, with weights proportional to the parents' fitness values.

# Arguments
- `p1`: The first parent, represented as a tuple containing the parent's index (`idx`), parameter vector (`x`), and fitness value (`F`).
- `p2`: The second parent, similarly represented.

# Returns
- `o`: The offspring, structured as a tuple containing its index (`idx`), parameter vector (`x`), and fitness value (`F`).

# Method
1. Extracts the indices, parameter vectors, and fitness values of the two parents.
2. Calculates the weighting factor `theta` based on the relative fitness values of the parents. The parent with higher fitness has a greater influence on the offspring's attributes.
3. Computes the offspring's parameter vector `x0` as a weighted average of the parents' parameter vectors, using `theta`.
4. Initializes the offspring's index (`idx`) to `-1` and its fitness value (`F`) to `-1`. The index is set to `-1` to indicate that the offspring is not a part of the current generation, distinguishing it from its parents. The fitness value is also set to `-1`, an impossible value under the fitness formula, signifying that the offspring's fitness is yet to be evaluated, typically after mutation operations.

# Usage Example
```julia
# Define two parents
p1 = (idx=1, x=[1.0, 2.0], F=5.0)
p2 = (idx=2, x=[2.0, 3.0], F=3.0)

# Perform crossover to generate an offspring
offspring = crossover(p1, p2)

# offspring's parameter vector is a weighted average of p1 and p2's vectors
```
"""
function crossover(p1, p2)
    
    idx1, x1, F1 = @unpack idx, x, F = p1
    idx2, x2, F2 = @unpack idx, x, F = p2
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

    @show popAdded 

    if popAdded < p
        popAdded += 1 # add the child
        Xkp1[:, popAdded] = om.x
        Fkp1[popAdded] = om.F
        childrenAdded += 1
        myprintln(verbose, "Adding child")
    end

    if popAdded < p
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

    if popAdded < p
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

    # @show popAdded 
    @show parentsAdded

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