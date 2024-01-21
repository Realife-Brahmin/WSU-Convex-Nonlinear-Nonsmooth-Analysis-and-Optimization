include("./generateRandomSequence.jl")

using Random

"""
    generateRandomPermutations(n, p; seed::Int=0)

Generate `n` random permutations of numbers from 1 to `p`.

This function generates `n` random permutations of the numbers from 1 to `p` and arranges them as rows in a matrix.

# Arguments:
- `n::Int`: The number of random permutations to generate.
- `p::Int`: The range of numbers from 1 to `p` for the permutations.
- `seed::Int=0`: (Optional) The seed for the random number generator. Default is 0 for reproducibility. Must be non-negative.

# Returns:
- A matrix `B` of size `n` by `p` containing `n` random permutations.

# Example:
```julia
# Generate 5 random permutations of numbers from 1 to 10 with a seed of 42
permutations = generateRandomPermutations(5, 10, seed=42)
```
"""
function generateRandomPermutations(n, p;
    seed::Int=0)
    rng = MersenneTwister(seed)
    B = Matrix{Int}(undef, n, p)
    for i in 1:n
        permutation = randperm(rng, p)
        B[i, :] = permutation
    end
    return B
end

"""
    sampleSpaceLatinHypercube(n, p; seed::Int=0)

Generate a Latin Hypercube Sampled Space.

This function generates a Latin Hypercube sampled space of `n` points in `p` dimensions. It combines random permutations and random numbers to create the sampled space.

# Arguments:
- `n::Int`: The number of points in the Latin Hypercube.
- `p::Int`: The number of dimensions.
- `seed::Int=0`: (Optional) The seed for the random number generator. Default is 0 for reproducibility. Must be non-negative.

# Returns:
- A matrix `sampledSpace` of size `n` by `p`, representing the Latin Hypercube sampled space.

# Example:
```julia
# Generate a Latin Hypercube sampled space with 10 points in 3 dimensions and a seed of 123
lhs = sampleSpaceLatinHypercube(10, 3, seed=123)
```
"""
function sampleSpaceLatinHypercube(n, p;
    seed::Int=0)
    
    B = generateRandomPermutations(n, p, seed=seed)
    R = sampleSpaceRandom(n, p, seed=seed)
    sampledSpace = (1/p)*(B-R)
    return sampledSpace 
end