include("./generateHaltonSequence.jl")
include("./generateRandomSequence.jl")
include("./generateLatinHypercubeSequence.jl")
include("./plotSampledSpace.jl")

"""
    sampleSpace(n, p; method::String="Halton", discard::Int=0, seed::Int=1)

Generate a sample space using different sampling methods.

# Arguments
- `n::Int`: The dimensionality of each point in the sample space.
- `p::Int`: The number of points to sample.
- `method::String`: The sampling method to use. Can be "Halton", "Latin Hypercube", or "Random". Default is "Halton".
- `discard::Int`: The number of initial points to discard. This is used with the "Halton" method. Default is 0.
- `seed::Int`: The seed for random number generation. This is used with the "Latin Hypercube" and "Random" methods. Default is 1.

# Returns
- `sampledSpace`: A matrix or array representing the sampled space.

# Note
This function calls `sampleSpaceHalton`, `sampleSpaceLatinHypercube`, or `sampleSpaceRandom` depending on the specified `method`.

# Errors
Throws an error if an unsupported `method` is provided.

# Examples
```julia
sampleSpace(2, 100)  # Samples 100 points in 2 dimensions using the Halton sequence
sampleSpace(3, 50, method="Latin Hypercube", seed=123)  # Samples 50 points in 3 dimensions using Latin Hypercube sampling with a specific seed
```
"""
function sampleSpace(n, p; 
    method::String="Halton",
    discard::Int=0,
    seed::Int=1)

    if method == "Halton"
        sampledSpace = sampleSpaceHalton(n, p, discard=discard)
    elseif method == "Latin Hypercube"
        sampledSpace = sampleSpaceLatinHypercube(n, p, seed=seed)
    elseif method == "Random"
        sampledSpace = sampleSpaceRandom(n, p, seed=seed)
    else
        @error "floc"
    end

    return sampledSpace
end