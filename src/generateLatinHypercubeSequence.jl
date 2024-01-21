include("./generateRandomSequence.jl")

using Random

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

function sampleSpaceLatinHypercube(n, p;
    seed::Int=0)
    
    B = generateRandomPermutations(n, p, seed=seed)
    R = sampleSpaceRandom(n, p, seed=seed)
    sampledSpace = (1/p)*(B-R)
    return sampledSpace 
end