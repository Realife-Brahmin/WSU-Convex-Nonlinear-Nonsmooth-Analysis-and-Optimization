using Random

"""
    sampleSpaceRandom(n, p; seed::Int=1)

Generate a random matrix of size n by p from a specified sample space.

This function generates random numbers within a sample space and returns them as a matrix of size n by p.

# Arguments:
- `n::Int`: The number of rows in the random matrix.
- `p::Int`: The number of columns in the random matrix.
- `seed::Int=0`: (Optional) The seed for the random number generator. Default is 0 for reproducibility. Must be non-negative.

# Returns:
- A matrix of size n by p containing random numbers generated from the specified sample space.

Note: The `seed` argument allows you to control the reproducibility of the random numbers. Using the same seed will produce the same random matrix.

# Example:
```julia
julia> randomMatrix = sampleSpaceRandom(3, 4, seed=42)
3×4 Matrix{Float64}:
    0.366025  0.794213  0.237254  0.512295
    0.731301  0.390758  0.573887  0.763487
    0.316836  0.474608  0.22422   0.281135
```
"""
function sampleSpaceRandom(n, p;
    seed::Int=0)
    rng = MersenneTwister(seed)
    sampledSpace = rand(rng, n, p)
    return sampledSpace
end
