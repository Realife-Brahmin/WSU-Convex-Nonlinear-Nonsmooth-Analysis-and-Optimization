include("sampleSpace.jl")

"""
    createInitialSimplexFromOnePoint(x0; deviationFactor::Float64 = 0.5, method::String = "Halton", seed::Int = 1, discard::Int = 0) -> Array{Float64,2}

Create an initial simplex for optimization algorithms starting from a single point `x0`.

# Arguments
- `x0::Array{Float64,1}`: The initial point from which the simplex is generated. This point should be a vector in the parameter space over which optimization is performed.
- `deviationFactor::Float64`: A factor that determines how far the generated points of the simplex deviate from the initial point `x0`. Default is `0.5`.
- `method::String`: The method used to sample points in the parameter space. Default is `"Halton"`, which refers to the Halton sequence for generating quasi-random points.
- `seed::Int`: The seed for the random number generator, ensuring reproducibility of the sample space. Default is `1`.
- `discard::Int`: The number of initial points to discard from the generated sequence. This can be used to avoid the initial sequence values that might be less uniformly distributed. Default is `0`.

# Functionality
This function generates an n+1 simplex from a single initial point `x0` in an n-dimensional space. It samples `n` other points based on the specified method, deviation factor, seed, and discard parameters. These sampled points are then renormalized to lie between -1 and 1. The deviation of these points from `x0` is scaled by the `deviationFactor` to form the vertices of the simplex, which are then returned as an `(n, n+1)` matrix.

# Returns
- `X0::Array{Float64,2}`: An `(n, n+1)` matrix where `n` is the dimensionality of the initial point `x0`. The first column of `X0` is the initial point `x0`, and the subsequent columns are the vertices of the simplex generated around `x0`.

# Example
```julia
x0 = [1.0, 2.0] # Initial point in a 2D space
simplex = createInitialSimplexFromOnePoint(x0)
```
"""
function createInitialSimplexFromOnePoint(x0;
    deviationFactor::Float64 = 0.5,
    method::String = "Halton",
    seed::Int = 1,
    discard::Int = 0)

    n = length(x0)
    sampledSpace = sampleSpace(n, n, method=method, seed=seed, discard=discard) # Create n other points, which will join our initial x0 to form the n+1 simplex
    reSampledSpace = -ones(size(sampledSpace)) + 2*sampledSpace # renormalizes the sample space points to lie between -1 to 1
    delta_X0_npoints = deviationFactor * (x0 .* reSampledSpace)
    X0 = zeros(n, n+1)

    X0[1:n, 1] = x0
    X0[1:n, 2:n+1] = x0 .+ delta_X0_npoints
    
    return X0

end

# x0 = pr.x0;
# X0 = createInitialSimplexFromOnePoint(x0, deviationFactor=0.1);