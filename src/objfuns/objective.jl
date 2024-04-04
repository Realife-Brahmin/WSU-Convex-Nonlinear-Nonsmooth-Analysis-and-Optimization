using Base.Threads
using DataFrames
using LinearAlgebra
using Parameters

include("../helperFunctions.jl");
include("boxConstraintPenalty.jl");
include("../../alg.jl"); # include alg.jl from parent directory

"""
    generate_pr(functionName::Function, x0::Vector{Float64};
                data::Matrix{Float64} = Matrix{Float64}(undef, 0, 0),
                params::Vector{Float64} = Float64[])

Generate a named tuple `pr` representing a problem with the provided parameters.

# Arguments
- `functionName::Function`: The objective function for the problem.
- `x0::Vector{Float64}`: The initial point for the optimization problem.

# Keyword Arguments
- `data::Matrix{Float64}`: A matrix of data associated with the problem. Defaults to an empty matrix.
- `params`: Parameters for the objective function. Defaults to an empty vector.

# Returns
- A named tuple `pr` with fields:
    * `p`: A named tuple containing `params` and `data`.
    * `x0`: The initial point for the optimization problem.
    * `objective`: The objective function.
    * `alg`: The algorithm settings (of type `AlgorithmSettings`) to be used for the optimization problem. It is assumed to be a globally accessible instance of `AlgorithmSettings`.

# Example
```julia-repl
f(x) = sum(x.^2)
x0 = [1.0, 2.0, 3.0]
problem = generate_pr(f, x0)
```
"""
function generate_pr(functionName::Function,
    x0::Vector{Float64};
    problemType::String="Unconstrained",
    data=Matrix{Float64}(undef, 0, 0),
    params=Dict()
)

    p = Dict(:params => params, :data => data)

    println("Problem (pr::NamedTuple) generated")
    pr = (p=p, x0=x0, objective=functionName, alg=alg, problemType=problemType)
    return pr
end