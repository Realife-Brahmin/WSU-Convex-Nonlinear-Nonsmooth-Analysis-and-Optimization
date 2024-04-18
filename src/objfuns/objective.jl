using Base.Threads
using DataFrames
using LinearAlgebra
using Parameters

include("../helperFunctions.jl");
# include("activeSetQuadraticProgramming")
include("boxConstraintPenalty.jl");
include("../equalityConstrainedQP.jl");
include("../../alg.jl"); # include alg.jl from parent directory
include("../types.jl")

"""
    generate_pr(functionName::Function, x0::Vector{Float64};
                problemType::String="Unconstrained",
                method::String="QuasiNewton",
                data=Matrix{Float64}(undef, 0, 0),
                params=Dict(),
                objectiveString::String="undefinedFunction",
                verbose::Bool=true) -> NamedTuple

Constructs and returns a named tuple representing an optimization problem, configured according to specified parameters and settings. This representation includes all necessary details such as the objective function, initial guess, algorithm settings, and additional problem-specific parameters.

# Arguments
- `functionName::Function`: The function to be minimized or maximized. This function should take a vector `x` as input and return a scalar value.
- `x0::Vector{Float64}`: The initial guess for the optimizer, provided as a vector.

# Keyword Arguments
- `problemType::String="Unconstrained"`: Specifies the type of optimization problem (e.g., "Unconstrained", "Constrained").
- `method::String="QuasiNewton"`: Specifies the optimization method to be used (e.g., "QuasiNewton", "GradientDescent").
- `data::Matrix{Float64}`: An optional matrix of additional data needed by `functionName` for evaluating the objective. Defaults to an empty matrix if unspecified.
- `params::Dict`: A dictionary of additional parameters that might be required by the objective function or the optimization algorithm.
- `objectiveString::String="undefinedFunction"`: A string representation of the objective function's name for logging or identification purposes. Defaults to the name of `functionName` if not specified.
- `verbose::Bool=true`: Controls whether to print detailed messages about the problem generation process.

# Returns
- `NamedTuple`: A named tuple containing the following fields:
    - `p`: A dictionary containing `params` and `data`.
    - `x0`: The initial guess vector.
    - `objective`: The objective function.
    - `alg`: Algorithm settings created based on `problemType` and `method`.
    - `problemType`: The type of problem.
    - `objectiveString`: The string representation of the objective function.

# Example
```julia
function myObjective(x)
    return sum(x.^2)
end

x0 = [1.0, 2.0, 3.0]
params = Dict("lambda" => 0.01)
pr = generate_pr(myObjective, x0, params=params, method="GradientDescent", verbose=true)
```
"""
function generate_pr(functionName::Function,
    x0::Vector{Float64};
    problemType::String="Unconstrained",
    method::String="QuasiNewton",
    data=Matrix{Float64}(undef, 0, 0),
    params=Dict(),
    objectiveString::String="undefinedFunction",
    verbose::Bool=true,
)

    p = Dict(:params => params, :data => data)

    if objectiveString == "undefinedFunction"
        objectiveString = string(functionName)
    elseif objectiveString != string(functionName)
        myprintln(verbose, "A wrapper function $(string(functionName)) will be used to solve for the original problem of $(objectiveString)")
    end

    alg = create_algorithm_settings(problemType=problemType, method=method)

    pr = (p=p, x0=x0, objective=functionName, alg=alg, problemType=problemType, objectiveString=objectiveString)

    myprintln(verbose, "Problem (pr::NamedTuple) generated")

    return pr
end