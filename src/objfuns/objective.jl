using Base.Threads
using DataFrames
using LinearAlgebra
using Parameters

include("../helperFunctions.jl");
include("boxConstraintPenalty.jl");
include("equalityConstrainedQP.jl");
include("../../alg.jl"); # include alg.jl from parent directory
include("../types.jl")

"""
    generate_pr(functionName, x0; problemType="Unconstrained", data=undef, params=Dict(), objectiveString="undefinedFunction", verbose=true)

Generate a problem representation (`pr`) as a `NamedTuple` for optimization or problem-solving algorithms.

# Arguments
- `functionName::Function`: The main function to be optimized or solved. This function should take `x0` as input and return a scalar value representing the objective or cost.
- `x0::Vector{Float64}`: The initial guess or starting point for the optimization process.

# Keyword Arguments
- `problemType::String="Unconstrained"`: Describes the type of problem. Default is "Unconstrained". Could be set to other values like "Constrained", depending on the problem.
- `data::Matrix{Float64}`: Additional data required by `functionName` for evaluating the objective. Defaults to an empty `Matrix{Float64}(undef, 0, 0)`.
- `params::Dict`: A dictionary of parameters that `functionName` might need to utilize. Defaults to an empty dictionary.
- `objectiveString::String="undefinedFunction"`: A string representation of the objective function's name. If not provided, it defaults to the name of `functionName`.
- `verbose::Bool=true`: Controls the verbosity of the function. If `true`, the function will print messages about the problem generation process.

# Returns
- `NamedTuple`: Contains the generated problem representation with keys:
    - `p`: A dictionary containing `params` and `data`.
    - `x0`: The initial guess for the optimization.
    - `objective`: The `functionName` provided as the objective function.
    - `alg`: Undefined in the provided signature but presumably intended for algorithm specification.
    - `problemType`: The type of problem.
    - `objectiveString`: The string representation of the objective function.

# Notes
- Ensure that `functionName` is compatible with the structure of `x0` and uses `data` and `params` correctly if they are provided.
- The `alg` key in the returned `NamedTuple` is mentioned in the documentation but not defined in the provided function signature. Make sure to define or remove it as necessary.
- The function automatically assigns the `functionName` to `objectiveString` if the latter is not explicitly provided, ensuring consistency in problem representation.

# Example
```julia
function myObjective(x)
    return sum(x.^2)
end

x0 = [1.0, 2.0, 3.0]
pr = generate_pr(myObjective, x0, problemType="Unconstrained", verbose=true)
```
"""
function generate_pr(functionName::Function,
    x0::Vector{Float64};
    problemType::String="Unconstrained",
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

    pr = (p=p, x0=x0, objective=functionName, alg=alg, problemType=problemType, objectiveString=objectiveString)

    myprintln(verbose, "Problem (pr::NamedTuple) generated")

    return pr
end