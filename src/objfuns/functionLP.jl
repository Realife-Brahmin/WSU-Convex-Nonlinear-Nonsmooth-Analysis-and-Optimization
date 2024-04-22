include("objective.jl")
include("../solveLP.jl")

using HiGHS
using JuMP
using Parameters

"""
    LPObjectiveFunction(x::Vector{Float64}, pDict; verbose=false, log=true, getGradientToo=true)

Evaluate the objective function and optionally its gradient for a linear programming problem.

# Arguments
- `x::Vector{Float64}`: The decision variable vector at which to evaluate the function.
- `pDict::Dict`: A dictionary containing problem-specific data with key `:c` for the objective function's coefficients.

# Keyword Arguments
- `verbose::Bool=false`: Flag to enable verbose output.
- `log::Bool=true`: Flag to enable logging of the function evaluations.
- `getGradientToo::Bool=true`: Flag to indicate whether the gradient should also be returned.

# Returns
- `f::Float64`: The evaluated objective function value at `x`.
- `g::Vector{Float64}`: (Optional) The gradient of the objective function at `x`.

# Description
This function computes the objective function value of a linear programming problem, which is simply the inner product of the coefficient vector `c` and the decision variable vector `x`. If the `getGradientToo` flag is set to `true`, it also returns the gradient, which in the case of a linear objective function, is constant and equal to `c`.

# Notes
- If `getGradientToo` is `false`, only the objective function value is returned.
- This function will throw an error if neither of the conditions for `getGradientToo` is met.

# Examples
```julia
# Define coefficients and decision variables
c = [1.0, 2.0, 3.0]
x = [4.0, 5.0, 6.0]
pDict = Dict(:c => c)

# Evaluate objective function value and gradient
f, g = LPObjectiveFunction(x, pDict, getGradientToo=true)

# The output will be the objective function value and the gradient vector
```
"""
function LPObjectiveFunction(x::Vector{Float64}, 
    pDict;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)
    
    n = length(x)
    @unpack c = pDict
    f = transpose(x)*c

    if !getGradientToo
        return f
    elseif getGradientToo
        g = c
        return f, g
    else
        @error "floc"
    end
    
end
