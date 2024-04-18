include("objective.jl")

"""
    QPObjectiveFunction(x::Vector{Float64}, pDict; verbose=false, log=true, getGradientToo=true) -> Union{Float64, Tuple{Float64, Vector{Float64}}}

Evaluates the objective function for a quadratic programming (QP) problem at a given point `x`, and optionally computes the gradient. The function is designed to work with quadratic forms, making it suitable for applications in optimization, control theory, and economics where quadratic models are prevalent.

# Arguments
- `x::Vector{Float64}`: The point at which the objective function (and possibly the gradient) is to be evaluated.
- `pDict`: A dictionary containing parameters `G` and `c`:
    - `G`: The symmetric matrix representing the quadratic part of the objective function.
    - `c`: The vector representing the linear part of the objective function.

# Keyword Arguments
- `verbose::Bool=false`: If set to `true`, enables detailed output useful for debugging or monitoring the computation. 
- `log::Bool=true`: If `true`, indicates that logging should be enabled, though specifics of the logging mechanism are not implemented in the snippet provided.
- `getGradientToo::Bool=true`: Controls whether the gradient of the objective function is also computed and returned. If `false`, only the objective function value is returned.

# Returns
- If `getGradientToo` is `false`, returns the value of the objective function `f` at `x`.
- If `getGradientToo` is `true`, returns a tuple `(f, g)`, where `f` is the objective function value and `g` is the gradient at `x`.

# Example
```julia
# Define the problem parameters
G = [2 0; 0 2]
c = [-1; -2]

# Create the parameters dictionary
paramsDict = Dict(:params => Dict(:G => G, :c => c))

# Point at which to evaluate the function
x = [1.0, 1.0]

# Calculate only the objective function value
f = QPObjectiveFunction(x, paramsDict, getGradientToo=false)

# Calculate both the objective function value and its gradient
f, g = QPObjectiveFunction(x, paramsDict)
```
"""
function QPObjectiveFunction(x::Vector{Float64},
    pDict;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)

    n = length(x)
    @unpack G, c = pDict[:params]
    f = 1 // 2 * transpose(x) * G * x + transpose(x) * c

    if !getGradientToo
        return f
    elseif getGradientToo
        g = G*x + c
        return f, g
    else
        @error "floc"
    end

end