include("objective.jl")

"""
    QPObjectiveFunction(x::Vector, pDict; verbose=false, log=true, getGradientToo=false)

Compute the objective function for a quadratic programming problem, optionally returning the gradient as well.

# Arguments
- `x::Vector`: The decision variable vector at which the objective and possibly the gradient are evaluated.
- `pDict::Dict`: A dictionary containing the necessary parameters for the quadratic form:
    - `G`: The symmetric matrix in the quadratic term of the objective function.
    - `c`: The coefficient vector for the linear term.
    - `c0`: Optional scalar representing the constant term (default is 0.0).

# Keyword Arguments
- `verbose::Bool=false`: If set to true, enables additional console output (currently not implemented).
- `log::Bool=true`: If set to true, logs function calls and results (currently not implemented).
- `getGradientToo::Bool=false`: If set to true, the function also returns the gradient of the objective function at `x`.

# Returns
- `f`: The computed value of the quadratic objective function at `x`.
- `g` (optional): The gradient of the objective function at `x`, returned if `getGradientToo` is true.

# Notes
- Ensure `pDict` contains all required keys (`G`, `c`, and optionally `c0`). If `c0` is not provided, it defaults to 0.0.
- The function throws an error if any required data is missing or if an inconsistency in input types is detected.

# Example
```julia
# Define the QP problem parameters
G = [2 0; 0 2]
c = [1; 1]
c0 = 0.5
pDict = Dict(:G => G, :c => c, :c0 => c0)

# Decision variable vector
x = [3; 4]

# Compute the objective function value and gradient
f, g = QPObjectiveFunction(x, pDict)
println("Objective value: ", f)
println("Gradient: ", g)
```
"""
function QPObjectiveFunction(x::Vector,
    pDict;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=false)

    n = length(x)
    @unpack G, c = pDict
    if haskey(pDict, :c0)
        c0 = pDict[:c0]
    else
        c0 = 0.0
    end

    f = 1 // 2 * transpose(x) * G * x + transpose(x) * c + c0

    if !getGradientToo
        return f
    elseif getGradientToo
        g = G*x + c
        return f, g
    else
        @error "floc"
    end

end