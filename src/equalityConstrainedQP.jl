"""
    equalityConstrainedQP(w::Vector{Float64}, pDict; verbose=false, log=true, getGradientToo=false) -> Union{Number, Tuple}

Evaluates the objective function of an equality-constrained quadratic programming (QP) problem at a given point `w`, with an option to also return the gradient.

This function specifically addresses the standard form of a quadratic problem:
minimize    1/2 * w^T * G * w + w^T * c
subject to  A * w = b

# Arguments
- `w::Vector{Float64}`: The point at which the QP objective function and possibly its gradient are to be evaluated.
- `pDict`: A dictionary containing problem-specific parameters with keys:
    - `G`: The quadratic term matrix in the objective function.
    - `A`: The matrix representing linear equality constraints.
    - `c`: The linear term vector in the objective function.
    - `b`: The right-hand side vector of the equality constraints.

# Keyword Arguments
- `verbose::Bool=false`: If `true`, enables printing of function-specific messages. Currently not utilized within the function.
- `log::Bool=true`: If `true`, activates logging of execution details. Currently not utilized within the function but can be implemented for logging purposes.
- `getGradientToo::Bool=false`: Controls whether the gradient of the objective function with respect to `w` should also be returned. Defaults to `false`.

# Returns
- If `getGradientToo` is `false`, returns the value of the objective function (`f`) at `w`.
- If `getGradientToo` is `true`, returns a tuple (`f`, `g`), where `f` is the value of the objective function at `w`, and `g` is the gradient of the objective function at `w`.

# Example
```julia
# Define problem parameters
G_matrix = [2 0; 0 2]
A_matrix = [1 1]
c_vector = [-2; -5]
b_vector = [0]

# Construct parameters dictionary
paramsDict = Dict(:params => Dict(:G => G_matrix, :A => A_matrix, :c => c_vector, :b => b_vector))

# Point to evaluate
w = [1.0, 3.0]

# Calculate only the objective function value
fval = equalityConstrainedQP(w, paramsDict)

# Calculate both the objective function value and its gradient
fval, grad = equalityConstrainedQP(w, paramsDict, getGradientToo=true)
```
"""
function equalityConstrainedQP(w::Vector{Float64},
    pDict;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=false) # yes, false by default

    @unpack G, A, c, b = pDict[:params]
    f = 1//2*transpose(w)*G*w + transpose(w)*c

    if getGradientToo
        g = G*w + c
        return f, g
    else
        return f
    end
    
end

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

"""
    normalizeLinearConstraints(B, d) -> (Matrix, Vector)

Normalize the rows of matrix `B` and the elements of vector `d` based on the norm of each row in `B`.

This function calculates the Euclidean norm of each row in `B` and divides each element in the row by its norm. Similarly, each corresponding element in `d` is divided by the norm of its respective row in `B`. This normalization process ensures that each constraint represented by a row in `B` and the corresponding element in `d` is scaled to have a unit norm.

# Arguments
- `B`: A matrix where each row represents a linear constraint's coefficients.
- `d`: A vector where each element corresponds to a linear constraint's right-hand side value.

# Returns
- `B_normalized`: A matrix where each row of `B` has been normalized.
- `d_normalized`: A vector where each element of `d` has been normalized according to its respective row's norm in `B`.

# Example
```julia
B = [3 4; 1 2]
d = [5; 2]
B_normalized, d_normalized = normalizeLinearConstraints(B, d)
```
"""
function normalizeLinearConstraints(B, d)
    norms = norm.(eachrow(B))
    B_normalized = B ./ norms
    d_normalized = d ./ norms
    return B_normalized, d_normalized
end
