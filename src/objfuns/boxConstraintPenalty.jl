"""
    boxConstraintPenalty(x::Vector, box::Dict)

Calculates the penalty for violating box constraints in an optimization problem. The penalty is computed based on the distance of the decision variables from their specified bounds if they are outside those bounds.

# Arguments
- `x`: A vector of decision variables.
- `box`: A dictionary containing the box constraint parameters. It must have the following keys:
    - `indices`: The indices of the decision variables that have box constraints.
    - `lbs`: A vector of lower bounds for each constrained variable.
    - `ubs`: A vector of upper bounds for each constrained variable.
    - `barrier`: The barrier coefficient, a positive scalar that determines the penalty's steepness for violating the constraints.

# Returns
- `P`: The total penalty for all variables outside their box constraints. This value is 0 if all variables are within their bounds.

# Examples
```julia
x = [1.5, -0.5, 2.5]
box = Dict(
    :indices => [1, 2, 3],
    :lbs => [1.0, 0.0, 2.0],
    :ubs => [2.0, 1.0, 3.0],
    :barrier => 100.0
)
penalty = boxConstraintPenalty(x, box)
```
"""
function boxConstraintPenalty(x::Vector, box::Dict)

    @unpack indices, lbs, ubs, barrier = box

    P = 0.0
    for constraintNum ∈ eachindex(indices)
        i = indices[constraintNum]
        P += barrier * (min(0, x[i] - lbs[constraintNum])^2)
        P += barrier * (min(0, ubs[constraintNum] - x[i])^2)
    end
    return P
end

"""
    boxConstraintGradient(x::Vector, box::Dict)

Computes the gradient of the penalty function for box constraints at a given point. This function is useful for gradient-based optimization algorithms, enabling them to account for the penalty gradient due to box constraints.

# Arguments
- `x`: A vector of decision variables.
- `box`: A dictionary containing the box constraint parameters, structured similarly to `boxConstraintPenalty`.

# Returns
- `gP`: The gradient of the penalty with respect to the decision variables. Note: This implementation returns the sum of gradients across all constrained variables, not a gradient vector. To be used correctly in optimization, this needs adjustment to return a gradient vector corresponding to each decision variable.

# Examples
```julia
x = [1.5, -0.5, 2.5]
box = Dict(
    :indices => [1, 2, 3],
    :lbs => [1.0, 0.0, 2.0],
    :ubs => [2.0, 1.0, 3.0],
    :barrier => 100.0
)
gradient = boxConstraintGradient(x, box)
```
"""
function boxConstraintGradient(x::Vector, box::Dict)

    @unpack indices, lbs, ubs, barrier = box

    gP = 0.0

    for constraintNum ∈ eachindex(indices)
        i = indices[constraintNum]
        gP += -2 * barrier * min(0, x[i] - lbs[constraintNum])
        gP += -2 * barrier * min(0, ubs[constraintNum] - x[i])
    end

    return gP
end