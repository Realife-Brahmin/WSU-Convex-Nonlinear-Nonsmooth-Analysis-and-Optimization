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