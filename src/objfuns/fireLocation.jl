include("objective.jl")
include("../equalityConstrainedQP.jl")

using Parameters

"""
    fireLocation(x::Vector{Float64}, pDict; verbose=false, log=true, getGradientToo=true)

Calculate the squared residuals of fire location estimates from ranger stations and optionally return the gradient of the residuals.

This function evaluates the sum of squares of residuals from predicted fire locations based on inputs from ranger stations. It's designed to support optimization algorithms, particularly those requiring objective function values and gradients.

# Arguments
- `x::Vector{Float64}`: A vector representing the estimated coordinates of the fire.
- `pDict`: A dictionary containing problem-specific parameters. Must include:
    - `B`: The matrix representing the directional vectors from ranger stations towards the fire.
    - `d`: The vector representing distances or directional magnitudes towards the fire from each ranger station.

# Keyword Arguments
- `verbose::Bool=false`: If `true`, enables the printing of function-specific messages.
- `log::Bool=true`: If `true`, logs execution details. Currently not used within the function but can be implemented for logging purposes.
- `getGradientToo::Bool=true`: If `true`, the function also returns the gradient of the objective function with respect to `x`.

# Returns
- If `getGradientToo` is `false`, returns the sum of squares of residuals (`f`).
- If `getGradientToo` is `true`, returns a tuple (`f`, `g`) where `f` is the sum of squares of residuals, and `g` is the gradient of `f` with respect to `x`.

# Context
This function is utilized in the context of locating a fire given observations from ranger stations. These stations, at known coordinates, observe the fire at specific angles, forming a system of equations `A*x = b`, where `x` represents the fire's location. This function computes the total residual squares (`f`) for an estimated point `x` from the respective 'lines' determined by the observations and the stations' locations. It's particularly useful as part of an Extended Constrained Quadratic Programming (ECQP) solver or other optimization frameworks that require an objective function and its gradient.

# Example
```julia
# Assuming B and d are defined based on ranger station observations
x_estimate = [1.0, 2.0]  # Initial estimate of the fire's location
pDict = Dict(:B => B_matrix, :d => d_vector)

# Calculate only the objective function value
fval = fireLocation(x_estimate, pDict, getGradientToo=false)

# Calculate both the objective function value and its gradient
fval, grad = fireLocation(x_estimate, pDict)
```
"""
function fireLocation(x::Vector{Float64}, 
    pDict;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)
    
    n = length(x)
    @unpack B, d = pDict[:params]
    r = B*x - d
    f = 1//2 *sum(r.^2)

    if !getGradientToo
        return f
    elseif getGradientToo
        g = transpose(B)*r
        return f, g
    else
        @error "floc"
    end
    
end

angles_deg = [95, 84, 75, 63]
m_line = tand.(90 .- angles_deg)
coords = [(10, 12), (3, 9), (2, 3), (8, 1)]
coords_x, coords_y = first.(coords), last.(coords)
c_line = coords_y - m_line.*coords_x
B_unnormalized = hcat(-m_line, -ones(4))

m, n = size(B_unnormalized)
d_unnormalized = c_line

normalizeConstraints = false
# normalizeConstraints = true

if normalizeConstraints
    B, d = normalizeLinearConstraints(B_unnormalized, d_unnormalized)
else
    B, d =  B_unnormalized, d_unnormalized
end

x0 = B \ d
x0 = myfill(x0, 1) # the original x0 is basically the perfect solution. Instead, for checking the ECQP solver, let's take something else, making sure there are no shape errors of course.
r0 = B * x0 - d
w0 = vcat(x0, r0)

# params = Dict(:B=>B, :d=>d);

# objective = fireLocation;

# pr = generate_pr(objective, x0, params=params, problemType="ECQP")

#test run
# f0 = fireLocation(x0, params, getGradientToo=false)

# Henceforth, x, n, m lose their original meanings, in favour of w, n+m, n+m
G = vcat(hcat(zeros(n, n), zeros(n, m)), hcat(zeros(m, n), I(m)));
A = hcat(B, -I(m))
b = d
c = zeros(n+m)

pECQP = Dict(:G=>G, :c=>c, :A=>A, :b=>b )

objective = equalityConstrainedQP
objectiveOriginal = fireLocation
objectiveString = string(objectiveOriginal)
params = pECQP

pr = generate_pr(objective, w0, params=params, problemType="ECQP"; objectiveString=objectiveString)

# f0 = equalityConstrainedQP(w0, pECQP)

# solState = SolStatePGCGType(w0, G, c, A, fk=f0)
