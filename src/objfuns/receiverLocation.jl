using CSV
using DataFrames
using Parameters
using Statistics

include("objective.jl")

"""
    receiverLocation(x::Vector{Float64}, p; verbose::Bool = false, log::Bool = true, getGradientToo::Bool = true)

Calculates an optimal receiver location in a multi-transmitter environment. This function aims to find a balance between minimizing the distance from the mean position of transmitters and avoiding too close proximity to any individual transmitter.

# Parameters
- `x`: Input vector representing the receiver location in `Vector{Float64}`.
- `p`: Parameters including transmitter locations `P`, mean position `pmean`, and weight vector `μ`.

# Keyword Arguments
- `verbose` (default `false`): If `true`, prints additional information.
- `log` (default `true`): If `true`, enables logging for debugging or analysis.
- `getGradientToo` (default `true`): If `true`, computes and returns the gradient along with the objective function value.

# Returns
- Tuple `(f, g)` where `f` is the objective function value and `g` is the gradient, if `getGradientToo` is `true`.
- Only the objective function value `f` if `getGradientToo` is `false`.

# Details
The function computes the objective function value `f`, which represents a compromise between being near the mean position and not too close to any transmitter. The gradient `g` is calculated for optimization purposes. The function uses parameters `P`, `pmean`, and `μ` to adjust the calculations based on the transmitter locations and weights.

# Special Cases
- Returns `f(x) = ∞` and `∇f(x) = ∞` if `x` equals any transmitter location `pk` for robustness in optimization scenarios.
"""
function receiverLocation(
    x::Vector{Float64}, 
    p;
    verbose::Bool = false,
    log::Bool = true,
    getGradientToo::Bool = true)
    # function implementation
end

function receiverLocation(
    x::Vector{Float64}, 
    p;
    verbose::Bool = false,
    log::Bool = true,
    getGradientToo::Bool = true)

    # params = p[:params]
    params = p

    @unpack P, pmean, μ = params # P is n x m, pmean is n-vector, μ is m-vector

    n = length(x)
    m = length(μ)

    f = 0.5*sum((x - pmean).^2)
    for k = 1:m
        f -= 1/(2*m)*μ[k]*Base.log(sum((x - P[:, k]).^2))
    end

    if getGradientToo
        g = zeros(n)
        for i = 1:n
            g[i] = (x[i] - pmean[i])
            for k = 1:m
                squaredDistanceFrom_pk = sum( (x-P[:, k]).^2 )
                if squaredDistanceFrom_pk == 0
                    g[i] = Inf
                else
                    g[i] -= (1/m)*μ[k]*(1/squaredDistanceFrom_pk) *(x[i] - P[i, k])
                end
            end
        end
        return f, g
    else
        return f
    end
end

rawDataFolder = "rawData/"
datasetName = "SD05"
datasetName = "SD10" 
ext = ".csv"

filename = rawDataFolder * datasetName * ext
df0 = CSV.File(filename, header=false) |> DataFrame
df = Matrix(df0)
# println(size(df))
m = size(df, 2)
n = size(df, 1) - 1
μ = df[1, 1:m]
P = df[2:n+1, 1:m]

pmean = vec(mean(P, dims=2))
x0 = Float64.(pmean)
x0 = Float64.(ones(n))
@show trans1 = rand(1:Int(floor(m/2)))
@show trans2 = rand(Int(floor(m/2+1)):m)
x0 = mean( [ P[:, trans1], P[:, trans2] ] )

params = Dict(:pmean => pmean, :P => P, :μ => μ, :datasetName => datasetName)
objective = receiverLocation;

pr = generate_pr(objective, x0, params=params)

# obj = pr.objective
# x0 = rand(n)
# Testing for f, g values for x = p_k for some k ∈ 1:m
# x0 = P[:, rand(1:m)]
# x0 = pmean
# f, g = obj(x0, pr.p)
# f = obj(x0, pr.p, getGradientToo=false)

# min_value = Inf
# min_x0 = nothing

# for i in 1:10000
#     global min_value, min_x0, f, x0
#     x0 = rand(n)  # Random starting point
#     f = obj(x0, pr.p, getGradientToo=false)

#     if f < min_value
#         min_value = f
#         min_x0 = x0
#     end
# end

# println("Minimum value: ", min_value)
# println("Corresponding x0: ", min_x0)