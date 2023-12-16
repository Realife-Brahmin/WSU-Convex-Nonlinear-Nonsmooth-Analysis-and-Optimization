using CSV
using DataFrames
using Parameters
using Statistics

include("objective.jl")

function receiverLocation(
    x::Vector{Float64}, 
    p,
    verbose::Bool = false,
    log::Bool = true,
    getGradientToo::Bool = true)

    params = p[:params]
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
# datasetName = "SD10" 
ext = ".csv"

filename = rawDataFolder * datasetName * ext
df0 = CSV.File(filename, header=false) |> DataFrame
df = Matrix(df0)
println(size(df))
m = size(df, 2)
n = size(df, 1) - 1
μ = df[1, 1:m]
P = df[2:n+1, 1:m]

pmean = vec(mean(P, dims=2))
x0 = Float64.(pmean)

params = Dict(:pmean => pmean, :P => P, :μ => μ, :datasetName => datasetName)
objective = receiverLocation;

pr = generate_pr(objective, x0, params=params)

obj = pr.objective

# Testing for f, g values for x = p_k for some k ∈ 1:m
# x0 = P[:, rand(1:m)]
# x0 = pmean
# f, g = obj(x0, pr.p)