using CSV
using DataFrames
using Parameters

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
        f -= 1/(2*m)*μ[k]*ln(sum((x - P[:, k]).^2))
    end

    if getGradientToo
        g = zeros(n)
        for i = 1:n
            g[i] = (x[i] - pmean[i])
            for k = 1:m
                g[i] -= (1/m)*μ[k]*(1/( sum( (x-P[:, k]).^2)) )*(x[i] - P[i, k])
            end
        end
        return f, g
    else
        return f
    end
end

rawDataFolder = "rawData/"
datasetName = "SD05"
ext = ".csv"
filename = rawDataFolder * datasetName * ext
df = CSV.File(filename, header=false) |> DataFrame
# rename!(df, [:x, :y])


# # x0 = Float64.(d)

# params = Dict(:pmean => pmean, :P => P, :mu => μ, :datasetName => datasetName)
# objective = receiverLocation;

# pr = generate_pr(objective, x0, params=params)

