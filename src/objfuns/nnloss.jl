using CSV
using DataFrames

include("objective.jl")

function nnloss(w::Vector{Float64}, 
    p;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)

    trainData = p.trainData
    trainClass = p.trainClass
    testData = p.testData
    testClass = p.testClass
    
    dims = p.dims
    classify = p.classify


    data = p.data
    M = size(data, 1)
    nf = size(data, 2) - 1
    y = data[:, nf+1]
    xf = data[:, 1:nf]
    t = xf

    n = length(x) # note that df's x (time) is different from function parameter x
    f = 0.0
    g = zeros(n)
    
    if getGradientToo
        return f, g
    else
        return f
    end
end


rawDataFolder = "rawData/"
filename = rawDataFolder * "Indian Liver Patient Dataset (ILPD).csv"
df = CSV.File(filename) |> DataFrame
rename!(df, [:x, :y])
data = Matrix(df)
# x0 = [13.8, 8.3, 0.022, 1800, 900, 4.2]

objective = nnloss;

pr = generate_pr(objective, x0, data=data)