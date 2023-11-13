using CSV
using DataFrames
import MLJ as MLJ

include("objective.jl")
include("preprocessLiverData.jl")

function nnloss(w::Vector{Float64}, 
    p;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)

    X = p[:trainData]
    y = p[:classData]
    
    dims = p[:dims]
    classify = p[:classify]
    

    nf, M = size(X)

    n = length(w) # this time w is the main hero

    f = 0.0
    g = zeros(n)
    
    if getGradientToo
        return f, g
    else
        return f
    end
end

dims = [10, 10, 10, 1]

w0 = construct_w0_from_dims(dims)

classify = true
df = preprocessLiverData()


# DF.describe(df)

MLJ.schema(df)

# vscodedisplay(df)

df_train, df_test = MLJ.partition(df, 0.7, rng=123);
yTrain, XTrain = MLJ.unpack(df_train, ==(:Diagnosis));
trainData = Matrix(XTrain)
classData = Vector(yTrain)

params = Dict(:trainData=>trainData, :classData=>classData, :dims=>dims, :classify=>classify)



objective = nnloss;

pr = generate_pr(objective, w0, params=params)