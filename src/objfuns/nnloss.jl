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

import MLJ as MLJ
# using  CSV
# import DataFrames as DF

rawDataFolder = "rawData/";
filename = rawDataFolder * "Indian Liver Patient Dataset (ILPD).csv";

header = ["Age", "Gender", "TB", "DB", "Alkphos", "Sgpt", "Sgot", "TP", "ALB", "A/G Ratio", "Diagnosis"];

df00 = CSV.File(filename, header=header) |> DataFrame;


# df0 = DF.dropmissing(df00);
df0 = dropmissing(df00);

df = deepcopy(df0);
# Convert "Gender" to 0 and 1
df.Gender = map(gender -> gender == "Male" ? 0 : 1, df.Gender);

# Convert "Diagnosis" from 1 to 1/3 and 2 to 2/3
df.Diagnosis = map(diagnosis -> diagnosis == 1 ? 1/3 : 2/3, df.Diagnosis);

# Normalize the first 10 columns
for col in names(df)[1:10]  # Adjust the 10 here if the number of columns to normalize changes
    min_val, max_val = minimum(df[!, col]), maximum(df[!, col]);
    df[!, col] = map(x -> (x - min_val) / (max_val - min_val), df[!, col]);
end

# DF.describe(df)

MLJ.schema(df)

# vscodedisplay(df)

df_train, df_test = MLJ.partition(df, 0.7, rng=123);

yTrain, xTrain = MLJ.unpack(df_train, ==(:Diagnosis));

data = Matrix(df)


# x0 = [13.8, 8.3, 0.022, 1800, 900, 4.2]

objective = nnloss;

pr = generate_pr(objective, x0, data=data)