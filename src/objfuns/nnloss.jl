using CSV
using DataFrames
import MLJ as MLJ

include("activation.jl")
include("construct_w0_from_dims.jl")
include("objective.jl")
include("preprocessLiverData.jl")

function nnloss(w::Vector{Float64}, 
    p;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)

    # p = p.params
    # params = p[:params]
    params = p
    trainData = params[:trainData]
    if size(trainData, 2) == 1
        X = reshape(trainData, length(trainData), 1)
    else
        X = trainData'
    end

    classData = params[:classData]

    if size(classData) == ()
        y = fill(classData, 1, 1)
    else
        y = reshape(classData, 1, length(classData))
    end
    
    dims = params[:dims]
    classify = params[:classify]
    
    nf, M = size(X, 1), size(X, 2)
    n = length(w)

    if M != length(y)
        @error "Incompatible sizes of X and y"
    elseif dims[1] != nf
        @error "Make sure first element of dims matches the number of used features from X"
    elseif dims[end] != 1
        @error "Currently only single class classification is supported. Make sure that dim[end] = 1."
    elseif n != getNumWeights(dims)
        @error "Mismatch between supplied w and the number of weights corresponding to dims."
    end

    d = length(dims) - 1

    f = 0.0

    Ws = Vector{}(undef, d)
    Ls = Vector{}(undef, d)
    
    weightsRetrieved = 0

    for layer = 1:d
        szW = ( dims[layer+1], dims[layer] )
        nW = prod(szW)
        firstIdx = weightsRetrieved + 1
        lastIdx = weightsRetrieved + nW
        W = reshape(w[firstIdx:lastIdx], szW)
        Ws[layer] = W
        weightsRetrieved += nW
    end

    t = Ws[1]*X
    Ls[1] = activation.(t)
    for layer = 2:d
        t = Ws[layer]*Ls[layer-1]
        Ls[layer] = activation.(t)
    end

    ypred = Ls[d]
    # f = (0.5/M)*sum((y - ypred).^2)
    f = 0.5*sum((y - ypred).^2)

    if getGradientToo && classify == false
        g = zeros(n)
        hs = Vector{VecOrMat{Float64}}(undef, d)
        Gs = Vector{VecOrMat{Float64}}(undef, d)

        gradientElementsInserted = 0

        t = ypred - y
        
        for layer = d:-1:1
            h = t.*( Ls[layer].*(1 .- Ls[layer]) )
            t = Ws[layer]'*h

            if layer == 1
                G = h*X'
            else
                G = h*Ls[layer-1]'
            end

            hs[layer] = h

            Gs[layer] = G
            nG = prod(size(G))
            lastIdx = n - gradientElementsInserted
            firstIdx = n - gradientElementsInserted - nG + 1
            g[firstIdx:lastIdx] = vec(G)

            gradientElementsInserted += nG
        end

        return f, g

    elseif !getGradientToo || classify == false
        return f
    
    elseif classify == true
        return ypred

    else
        @error "floc"
    end

end

classify = true
classify = false
# describe(df)

df = preprocessLiverData()

hiddenNodes = [10, 10]
# hiddenNodes = [5, 5]
# hiddenNodes = [8, 8]
# hiddenNodes = [8, 5]
# hiddenNodes = [8, 8, 5]
# hiddenNodes = [10, 10]

# dims = [10, 10, 10, 1]

# MLJ.schema(df)
# vscodedisplay(df)

df_train, df_test = MLJ.partition(df, 0.7, rng=123);
yTrain, XTrain = MLJ.unpack(df_train, ==(:Diagnosis));
trainData = Matrix(XTrain)
classData = Vector(yTrain)

nf = size(trainData, 2)
dims = vcat(nf, hiddenNodes, 1)

w0 = construct_w0_from_dims(dims)

yTest, XTest = MLJ.unpack(df_test, ==(:Diagnosis));
testData = Matrix(XTest);
testClassData = Vector(yTest);

params = Dict(:trainData=>trainData, :classData=>classData, :dims=>dims, 
    :classify=>classify, :testData=>testData, :testClassData=>testClassData)

objective = nnloss;

pr = generate_pr(objective, w0, params=params)

# f = nnloss(pr.x0, pr.p)