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

    p = p.params
    trainData = p[:trainData]
    if size(trainData, 2) == 1
        X = reshape(trainData, length(trainData), 1)
    else
        X = trainData'
    end

    classData = p[:classData]

    if size(classData) == ()
        y = fill(classData, 1, 1)
    else
        y = reshape(classData, 1, length(classData))
    end
    
    dims = p[:dims]
    classify = p[:classify]
    
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
        # @show size(W)
        Ws[layer] = W
        weightsRetrieved += nW
    end

    t = Ws[1]*X
    Ls[1] = activation.(t)
    # @show size(Ls[1])
    for layer = 2:d
        t = Ws[layer]*Ls[layer-1]
        Ls[layer] = activation.(t)
        # @show size(Ls[layer])
    end

    # ypred = vec(Ls[d])
    # @show size(y)
    ypred = Ls[d]
    f = (0.5/M)*sum((y - ypred).^2)

    if getGradientToo && classify == false
        g = zeros(n)
        hs = Vector{VecOrMat{Float64}}(undef, d)
        Gs = Vector{VecOrMat{Float64}}(undef, d)

        gradientElementsInserted = 0

        t = ypred - y
        # @show size(t)
        
        for layer = d:-1:1
            h = t.*( Ls[layer].*(1 .- Ls[layer]) )
            # @show size(h)
            t = Ws[layer]'*h
            # @show size(t)

            if layer == 1
                G = h*X'
            else
                G = h*Ls[layer-1]'
            end

            # @show size(G)
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

dims = [10, 10, 10, 1]
classify = true
classify = false


# describe(df)


w0 = construct_w0_from_dims(dims)
df = preprocessLiverData()

# MLJ.schema(df)

# vscodedisplay(df)

df_train, df_test = MLJ.partition(df, 0.7, rng=123);
yTrain, XTrain = MLJ.unpack(df_train, ==(:Diagnosis));
trainData = Matrix(XTrain)
classData = Vector(yTrain)

# trainData = Matrix(XTrain)[1, :]
# classData = Vector(yTrain)[1]

params = Dict(:trainData=>trainData, :classData=>classData, :dims=>dims, :classify=>classify)

objective = nnloss;

pr = generate_pr(objective, w0, params=params)

f = nnloss(pr.x0, pr.p)