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
    X = p[:trainData]' # note the transpose
    y = p[:classData]
    
    dims = p[:dims]
    classify = p[:classify]
    
    nf, M = size(X)
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

    Ws = Vector{VecOrMat{Float64}}(undef, d)
    Ls = Vector{VecOrMat{Float64}}(undef, d)
    
    weightsRetrieved = 0
    L = X

    for layer = 1:d
        szW = ( dims[layer], dims[layer+1] )
        nW = prod(szW)
        firstIdx = weightsRetrieved + 1
        lastIdx = weightsRetrieved + nW
        W = reshape(w[firstIdx:lastIdx], szW)
        @show size(W)
        Ws[layer] = W
        L = activation.(W'*L)
        Ls[layer] = L
        @show size(L)
        weightsRetrieved += nW
    end

    ypred = vec(Ls[d])
    f = (0.5/M)*sum((y - ypred).^2)

    if getGradientToo && classify == false
        g = zeros(n)
        hs = Vector{VecOrMat{Float64}}(undef, d)
        
        for layer = d:-1:1
            L = Ls[layer]
            @show size(L)
            @show size(y)
            @show size(h)
            h = (L-y).*L.*(1 .- L)
            if layer == 1
                # Ld = ypred
                # h = (Ld-y).*Ld.*(1 .- Ld)
                G = h*X'
            else
                Wp1 = Ws[layer+1]
                hp1 = hs[layer+1]
                h = Wp1'*hp1*L*(I-L)
                Lm1 = Ls[layer-1]
                G = h*Lm1
            end

            hs[layer] = h

            Gs[layer] = G
            nG = prod(size(G))
            lastIdx = n - gradientElementsInserted
            firstIdx = n - gradientElementsInserted - nG + 1
            g[firstIdx:lastIdx] = vec(G)
        end

        Gs = Vector{VecOrMat{Float64}}(undef, d)

        gradientElementsInserted = 0

        for layer = d:-1:1
            h = hs[layer]
            if layer == 1
                G = h*X'
            else
                Lm1 = Ls[layer-1]
                G = h*Lm1
            end

            # Gs[layer] = G
            # nG = prod(size(G))
            # lastIdx = n - gradientElementsInserted
            # firstIdx = n - gradientElementsInserted - nG + 1
            # g[firstIdx:lastIdx] = vec(G)
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
# classify = false


# describe(df)


w0 = construct_w0_from_dims(dims)
df = preprocessLiverData()

# MLJ.schema(df)

# vscodedisplay(df)

df_train, df_test = MLJ.partition(df, 0.7, rng=123);
yTrain, XTrain = MLJ.unpack(df_train, ==(:Diagnosis));
trainData = Matrix(XTrain)
classData = Vector(yTrain)

params = Dict(:trainData=>trainData, :classData=>classData, :dims=>dims, :classify=>classify)

objective = nnloss;

pr = generate_pr(objective, w0, params=params)

f = nnloss(pr.x0, pr.p)