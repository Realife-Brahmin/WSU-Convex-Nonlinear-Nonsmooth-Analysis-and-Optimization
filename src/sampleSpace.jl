include("./generateHaltonSequence.jl")
include("./generateRandomSequence.jl")
include("./generateLatinHypercubeSequence.jl")
include("./plotSampledSpace.jl")

function sampleSpace(n, p; 
    method::String="Halton",
    discard::Int=0,
    seed::Int=1)

    if method == "Halton"
        sampledSpace = sampleSpaceHalton(n, p, discard=discard)
    elseif method == "Latin Hypercube"
        sampledSpace = sampleSpaceLatinHypercube(n, p, seed=seed)
    elseif method == "Random"
        sampledSpace = sampleSpaceRandom(n, p, seed=seed)
    else
        @error "floc"
    end

    return sampledSpace
end