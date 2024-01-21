using Plots

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
        # @error "not implemented"
    else
        @error "floc"
    end

    return sampledSpace
end

n = 2
p = 62
# p = 5
method = "Random"
method = "Halton"
method = "Latin Hypercube"
discard = 20
seed = 1234
ss = sampleSpace(n, p, method=method, seed=seed, discard=discard)
# ss = sampleSpaceWithHalton(n, p, discard=discard)

plotSampledSpace(ss, method, view=true, savePlot=true)