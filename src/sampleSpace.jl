using Plots

include("./generateHaltonSequence.jl")
include("./generateRandomSequence.jl")

function sampleSpace(n, p; 
    method::String="Halton",
    discard::Int=0,
    seed::Int=1)

    if method == "Halton"
        sampledSpace = sampleSpaceHalton(n, p, discard=discard)
    elseif method == "Latin Hypercube"
        @error "not implemented"
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
method = "Random"
method = "Halton"
discard = 20
seed = 1234
ss = sampleSpace(n, p, method=method, seed=seed, discard=discard)
# ss = sampleSpaceWithHalton(n, p, discard=discard)
Plots.theme(:dao)
p1 = scatter(ss[1, :], ss[2, :],
        aspect_ratio=:equal,
        title="Sampled Space using $(method) method",
        xlims=[0.0, 1.0],
        ylims=[0.0, 1.0],
        label=:none);
display(p1)
ext = ".png"
filename = method*"_n_"*string(n)*"_p_"*string(p)*ext
fullAdd = joinpath("processedData", filename)
savefig(p1, fullAdd)