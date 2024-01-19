include("./generateHaltonSequence.jl")

function sampleSpace(n, p; 
    method::String="Halton",
    discard::Int=0)

    space = zeros(n, p)
    if method == "Halton"
        sampledSpace = sampleSpaceWithHalton(n, p, discard=discard)
    elseif method == "Latin Hypercube"
        @error "not implemented"
    elseif method == "Random"
        @error "not implemented"
    else
        @error "floc"
    end

    return sampledSpace
end

n = 2
p = 27
discard = 20
ss = sampleSpace(n, p, discard=discard)
# ss = sampleSpaceWithHalton(n, p, discard=discard)
Plots.theme(:dao)
p1 = scatter(ss[1, :], ss[2, :],
        aspect_ratio=:equal,
        xlims=[0.0, 1.0],
        ylims=[0.0, 1.0],
        label=:none);
display(p1)