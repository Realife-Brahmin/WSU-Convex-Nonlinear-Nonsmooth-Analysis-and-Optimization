include("objective.jl")

function QPObjectiveFunction(x::Vector{Float64},
    pDict;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)

    n = length(x)
    @unpack G, c = pDict
    f = 1 // 2 * transpose(x) * G * x + transpose(x) * c

    if !getGradientToo
        return f
    elseif getGradientToo
        g = G*x + c
        return f, g
    else
        @error "floc"
    end

end