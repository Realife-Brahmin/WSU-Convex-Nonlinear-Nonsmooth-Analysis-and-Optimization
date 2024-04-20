include("objective.jl")

function QPObjectiveFunction(x::DenseVector,
    pDict;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)

    n = length(x)
    @unpack G, c = pDict
    if haskey(pDict, :c0)
        c0 = pDict[:c0]
    else
        c0 = 0.0
    end

    f = 1 // 2 * transpose(x) * G * x + transpose(x) * c + c0

    if !getGradientToo
        return f
    elseif getGradientToo
        g = G*x + c
        return f, g
    else
        @error "floc"
    end

end