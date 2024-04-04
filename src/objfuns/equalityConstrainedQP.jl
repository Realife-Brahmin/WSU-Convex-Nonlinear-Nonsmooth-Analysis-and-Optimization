function equalityConstrainedQP(w::Vector{Float64},
    pDict;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=false) # yes, false by default

    @unpack G, A, c, b = pDict
    f = 1//2*transpose(w)*G*w + transpose(w)*c

    if getGradientToo
        g = G*w + c
        return f, g
    else
        return f
    end
    
end