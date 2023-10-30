function signalDenoise(x::Vector{Float64}, 
    p::FuncParam;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)

    n = length(x)
    f = zeros(Float64, n)
    if getGradientToo
        g = zeros(Float64, n)
        return f, g
    else
        return f
    end

    @error "forbidden loc"

end