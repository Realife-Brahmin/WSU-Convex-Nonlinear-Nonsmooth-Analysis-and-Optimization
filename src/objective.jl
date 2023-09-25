# module objective

using Base.Threads
using DataFrames
using Symbolics

include("helperFunctions.jl");

function dampedSHM(x::Vector{Float64}, 
    # p::NamedTuple{(:df, :params), Tuple{DataFrame, Vector{Float64}}};
    p::NamedTuple{(:data, :params), Tuple{Matrix{Float64}, Vector{Float64}}};
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)

    # df = p.df
    data = p.data
    M = size(data, 1)
    nf = size(data, 2) - 1
    y = data[:, nf+1]
    xf = data[:, 1:nf]
    t = xf
    # y = df.y
    # t = df.x
    n = length(x) # note that df's x (time) is different from function parameter x
    # M = length(y)
    A₀, A, τ, ω, α, ϕ = x
    f = 0.0;
    g = zeros(Float64, n)
    for k = 1:M
        tₖ = t[k]
        yₖ = y[k]
        expₖ = exp(-tₖ/τ)
        Sₖ = expₖ*sin((ω+α*tₖ)tₖ + ϕ)
        Cₖ = expₖ*cos((ω+α*tₖ)tₖ + ϕ)

        # ŷₖ = A₀+ A*exp(-tₖ/τ)sin((ω+α*tₖ)tₖ + ϕ)
        ŷₖ = A₀ + A*Sₖ
        Δyₖ = ŷₖ - yₖ
        f += (1/M)*Δyₖ^2
        if getGradientToo
            g += (2/M)* Δyₖ * [1, Sₖ, A*tₖ*(τ^-2)*Sₖ, A*tₖ*Cₖ, A*tₖ^2*Cₖ, A*Cₖ]
        end
    end
    
    if getGradientToo
        return f, g
    else
        return f
    end
end

# end