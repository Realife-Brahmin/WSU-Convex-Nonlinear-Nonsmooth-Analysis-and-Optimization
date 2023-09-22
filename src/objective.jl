# module objective

using Base.Threads
using DataFrames
using Symbolics

include("helperFunctions.jl");

function dampedSHM(x::Vector{Float64}, 
    p::NamedTuple{(:df, :params), Tuple{DataFrame, Vector{Float64}}};
    getGradientToo::Bool=true)

    df = p.df
    y = df.y
    t = df.x
    n = length(x) # note that df's x (time) is different from function parameter x
    M = length(y)
    # params = p.params
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

# function buildFunctions_DampedSHM()
    
#     @variables t, A₀, A, τ, ω, α, ϕ
#     x = [A₀, A, τ, ω, α, ϕ] 
#     f = A₀ + A*exp(-t/τ)sin((ω+α*t)t + ϕ)
#     fnum = build_function(f, x, t, expression=Val{false})
#     ∇f = Symbolics.gradient(f, x)
#     ∇fnum = build_function(∇f, x, t, expression=Val{false})
#     ∇fnum = ∇fnum[1] # Only taking the first function from the tuple

#     function dampedSHM(x::Vector{Float64},
#             p::Union{Float64, Vector{Float64}};
#             getGradientToo::Bool=true,
#             printSymbolicEquations::Bool=true,
#             verbose::Bool=false)
        
#         if printSymbolicEquations
#             println("Function: ", f)
#             if getGradientToo
#                 println("Gradient: ", ∇f)
#             end
#         end
        
#         if isa(p, Float64)
#             println("Okay, p is a single Float.")
#             t = p;
#         elseif isa(p, Vector{Float64})
#             t = p[1]
#         else
#             @error "Bad Condition"
#         end

#         value = fnum(x, t)
#         if getGradientToo
#             gradient = ∇fnum(x, t)
#             return value, gradient
#         else
#             return value
#         end
#     end
    
#     return dampedSHM
# end


# # Create your dampedSHM function with your desired settings
# dampedSHM = buildFunctions_DampedSHM()


function computeCost(pr::NamedTuple, xnow::Vector{Float64}; getGradientToo::Bool=true, verbose::Bool=false, log=true)
    
    df = pr.df
    p = pr.p
    np = length(p)
    nf = size(df, 2) - 1
    M = length(df.y)
    
    if !isempty(df)
        # Calculate SSE
        if getGradientToo
            fval_contributions, gradient_contributions = [ 
                begin
                    fval, gval = pr.objective(xnow, Vector{Float64}(df[i, 1:nf]))
                    -2/M * (df.y[i] - fval) * gval
                end
                for i in 1:M
            ]
            Fval = sum(fval_contributions)
            Gval = sum(gradient_contributions)
            return Fval, Gval
        else
            @show Vector{Float64}(df[1, 1:nf])
            residuals = [df.y[i] - pr.objective(xnow, Vector{Float64}(df[i, 1:nf]), getGradientToo=false) for i in 1:M]
            Fval = mean(residuals.^2)
            return Fval
        end
    else
        @error "Yet to write function for usage without data"
    end

end


# end