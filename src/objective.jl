# module objective

using Base.Threads
using DataFrames
using Symbolics

include("helperFunctions.jl");

function buildFunctions_DampedSHM()

    @variables t, A₀, A, τ, ω, α, ϕ
    x = [A₀, A, τ, ω, α, ϕ] 
    f = A₀ + A*exp(-t/τ)sin((ω+α*t)t + ϕ)
    fnum = build_function(f, x, t, expression=Val{false})
    ∇f = Symbolics.gradient(f, x)
    ∇fnum = build_function(∇f, x, t, expression=Val{false})
    ∇fnum = ∇fnum[1] # Only taking the first function from the tuple

    function dampedSHM(x::Vector{Float64},
            p::Union{Float64, Vector{Float64}};
            getGradientToo::Bool=true,printSymbolicEquations::Bool=false,
            verbose::Bool=false)
        
        t = p;
        if printSymbolicEquations
            println("Function: ", f)
            if getGradientToo
                println("Gradient: ", ∇f)
            end
        end
        
        value = fnum(x, t)
        if getGradientToo
            gradient = ∇fnum(x, t)
            return value, gradient
        else
            return value
        end
    end
    
    return dampedSHM
end


# Create your dampedSHM function with your desired settings
dampedSHM = buildFunctions_DampedSHM()


function computeCost(pr::NamedTuple, xnow::Vector{Float64}; getGradientToo::Bool=true, verbose::Bool=false, log=true)
    
    df = pr.df
    p = pr.p
    np = length(p)
    M = length(df.y)
    
    if !isempty(df)
        # Calculate SSE
        @show Vector{Float64}(df[1, 1:np])
        residuals = [df.y[i] - pr.objective(xnow, Vector{Float64}(df[i, 1:np]), getGradientToo=false) for i in 1:M]
        Fval = mean(residuals.^2)

        if getGradientToo
            gradient_contributions = [ 
                begin
                    fval, gval = pr.objective( xnow, Vector{Float64}(df[i, 1:np]))
                    -2/M * (df.V[i] - fval) * gval
                end
                for i in 1:M
            ]
            Gval = sum(gradient_contributions)
            return Fval, Gval
        else
            return Fval
        end
    else
        @error "Yet to write function for usage without data"
    end

end


# end