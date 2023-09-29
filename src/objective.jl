# module objective

using Base.Threads
using DataFrames
using Symbolics

include("helperFunctions.jl");
include("../alg.jl") # include alg.jl from parent directory

function generate_pr(functionName::String)
    data = Matrix{Float64}(undef, 0, 0)
    params = Vector{Float64}()
    x0 = Vector{Float64}()
    
    # Get the function reference from functionName
    if isdefined(Main, :objective)
        println("The objective function of name $(functionName) is already defined.")
    else
        objective = eval(Symbol(functionName))
        println("We'll be working with the $(functionName) function.")
    end

    if functionName == "dampedSHM" || functionName == "dampedSHM_Parallel"
        rawDataFolder = "rawData/"
        filename = rawDataFolder * "FFD.csv"
        df = CSV.File(filename) |> DataFrame
        rename!(df, [:x, :y])
        data = Matrix(df)
        x0 = [13.8, 8.3, 0.022, 1800, 900, 4.2]
        
    elseif functionName == "TestFunction1"
        x0 = randn(2)
        params = [2]

    elseif functionName == "TestFunction2"
        x0 = -2 .+ 2 .* rand(15)
        params = [1, -16, 5]

    elseif functionName == "TestFunction3"
        x0 = sort!(rand(10).^2, rev=true)
        params = [10]

    elseif functionName == "rosenbrock"
        x0 = 0.1:0.1:10
        params = [10]

    else
        @error "Unknown Function. If the function definition is known, please define data, params, x0 first!"
    end

    p = (params=params, data=data)

    println("Problem (pr::NamedTuple) generated")
    return (p=p, x0=x0, objective=objective, alg=alg)
end



function dampedSHM(x::Vector{Float64}, 
    p::NamedTuple{ (:params, :data), Tuple{ Vector{Float64}, Matrix{Float64} } };
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)

    data = p.data
    M = size(data, 1)
    nf = size(data, 2) - 1
    y = data[:, nf+1]
    xf = data[:, 1:nf]
    t = xf

    n = length(x) # note that df's x (time) is different from function parameter x
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