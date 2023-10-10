# module objective

using Base.Threads
using DataFrames

include("helperFunctions.jl");
include("../alg.jl") # include alg.jl from parent directory

FuncParam = NamedTuple{(:params, :data), Tuple{Vector{Float64}, Matrix{Float64}}}

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
        # x0 = randn(2)
        x0 = [0.5, 0.5]
        params = Float64.([2])

    elseif functionName == "TestFunction2"
        x0 = -2 .+ 2 .* rand(15)
        params = Float64.([1, -16, 5])

    elseif functionName == "TestFunction3"
        x0 = sort!(rand(10).^2, rev=true)
        params = Float64.([10])

    elseif functionName == "rosenbrock"
        x0 = collect(0.1:0.1:1)
        params = Float64.([10])

    elseif functionName == "rosenbrock2d"
        x0 = Float64.([0.1, 0.2])
        params = Float64[]

    elseif functionName == "rosenbrock2d_oscillatory"
        x0 = Float64.([0.1, 0.2])
        params = Float64[]

    elseif functionName == "sphere"
        x0 = Float64.([0, -100, 11, -23])
        params = Float64[]
        
    elseif functionName == "Rastrigin2d"
        # x0 = Float64.([0.01, 0.01])
        # x0 = Float64.([1, 1])
        x0 = Float64.([10, 12])
        params = Float64[]

    elseif functionName == "hardFunction1"
        x0 = Float64.([-2, -2])
        params = Float64[]
        
    else
        @error "Unknown Function. If the function definition is known, please define data, params, x0 first!"
    end

    p = (params=params, data=data)

    println("Problem (pr::NamedTuple) generated")
    return (p=p, x0=x0, objective=objective, alg=alg)
end

function dampedSHM(x::Vector{Float64}, 
    p::FuncParam;
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