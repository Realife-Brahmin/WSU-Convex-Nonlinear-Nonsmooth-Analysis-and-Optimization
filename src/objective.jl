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
    objective = eval(Symbol(functionName))
    println("We'll be working with the $(functionName) function.")

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
        x0 = Float64.([10, 12])
        params = Float64[]

    elseif functionName == "hardFunction1"
        x0 = Float64.([-2, -2])
        params = Float64[]

    elseif functionName == "drag"
        n = 1000
        x0 = Float64.(collect(LinRange(0.0, 1.0, n+2)[2:n+1]))
        params = Float64[]
    else
        @error "Unknown Function. If the function definition is known, please define data, params, x0 first!"
    end

    p = (params=params, data=data)

    println("Problem (pr::NamedTuple) generated")
    return (p=p, x0=x0, objective=objective, alg=alg)
end



# end