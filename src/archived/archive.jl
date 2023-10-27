    # Get the function reference from functionName
    # objective = eval(Symbol(functionName))
    # println("We'll be working with the $(functionName) function.")

    # if functionName == "rosenbrock2d_oscillatory"
    #     x0 = Float64.([0.1, 0.2])
    #     params = Float64[]

    # elseif functionName == "sphere"
    #     x0 = Float64.([0, -100, 11, -23])
    #     params = Float64[]
        
    # elseif functionName == "Rastrigin2d"
    #     x0 = Float64.([10, 12])
    #     params = Float64[]

    # elseif functionName == "hardFunction1"
    #     x0 = Float64.([-2, -2])
    #     params = Float64[]

    # else
    #     @error "Unknown Function. If the function definition is known, please define data, params, x0 first!"
    # end