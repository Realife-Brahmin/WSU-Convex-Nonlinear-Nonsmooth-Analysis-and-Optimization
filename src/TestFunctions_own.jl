function rosenbrock2d(x::Vector{Float64}, p::FuncParam; getGradientToo::Bool=true, verbose::Bool=false)
    # Ensure x is of length 2
    if length(x) != 2
        throw(ArgumentError("Input vector x should be of length 2 for the 2D Rosenbrock function."))
    end
    
    # Extract components
    x1, x2 = x
    
    # Compute function value
    f = (1 - x1)^2 + 100*(x2 - x1^2)^2

    # If gradient is not needed, return only the function value
    if !getGradientToo
        return f
    end
    
    # Compute gradient
    dfdx1 = -2*(1 - x1) - 400*x1*(x2 - x1^2)
    dfdx2 = 200*(x2 - x1^2)
    
    g = [dfdx1, dfdx2]

    # Optionally print the computed values
    if verbose
        println("Function value: $f")
        println("Gradient: $g")
    end

    return f, g
end

function sphere(x::Vector{Float64}, p::FuncParam;
    getGradientToo::Bool=true,
    verbose::Bool=false)

    f = dot(x, x)
    if getGradientToo
        g = 2x
        return f, g
    else
        return f
    end
end