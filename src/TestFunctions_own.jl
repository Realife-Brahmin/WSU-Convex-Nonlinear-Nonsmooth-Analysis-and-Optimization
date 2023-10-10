FuncParam = NamedTuple{(:params, :data), Tuple{Vector{Float64}, Matrix{Float64}}}

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

function Rastrigin2d(x::Vector{Float64}, p::FuncParam;
    getGradientToo::Bool=true,
    verbose::Bool=false)

    x1, x2 = x[1], x[2]
    f = 10 * 2 + (x1^2 - 10cos(2π*x1)) + (x2^2 - 10cos(2π*x2))
    if getGradientToo
        g = zeros(Float64, 2)
        g[1] = 2x1 + 20π*sin(2π*x1)
        g[2] = 2x2 + 20π*sin(2π*x2)
        return f, g
    else
        return f
    end
end

function rosenbrock2d_oscillatory(x::Vector{Float64}, p::FuncParam; getGradientToo::Bool=true, verbose::Bool=false)
    # Ensure x is of length 2
    if length(x) != 2
        throw(ArgumentError("Input vector x should be of length 2 for the 2D Rosenbrock function."))
    end
    
    # Extract components
    x1, x2 = x
    
    # Compute function value with oscillatory term
    f = (1 - x1)^2 + 100*(x2 - x1^2)^2 + sin(5π*x1)^2

    # If gradient is not needed, return only the function value
    if !getGradientToo
        return f
    end
    
    # Compute gradient with derivative of the oscillatory term
    dfdx1 = -2*(1 - x1) - 400*x1*(x2 - x1^2) + 2sin(5π*x1)*5π*cos(5π*x1)
    dfdx2 = 200*(x2 - x1^2)
    
    g = [dfdx1, dfdx2]

    # Optionally print the computed values
    if verbose
        println("Function value: $f")
        println("Gradient: $g")
    end

    return f, g
end

function rosenbrock2d_oscillatory(x::Vector{Float64}, p::FuncParam; getGradientToo::Bool=true, verbose::Bool=false)
    # Ensure x is of length 2
    if length(x) != 2
        throw(ArgumentError("Input vector x should be of length 2 for the 2D Rosenbrock function."))
    end
    
    # Extract components
    x1, x2 = x
    
    # Compute function value with oscillatory term
    f = (1 - x1)^2 + 100*(x2 - x1^2)^2 + sin(5π*x1)^2

    # If gradient is not needed, return only the function value
    if !getGradientToo
        return f
    end
    
    # Compute gradient with derivative of the oscillatory term
    dfdx1 = -2*(1 - x1) - 400*x1*(x2 - x1^2) + 2sin(5π*x1)*5π*cos(5π*x1)
    dfdx2 = 200*(x2 - x1^2)
    
    g = [dfdx1, dfdx2]

    # Optionally print the computed values
    if verbose
        println("Function value: $f")
        println("Gradient: $g")
    end

    return f, g
end

function hardFunction1(x::Vector{Float64}, p::FuncParam; getGradientToo::Bool=true, verbose::Bool=false)
    # Ensure x is of length 2
    if length(x) != 2
        throw(ArgumentError("Input vector x should be of length 2 for the 2D Rosenbrock function."))
    end
    
    # Extract components
    x1, x2 = x
    
    # Compute function value with oscillatory term
    f = (1.5-x1+x1*x2)^2 + (2.25-x1+x1*x2^2)^2 + (2.625-x1+x1*x2^3)^2

    # If gradient is not needed, return only the function value
    if !getGradientToo
        return f
    end
    
    # Compute gradient with derivative of the oscillatory term
    dfdx1 = 2*(1.5-x1+x1*x2)*(-1+x2) + 2*(2.25-x1+x1*x2^2)*(-1+x2^2)+ 2*(2.625-x1+x1*x2^3)*(-1+x2^3)
    dfdx2 = 2*(1.5-x1+x1*x2)*(x1) + 2*(2.25-x1+x1*x2^2)*(2*x1*x2) + 2*(2.625-x1+x1*x2^3)*(3*x1*x2^2)
    
    g = [dfdx1, dfdx2]

    # Optionally print the computed values
    if verbose
        println("Function value: $f")
        println("Gradient: $g")
    end

    return f, g
end
