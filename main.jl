include("setup.jl")

# global const JULIA_NUM_THREADS = 4;
println("You are currently using $(Threads.nthreads()) threads.")
println("Your machine has a total of $(Sys.CPU_THREADS) available threads.")

function quadratic1d(x::Vector{Float64}, p::FuncParam;
    getGradientToo::Bool=true,
    verbose::Bool=false)
    f = x[1]^2
    if getGradientToo
        g = 2x
        return f, g
    else
        return f
    end
end

# function rosenbrock2d(x::Vector{Float64}, p::FuncParam;
#     getGradientToo::Bool=true,
#     verbose::Bool=false)
#     f = (1 - x[1])^2 + 100(x[2] - x[1]^2)^2
#     if getGradientToo
#         g = zeros(Float64, 2)
#         g[1] = 2(1-x[1]) - 400*x[1]*(x[2] - x[1]^2)
#         g[2] = 200(x[2] - x[1]^2)
#         return f, g
#     else
#         return f
#     end
# end

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

functionName = "sphere"
# functionName = "rosenbrock2d"
# functionName = "quadratic1d";
# functionName = "dampedSHM";
# functionName = "TestFunction1";
# functionName = "TestFunction2";
# functionName = "TestFunction3";
# functionName = "rosenbrock";

pr = generate_pr(functionName);


# @error "Okay stop now."


# verbose = false
verbose = true;
logging = true;
profiling = false;
benchmarking = false;

@time res = optimize(pr, verbose=verbose)

showresults(res, pr=pr)
# plotresults(pr, res)



# ProfileView.view();

# Open a file in write mode
# f = open("./logging/profile_results.txt", "w")

# Redirect the output of Profile.print to the file
# Profile.print(f, mincount=5000)

# Close the file
# close(f)