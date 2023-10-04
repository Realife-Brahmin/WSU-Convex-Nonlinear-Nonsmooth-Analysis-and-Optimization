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

functionName = "quadratic1d";
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