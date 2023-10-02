include("setup.jl")

# global const JULIA_NUM_THREADS = 4;
println("You are currently using $(Threads.nthreads()) threads.")
println("Your machine has a total of $(Sys.CPU_THREADS) available threads.")

# functionName = "dampedSHM";
functionName = "TestFunction1";
# functionName = "TestFunction2";
# functionName = "TestFunction3";
# functionName = "rosenbrock";

pr = generate_pr(functionName);

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