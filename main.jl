include("setup.jl")

println("You are currently using $(Threads.nthreads()) threads.")
println("Your machine has a total of $(Sys.CPU_THREADS) available threads.")


# functionName = "sphere"
# functionName = "dampedSHM";
# functionName = "TestFunction1";
# functionName = "TestFunction2";
# functionName = "TestFunction3";
# functionName = "rosenbrock";

pr = generate_pr(functionName);

verbose = false
# verbose = true;
logging = true;
profiling = false;
benchmarking = false;

@time res = optimize(pr, verbose=verbose)

showresults(res, pr=pr)
# plotresults(pr, res)

# linesearchSW(pr, pr.x0, findDirection(pr, pr.objective(pr.x0, pr.p)[2]))