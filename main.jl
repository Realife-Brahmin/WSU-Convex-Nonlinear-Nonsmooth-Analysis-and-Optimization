include("setup.jl")

println("You are currently using $(Threads.nthreads()) threads.")
println("Your machine has a total of $(Sys.CPU_THREADS) available threads.")

# functionName = "dampedSHM";
functionName = "drag"
# functionName = "rosenbrock";
# functionName = "sphere";
# functionName = "TestFunction1";
# functionName = "TestFunction2";
# functionName = "TestFunction3";

pr = include("src/objfuns/"*String(functionName)*".jl")
# verbose = false
verbose = true;
verbose_ls = false;
# verbose_ls = true;
verbose_ls = verbose & verbose_ls
logging = true;
profiling = false;
benchmarking = false;

@time res = optimize(pr, verbose=verbose, verbose_ls=verbose_ls)

showresults(res)