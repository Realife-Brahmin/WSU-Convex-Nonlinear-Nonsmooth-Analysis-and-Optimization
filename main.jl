include("setup.jl")

println("You are currently using $(Threads.nthreads()) threads.")
println("Your machine has a total of $(Sys.CPU_THREADS) available threads.")

verbose = false
# verbose = true;
verbose_ls = false;
# verbose_ls = true;
# verbose_ls = verbose & verbose_ls
logging = true;
profiling = false;
# benchmarking = false;
benchmarking = true;

warmStart = true

# functionName = "dampedSHM";
# functionName = "drag"; functionName == "drag" ? verbose = false : verbose = verbose  
# functionName = "rosenbrock";
# functionName = "sphere";
# functionName = "TestFunction1";
# functionName = "TestFunction2";
# functionName = "TestFunction3";
functionName = "signalDenoise";

pr = include("src/objfuns/"*String(functionName)*".jl")

# Call the function with the initial problem setup
# res = @btime begin
@time begin
    res = warm_start_optimize(pr, verbose=verbose, verbose_ls=verbose_ls)
end

showresults(res)

plotresults(res, savePlot=true)