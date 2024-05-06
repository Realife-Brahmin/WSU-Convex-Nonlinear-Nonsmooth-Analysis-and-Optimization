include("src/setup.jl")

myprintln(true, "You are currently using $(Threads.nthreads()) threads.", color=:yellow)
myprintln(true, "Your machine has a total of $(Sys.CPU_THREADS) available threads.", color=:yellow)
myprintln(true, "*"^25)
# revive this branch from my notebook
verbose = false
verbose = true;
verbose_ls = false;
verbose_ls = true;
# verbose_ls = verbose & verbose_ls
logging = true; 
profiling = false;
# benchmarking = false;
benchmarking = true;

warmStart = true

# functionName = "dampedSHM";
# functionName = "drag"; functionName == "drag" ? verbose = false : verbose = verbose  
# functionName = "fireLocation";
# functionName = "nnloss"; functionName == "nnloss" ? verbose = false : verbose = verbose;
# functionName = "pathtime"
# functionName = "propellor"
# functionName = "receiverLocation"
# functionName = "rosenbrock";
# functionName = "signalDenoise";
# functionName = "sphere";
# functionName = "transmissionLines";
# functionName = "TestFunction1";
# functionName = "TestFunction2";
# functionName = "TestFunction3";

functionName = "alpTestFunction01";
# functionName = "asqpTestFunction1";
# functionName = "ecqpTestFunction1";
# functionName = "ecqpTestFunction2";
# functionName = "Rastrigin";

pr = include("src/objfuns/"*functionName*".jl")

# res = @btime begin
@time begin
    res = tailoredOptimize(pr, verbose=verbose, verbose_ls=verbose_ls)
end

showresults(res)

plotresults(res, savePlot=true)
