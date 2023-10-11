include("setup.jl")

println("You are currently using $(Threads.nthreads()) threads.")
println("Your machine has a total of $(Sys.CPU_THREADS) available threads.")


# functionName = "sphere"
functionName = "dampedSHM";
# functionName = "Rastrigin2d";
# functionName = "TestFunction1";
# functionName = "TestFunction2";
# functionName = "TestFunction3";
# functionName = "rosenbrock";
# functionName = "hardFunction1"
# functionName = "rosenbrock2d_oscillatory"

pr = generate_pr(functionName);


verbose = false
# verbose = true;
logging = true;
profiling = false;
benchmarking = false;

f0, ∇f0 = pr.objective(pr.x0, pr.p)
pk = findDirection(pr, pr.x0)
α0 = 1
ϕ0 = f0
ϕprime0 = ∇f0

strongWolfe(pr.objective, pr.x0, pr.p, pk, α0, ϕ0, ϕprime0; pr=pr)

# @time res = optimize(pr, verbose=verbose)

# showresults(res, pr=pr)



# linesearchSW(pr, pr.x0, findDirection(pr, pr.objective(pr.x0, pr.p)[2]))