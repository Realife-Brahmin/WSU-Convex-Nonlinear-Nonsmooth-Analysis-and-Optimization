include("setup.jl")

println("You are currently using $(Threads.nthreads()) threads.")
println("Your machine has a total of $(Sys.CPU_THREADS) available threads.")

# functionName = "dampedSHM";
functionName = "drag"
# functionName = "hardFunction1"
# functionName = "Rastrigin2d";
# functionName = "rosenbrock";
# functionName = "rosenbrock2d_oscillatory"
# functionName = "sphere"
# functionName = "TestFunction1";
# functionName = "TestFunction2";
# functionName = "TestFunction3";

pr = generate_pr(functionName);


verbose = false
# verbose = true;
logging = true;
profiling = false;
benchmarking = false;

# f0, ∇f0 = pr.objective(pr.x0, pr.p)
# pk = findDirection(pr, pr.x0)
# α0 = 1
# ϕ0 = f0
# ϕprime0 = ∇f0

# strongWolfe(pr.objective, pr.x0, pr.p, pk, α0, ϕ0, ϕprime0; pr=pr)

@time res = optimize(pr, verbose=verbose)

showresults(res)



# linesearchSW(pr, pr.x0, findDirection(pr, pr.objective(pr.x0, pr.p)[2]))