include("setup.jl")

println("You are currently using $(Threads.nthreads()) threads.")
println("Your machine has a total of $(Sys.CPU_THREADS) available threads.")


# functionName = "sphere"
# functionName = "dampedSHM";
# functionName = "Rastrigin2d";
# functionName = "TestFunction1";
# functionName = "TestFunction2";
# functionName = "TestFunction3";
# functionName = "rosenbrock";
functionName = "hardFunction1"
# functionName = "rosenbrock2d_oscillatory"

pr = generate_pr(functionName);

verbose = false
# verbose = true;
logging = true;
profiling = false;
benchmarking = false;

@time res = optimize(pr, verbose=verbose)

showresults(res, pr=pr)

# Generate data for plotting
x1_range = -2:0.05:2
x2_range = -2:0.05:2

f_values = [hardFunction1([x1, x2], empty_FuncParam(), getGradientToo=false) for x1 in x1_range, x2 in x2_range]
# f_values = [Rastrigin2d([x1, x2], empty_FuncParam(), getGradientToo=false) for x1 in x1_range, x2 in x2_range]
# f_values = [sphere([x1, x2], empty_FuncParam(), getGradientToo=false) for x1 in x1_range, x2 in x2_range]

# Plot
contour(x1_range, x2_range, f_values, legend=true, xlabel="x1", ylabel="x2")
# contour(x1_range, x2_range, f_values, legend=true, xlabel="x1", ylabel="x2", title="Rastrigin2d")

# surface(x1_range, x2_range, f_values, xlabel="x1", ylabel="x2", zlabel="f(x1,x2)", title="3D Plot of rosenbrock2d_oscillatory")

# contour(x1_range, x2_range, (x1, x2) -> x1^2 + x2^2)
# plotresults(pr, res)

# linesearchSW(pr, pr.x0, findDirection(pr, pr.objective(pr.x0, pr.p)[2]))