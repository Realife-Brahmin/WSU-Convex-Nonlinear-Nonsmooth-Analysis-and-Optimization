# main.jl
using Pkg
Pkg.activate(".")

using Base.Threads
using BenchmarkTools
using CSV
using DataFrames
using Latexify
using LaTeXStrings
using LinearAlgebra
using Parameters
using Plots
using Printf
using Profile
using ProfileView
using Revise
using Symbolics

include("src/display.jl");
include("src/helperFunctions.jl");
include("src/objective.jl");
include("src/objectiveParallel.jl");
include("src/optimize.jl");
include("src/optimizeParallel.jl");
include("src/TestFunctions.jl");
include("src/utilities.jl"); 

# global const JULIA_NUM_THREADS = 4;

functionName = "dampedSHM"
# functionName = "dampedSHM_Parallel"
# functionName = "TestFunction1"

if functionName == "dampedSHM" || functionName == "dampedSHM_Parallel"
        rawDataFolder = "rawData/";
        filename = rawDataFolder*"FFD.csv";
        df = CSV.File(filename) |> DataFrame;
        rename!(df, [:x, :y]);
        data = Matrix(df)

        params = []
        x0 = [13.8, 8.3, 0.022, 1800, 900, 4.2];

        scatter_voltage_vs_time(df)
elseif functionName == "TestFunction1"
        data = []
        params = [2]
        x0 = randn(2)
elseif functionName == "TestFunction2"
        data = []
        params = [1, -16, 5]
        x0 = -2 .+ 2 .* rand(15)
elseif functionName == "TestFunction3"
        data = []
        params == [10]
        x0 = sort!(rand(10).^2, rev=true)
elseif functionName = "rosenbrock"
        data = []
        params == [10]
        x0 = 0.1:0.1:10
else
        @error "Unkown Function. If the function definition is known, 
        please define data, params, x0 first!"
end

verbose = false
# verbose = true;
logging = true;

if logging
        if !isdir("./logging")
                println("Creating logging directory since it doesn't exist.")
                mkdir("./logging")
        else 
                println("No need to create logging directory, it already exists.")
                initialize_logging(overwrite=true)
        end
end

alg = (method = "GradientDescent",
        maxiter = Int(1e5),
        ngtol = 1e-10,
        dftol = 1e-12,
        dxtol = 1e-10,
        lambda = 1,
        lambdaMax = 100,
        linesearch = "Armijo",
        # linesearch = "StrongWolfe",
        c1 = 1e-4, # Pg 33 (3.1 Step Length)
        c2 = 0.9,
        progress = 50);

if isdefined(Main, :obj)
        println("The obj function of name $(nameof(obj)) is already defined.")
else
        const obj = eval(Symbol(functionName))
        println("We'll be working with the $(nameof(obj)) function.")
end

p = (data=data, params=params)

pr = (alg=alg, objective=obj, p=p, x0=x0);

println("You are currently using $(Threads.nthreads()) threads.")
println("Your machine has a total of $(Sys.CPU_THREADS) available threads.")

# @profile begin
# @btime begin
@time begin
        res = optimize(pr, verbose=verbose, itrStart=7)
end

showresults(res)
# plotresults(pr, res)


# ProfileView.view();

# Open a file in write mode
# f = open("./logging/profile_results.txt", "w")

# Redirect the output of Profile.print to the file
# Profile.print(f, mincount=5000)

# Close the file
# close(f)