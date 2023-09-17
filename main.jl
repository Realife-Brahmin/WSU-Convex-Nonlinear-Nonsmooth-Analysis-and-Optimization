# main.jl
using BenchmarkTools
using CSV
using DataFrames
using Latexify
using LaTeXStrings
using LinearAlgebra
using Plots
using Profile
using ProfileView
using Revise
using Symbolics

# include("src/initializer.jl");
include("src/helperFunctions.jl");
include("src/objective.jl");
include("src/plotter.jl");
include("src/utilities.jl");

rawDataFolder = "rawData/";
filename = rawDataFolder*"FFD.csv";
df = CSV.File(filename) |> DataFrame;
rename!(df, [:t, :V]);

logging = true

if logging
        if !isdir("./logging")
                println("Creating logging directory since it doesn't exist.")
                mkdir("./logging")
        else 
                println("No need to create logging directory, it already exists.")
                initialize_logging(overwrite=true)
        end
end
scatter_voltage_vs_time(df)

alg = (method = "GradientDescent",
        maxiter = 200,
        ngtol = 1e-10,
        dftol = 1e-12,
        dxtol = 1e-10,
        lambda = 1,
        lambdaMax = 100,
        # linesearch = "Armijo",
        linesearch = "StrongWolfe",
        c1 = 1e-4, # Pg 33 (3.1 Step Length)
        c2 = 0.9,
        progress = 50);

functionName = "dampedSHM";

x0 = [13.8, 8.3, 0.022, 1800, 900, 4.2];

pr = (objective=functionName, x0=x0, alg=alg, df=df);


dftol = pr.alg.dftol;
progress = pr.alg.progress;
maxiter = pr.alg.maxiter;

fnext = 1e10;
fₖ = computeCost(pr, x0, getGradientToo=false);
x = pr.x0;
n = length(x);
itr = 1;
fvals, αvals = [zeros(Float64, maxiter) for _ in 1:2];
backtrackVals = zeros(Int64, maxiter, 1);
xVals = zeros(Float64, n, maxiter);

println("Begin with the solver:")
@profile begin
        while abs(fnext-fₖ) ≥ dftol && itr ≤ maxiter
                global fₖ, x, fnext, itr 
                printOrNot = (itr%progress==0)
                myprintln(printOrNot, "Iteration $(itr):", log=true)
                fₖ, ∇fₖ = computeCost(pr, x)
                myprintln(printOrNot, fₖ, log=true)
                pₖ = findDirection(pr, ∇fₖ)
                α, x, fnext, backtrackNum = linesearch(pr, x, pₖ, verbose=true)
                fvals[itr] = fnext
                αvals[itr] = α
                backtrackVals[itr] = backtrackNum
                xVals[:, itr] = x
                itr += 1
        end
end
if itr > maxiter
        @warn ("Failed to converge despite $(maxiter) iterations!")
else
        println("Convergence achieved in $(itr) iterations 😄")
        fvals = fvals[1:itr]
        αvals = αvals[1:itr]
        backtrackVals = backtrackVals[1:itr]
end

# For testing linesearch
# fₖ, ∇fₖ = computeCost(pr, x0)
# pₖ = findDirection(pr, ∇fₖ)
# linesearch(pr, x0, pₖ, verbose=true)

ProfileView.view();
# pₖ = findDirection(pr, g)
# α = linesearch(pr, x0, pₖ, verbose=true)