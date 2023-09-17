# main.jl
using BenchmarkTools
using CSV
using DataFrames
using Latexify
using LaTeXStrings
using LinearAlgebra
using Plots
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

scatter_voltage_vs_time(df)

alg = (method = "GradientDescent",
        maxiter = 2000,
        ngtol = 1e-10,
        dftol = 1e-12,
        dxtol = 1e-10,
        lambda = 1,
        lambdaMax = 100,
        linesearch = "Armijo",
        c1 = 1e-4, # Pg 33 (3.1 Step Length)
        c2 = 0.9,
        progress = 50);

functionName = "dampedSHM";

x0 = [13.8, 8.3, 0.022, 1800, 900, 4.2];

pr = (objective=functionName, x0=x0, alg=alg, df=df);


dftol = pr.alg.dftol
progress = pr.alg.progress;
maxiter = pr.alg.maxiter;

fnext = 1e10;
fₖ = computeCost(pr, x0, getGradientToo=false);
x = pr.x0;
n = length(x)
itr = 1
fvals, αvals = [zeros(Float64, maxiter) for _ in 1:2]
backtrackVals = zeros(Int64, maxiter, 1);
xVals = zeros(Float64, n, maxiter);
# @profile begin
        while abs(fnext-fₖ) ≥ pdftol && itr ≤ maxiter 
                printOrNot = (itr%progress==0)
                myprintln(printOrNot, "Iteration $(itr):")
                fₖ, ∇fₖ = computeCost(pr, x)
                myprintln(printOrNot, fₖ)
                pₖ = findDirection(pr, ∇fₖ)
                α, x, fnext, backtracks = linesearch(pr, x, pₖ, verbose=false)
                fvals[itr] = fnext
                αvals[itr] = α
                backtrackVals[itr] = backtracks
                xVals[:, itr] = x
                itr += 1
        end
# end
if itr > maxiter
        @error("Failed to converge despite $(maxiter) iterations!")
else
        fvals = fvals[1:itr]
        αvals = αvals[1:itr]
        backtrackVals = backtrackVals[1:itr]
end

# pₖ = findDirection(pr, g)
# α = linesearch(pr, x0, t0, pₖ, verbose=true)





