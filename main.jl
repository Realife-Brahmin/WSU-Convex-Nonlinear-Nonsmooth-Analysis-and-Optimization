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
include("src/display.jl");
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
        maxiter = 10000,
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

functionName = "dampedSHM"

if isdefined(Main, :obj)
        println("The obj function of name $(nameof(obj)) is already defined.")
else
        const obj = eval(Symbol(functionName))
        println("We'll be working with the $(nameof(obj)) function.")
end
        

x0 = [13.8, 8.3, 0.022, 1800, 900, 4.2];

pr = (objective=obj, x0=x0, alg=alg, df=df);


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
xvals = zeros(Float64, n, maxiter);

println("Begin with the solver:")
@profile begin
        while abs(fnext-fₖ) ≥ dftol && itr ≤ maxiter
                global fₖ, x, fnext, itr 
                printOrNot = (itr%progress==0)
                myprintln(printOrNot, "Iteration $(itr):", log=true)
                fₖ, ∇fₖ = computeCost(pr, x)
                pₖ = findDirection(pr, ∇fₖ)
                α, x, fnext, backtrackNum = linesearch(pr, x, pₖ, verbose=printOrNot, itrStart=7)
                fvals[itr] = fnext
                αvals[itr] = α
                backtrackVals[itr] = backtrackNum
                xvals[:, itr] = x
                itr += 1
        end
end
if itr > maxiter
        converged = false
        statusMessage = "Failed to converge despite $(maxiter) iterations! 😢"
        @warn statusMessage
else
        converged = true
        statusMessage = "Convergence achieved in $(itr) iterations 😄"
        println(statusMessage)
        # fvals = fvals[1:itr]
        # αvals = αvals[1:itr]
        # backtrackVals = backtrackVals[1:itr]
        # xvals = xvals[1:itr]
        # truncating arrays as they weren't filled to capacity
        fvals, αvals, backtrackVals, xvals = [arr[1:itr] for arr in (fvals, αvals, backtrackVals, xvals)]

end

res = (converged=converged, statusMessage=statusMessage, fvals=fvals, αvals=αvals, backtrackVals=backtrackVals, xvals=xvals)

showresults(res)

# For testing linesearch
# fₖ, ∇fₖ = computeCost(pr, x0);
# pₖ = findDirection(pr, ∇fₖ);
# linesearch(pr, x0, pₖ, verbose=true, itrStart=7);

ProfileView.view();
# pₖ = findDirection(pr, g)
# α = linesearch(pr, x0, pₖ, verbose=true)