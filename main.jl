# main.jl
using BenchmarkTools
using CSV
using DataFrames
using Latexify
using LaTeXStrings
using LinearAlgebra
using Plots
using Revise
using Statistics # for estimating x0, I'm too lazy to code it myself
using Symbolics

include("src/initializer.jl");
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
        maxiter = 200,
        ngtol = 1e-8,
        dxtol = 1e-8,
        lambda = 1,
        lambdaMax = 100,
        linesearch = "Armijo",
        c1 = 1e-4, # Pg 33 (3.1 Step Length)
        c2 = 0.9,
        progress = 10);

functionName = "dampedSHM"

pr = (objective=functionName, x0=x0, alg=alg)

# pr = (alg=alg);
pr == alg # returns true
method = pr.alg.method
# f, ∇f, fnum, ∇fnum, x = objFun(df);
# x0_mine = estimate_x0(df, x)

x0 = [13.8, 8.3, 0.022, 1800, 900, 4.2]
# fnum([x0[1], 1, 1, 1, 1, 1])
# fnum(x0_mine) # should be as close to 0.00 as possible
# fnum(x0_Good)
# lol it got worse after I inserted a more sensible value of A?

t0 = df.t[1]
f, g = dampedSHM(x0, t0)
pk = findDirection(pr, g)
α = linesearch(pr, x0, t0, pk, verbose=true)





