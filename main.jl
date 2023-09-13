# main.jl
using BenchmarkTools
using CSV
using DataFrames
using Latexify
using LaTeXStrings
using LinearAlgebra
using Plots
using Revise
using Symbolics

include("src/utilities.jl");
include("src/initializer.jl");
include("src/objective.jl");
include("src/plotter.jl");

rawDataFolder = "rawData/";
filename = rawDataFolder*"FFD.csv";
df = CSV.File(filename) |> DataFrame;
rename!(df, [:t, :V]);

scatter_voltage_vs_time(df)

@btime f, ∇f, fnum, ∇fnum, x = objFun(df);

x0 = estimate_x0(df, x)
# fnum([x0[1], 1, 1, 1, 1, 1])
fnum(x0) # should be close to A₀ ≈ 13.76
# lol it got worse after I inserted a more sensible value of A?






