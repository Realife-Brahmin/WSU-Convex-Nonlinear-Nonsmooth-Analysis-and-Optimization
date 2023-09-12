# main.jl
using BenchmarkTools
using CSV
using DataFrames
using LinearAlgebra
using Revise
using Symbolics

include("src/utilities.jl")
include("src/initializer.jl");
include("src/objective.jl");

rawDataFolder = "rawData/";
filename = rawDataFolder*"FFD.csv";
df = CSV.File(filename) |> DataFrame;
rename!(df, [:t, :V]);

f, ∇f, fnum, ∇fnum, x = objFun(df);

x0 = estimate_x0(df, x)
# fnum([A₀₀, 1, 1, 1, 1, 1])
fnum(x0)