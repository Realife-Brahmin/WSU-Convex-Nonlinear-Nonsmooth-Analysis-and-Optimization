using Pkg
# Pkg.activate(".")
Pkg.instantiate()

using Base.Threads
using BenchmarkTools
using CSV
using DataFrames
using LaTeXStrings
using LinearAlgebra
using Parameters
using Plots
using Printf
using Revise

include("../alg.jl") # include alg.jl from root directory

include("objfuns/objective.jl");
include("display.jl");
include("helperFunctions.jl");

initializeLogFile()

include("findDirection.jl");

include("optimize.jl");
include("optimize2.jl");

include("plotting/plotresults.jl");
include("types.jl");
include("utilities.jl"); 