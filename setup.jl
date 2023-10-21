using Pkg
Pkg.activate(".")

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

include("alg.jl") # include alg.jl from root directory

include("src/dampedSHM.jl");
include("src/display.jl");
include("src/drag.jl");
include("src/helperFunctions.jl");
include("src/findDirection.jl");
include("src/LineSearchAlgos.jl");

include("src/objective.jl");
include("src/optimize.jl");
include("src/TestFunctions.jl");
include("src/TestFunctions_own.jl");
include("src/utilities.jl"); 