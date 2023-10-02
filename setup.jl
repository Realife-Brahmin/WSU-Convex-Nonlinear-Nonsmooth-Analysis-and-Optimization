using Pkg
Pkg.activate(".")

using Base.Threads
using BenchmarkTools
using CSV
using DataFrames
using Latexify
using LaTeXStrings
using LinearAlgebra
using Optim
using Parameters
using Plots
using Printf
using Profile
using ProfileView
using Revise

include("src/display.jl");
include("src/helperFunctions.jl");
include("src/objective.jl");
include("src/objectiveParallel.jl");
include("src/optimize.jl");
include("src/optimizeParallel.jl");
include("src/TestFunctions.jl");
include("src/utilities.jl"); 