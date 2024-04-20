using Pkg
Pkg.instantiate()

using Base.Threads
import BenchmarkTools as BT
using CSV
using DataFrames
using HiGHS
using JuMP
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