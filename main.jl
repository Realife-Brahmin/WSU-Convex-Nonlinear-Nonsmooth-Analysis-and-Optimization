# main.jl
using BenchmarkTools
using CSV
using DataFrames
using LinearAlgebra
using Revise
using Symbolics

include("src/initializer.jl")
include("src/objective.jl")

rawDataFolder = "rawData/";
filename = rawDataFolder*"FFD.csv";
df = CSV.File(filename) |> DataFrame;
rename!(df, [:t, :V]);

# function objFun(df::DataFrame;getGradientToo=true)
#     N = size(df, 1)
#     # N = size(df, 1)
#     @variables t, A₀, A, τ, ω, α, ϕ
#     x = [A₀, A, τ, ω, α, ϕ] 
#     fsym = A₀ + A*exp(-t/τ)sin((ω+α*t)t + ϕ)
#     f = sum([(df.V[i] - substitute(fsym, t => df.t[i]))^2 for i in 1:N])

#     if getGradientToo
#         # Calculate the gradient for each squared difference
#         gradients = [Symbolics.gradient(f[i], x) 
#         for i ∈ eachindex(f)]   
#         # Sum up the gradients to get the overall gradient
#         ∇f = sum(gradients)
#         return f, ∇f
#     else
#         return f
#     end
# end

f, ∇f = objFun(df);

fnum = build_function(f, x, expression=Val{false})
∇fnum = build_function(∇f, x, expression=Val{false})

A₀₀ = mean(df.V)
fnum([A₀₀, 1, 1, 1, 1, 1])