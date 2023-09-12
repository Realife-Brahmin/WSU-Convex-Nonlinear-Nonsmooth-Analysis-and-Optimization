# main.jl
using BenchmarkTools
using CSV
using DataFrames
using LinearAlgebra
using Revise
using Symbolics

rawDataFolder = "rawData/";
filename = rawDataFolder*"FFD.csv";
df = CSV.File(filename) |> DataFrame;
rename!(df, [:t, :V])

function objFun(df;getGradientToo=true)
        @variables A₀, A, τ, ω, α, ϕ
        x = [A₀, A, τ, ω, α, ϕ] 
        fsym = A₀ + A*exp(-t/τ)sin((ω+α*t)t + ϕ)
        f = (df.V - substitute.(fsym, t => df.t)).^2
        ∇f = Symbolics.gradient.(f, x)
        return f, ∇f
end

f, ∇f = objFun(df);

fnum = build_function(f, x, expression=Val{false})
∇fnum = build_function(∇f, x, expression=Val{false})

