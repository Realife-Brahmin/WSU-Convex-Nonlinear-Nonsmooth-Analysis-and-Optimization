using CSV
using DataFrames

# include("objective.jl")
# include("dampedSHM.jl")

# rawDataFolder = "rawData/"
# filename = rawDataFolder * "FFD.csv"
# df = CSV.File(filename) |> DataFrame
# rename!(df, [:x, :y])
# data = Matrix(df)

function dampedSHM_signal(x::Vector{Float64}, k::Int,
    p)

    data = p.data
    M = size(data, 1)
    nf = size(data, 2) - 1
    y = data[:, nf+1]
    xf = data[:, 1:nf]
    t = xf

    n = length(x) # note that df's x (time) is different from function parameter x
    A₀, A, τ, ω, α, ϕ = x
    f = 0.0;
    g = zeros(Float64, n)

    tₖ = t[k]
    yₖ = y[k]
    expₖ = exp(-tₖ/τ)
    Sₖ = expₖ*sin((ω+α*tₖ)tₖ + ϕ)
    Cₖ = expₖ*cos((ω+α*tₖ)tₖ + ϕ)

    ŷₖ = A₀ + A*Sₖ
    f = ŷₖ
    
    return f
end

# x0 = [13.8, 8.3, 0.022, 1800, 900, 4.2]

# objective = dampedSHM;

# pr = generate_pr(objective, x0, data=data)

y0 = [dampedSHM_signal(pr.x0, k, pr.p) for k ∈ 1:200]
yopt = [dampedSHM_signal(res.xvals[:, end], k, pr.p) for k ∈ 1:200]
