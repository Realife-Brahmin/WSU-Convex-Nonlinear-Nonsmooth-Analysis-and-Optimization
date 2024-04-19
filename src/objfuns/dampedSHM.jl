include("objective.jl")
include("boxConstraintPenalty.jl")

function dampedSHM(x::Vector{Float64}, 
    p;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)

    # data = p[:params][:data]
    data = p[:data]
    M = size(data, 1)
    nf = size(data, 2) - 1
    y = data[:, nf+1]
    xf = data[:, 1:nf]
    t = xf

    n = length(x) # note that df's x (time) is different from function parameter x
    # mags = p[:params][:mags]
    mags = p[:mags]
    # A₀, A, τ, ω, α, ϕ = x.*mags
    x = x .* mags
    A₀, A, τ, ω, α, ϕ = x
    
    # constraints = p[:params][:constraints]
    constraints = p[:constraints]
    box = constraints[:box]

    f = 0.0;
    g = zeros(Float64, n)
    for k = 1:M
        tₖ = t[k]
        yₖ = y[k]
        expₖ = exp(-tₖ/τ)
        Sₖ = expₖ*sin((ω+α*tₖ)tₖ + ϕ)
        Cₖ = expₖ*cos((ω+α*tₖ)tₖ + ϕ)

        ŷₖ = A₀ + A*Sₖ
        Δyₖ = ŷₖ - yₖ
        f += (1/M)*Δyₖ^2

        P = boxConstraintPenalty(x, box)

        f += P

        if getGradientToo
            g += (2/M)* Δyₖ * [1, Sₖ, A*tₖ*(τ^-2)*Sₖ, A*tₖ*Cₖ, A*tₖ^2*Cₖ, A*Cₖ]

            gP = boxConstraintGradient(x, box)

            g .+= gP
        end
    end
    
    if getGradientToo
        return f, g
    else
        return f
    end
end

rawDataFolder = "rawData/"
filename = rawDataFolder * "FFD.csv"
df = CSV.File(filename) |> DataFrame
rename!(df, [:x, :y])
data = Matrix(df)
x00 = [13.8, 8.3, 0.022, 1800, 900, 4.2]

indices = [] # no constraint
# indices = [1, 3, 5]
# indices = collect(1:n)
# lazily automatically generating lower and upper bound values for the constrained decision variables from their own x0 value.
# lbs = x00[indices]*0.80
lbs = x00[indices]*0.95
# lbs = x00[indices] * 0.99
# ubs = x00[indices]*1.20
ubs = x00[indices] * 1.05
# ubs = x00[indices] * 1.01

barrier = 1e8
box = Dict(:indices=>indices, :lbs=>lbs, :ubs=>ubs, :barrier=>barrier)
constraints = Dict(:box => box)

n = length(x00)
mags = ones(n)
x0 = x00./mags;
params = Dict(:data=>data, :x00=>x00, :mags=>mags, :constraints=>constraints);

objective = dampedSHM;

pr = generate_pr(objective, x0, params=params)