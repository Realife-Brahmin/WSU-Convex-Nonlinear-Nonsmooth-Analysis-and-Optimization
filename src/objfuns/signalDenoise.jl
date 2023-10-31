using CSV
using DataFrames

include("objective.jl")

FuncParam = NamedTuple{(:params, :data), Tuple{Vector{Float64}, Matrix{Float64}}}

function signalDenoise(x::Vector{Float64}, 
    p;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)

    n = length(x)
    params = p.params
    d = p.data
    if ~isempty(params) && length(params) == 5
        p, α, β, L, R = params
    elseif length(params) == 3
        L = d[1]
        R = d[end]
        p, α, β = params
    else
        @error "Couldn't define parameters correctly."
    end

    xfull = vcat(L, d, R)
    

    f = (1/2)*sum( (x-d).^2 )
    xdiff = diff(xfull)
    if mod(p, 2) == 0
        f += (α/p)*sum(xdiff.^p)
    else
        f += (α/p)*sum( (xdiff.^2 + myfill(xdiff, β).^2).^(p/2) )
    end
    
    if getGradientToo
        g = x - d
        if mod(p, 2) == 0
            for i = 1:n
                g[i] += α * (-xdiff[i+1]^(p-1) + xdiff[i]^p)
            end 
        else
            for i = 1:n
                g[i] += α *( xdiff[i+1]*(-xdiff[i+1]^2 + β^2)^(p/2-1) + xdiff[i]*(xdiff[i]^2 + β^2)^(p/2-1))
            end
        end
        return f, g
    else
        return f
    end

    @error "forbidden loc"

end

rawDataFolder = "rawData/"
filename = rawDataFolder * "FFD.csv"
df = CSV.File(filename) |> DataFrame
rename!(df, [:x, :y])


# data = vec(df.y) # dampedSHM data
data = [1.00, 1.01, 1.05, 0.93, 0.96, 1.10]; 
# data = [1.00, 1, 1, 1, 1]
data = [1, 2, 3]
x0 = Float64.(data)

# p = 0.5
# p = 1
p = 2

# α = 0
# α = 0.5
α = 1

β = 1e-5

L = Float64[]
R = Float64[]

params = vcat([p, α, β], L, R)
objective = signalDenoise;

pr = generate_pr(objective, x0, data=data, params=params)

# signalDenoise(pr.x0, pr.p)