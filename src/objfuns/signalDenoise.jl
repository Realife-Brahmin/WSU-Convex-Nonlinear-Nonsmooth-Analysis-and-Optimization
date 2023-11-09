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
    if length(params) == 5
        p, alpha, beta, L, R = params
    elseif length(params) == 3
        L = d[1]
        R = d[end]
        p, alpha, beta = params
    else
        @error "Couldn't define parameters correctly."
    end

    xfull = vcat(L, x, R)
    

    f = (1/2n)*sum( (x-d).^2 )
    @show xdiff = diff(xfull)
    if mod(p, 2) == 0
        # println("Even integral p value.")
        f += (alpha/(p*n)) * sum(xdiff.^p)
    else
        println("NOT even integral p value.")
        f += (alpha/(p*n)) * sum( (xdiff.^2 + myfill(xdiff, beta).^2).^(p/2) )
    end
    
    if getGradientToo
        g = x - d
        if mod(p, 2) == 0
            # println("Even integral p value.")
            for i = 1:n
                g[i] += alpha/n * (
                    xdiff[i]^(p-1)
                    -1*(xdiff[i+1])^(p-1) 
                )
            end 
        else
            println("NOT even integral p value.")
            for i = 1:n
                
                # @show xdiff[i], xdiff[i+1]
                g[i] += alpha/n *( 
                        xdiff[i]*(xdiff[i]^2 + beta^2)^(p/2-1)
                        -xdiff[i+1]*(xdiff[i+1]^2 + beta^2)^(p/2-1)
                )
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
# data = [1.00, 1.01, 1.05, 0.93, 0.96, 1.10]; 
# data = Float64.(abs.(rand(Int8, 15)))
# data = [1.00, 1, 1, 1, 1]
data = [1, 2, 3]
# data = [0.01, 0.02, 0.03]
x0 = Float64.(data)

# p = 0.5
# p = 1
p = 2

# alpha = 0
alpha = 0.5
# alpha = 1
# alpha = 100

beta = 1e-5

L = Float64[]
R = Float64[]

params = vcat([p, alpha, beta], L, R)
objective = signalDenoise;

pr = generate_pr(objective, x0, data=data, params=params)

# signalDenoise(pr.x0, pr.p)