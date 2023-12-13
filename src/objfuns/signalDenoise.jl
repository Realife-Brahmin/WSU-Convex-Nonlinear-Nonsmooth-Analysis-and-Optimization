using CSV
using DataFrames

include("objective.jl")

function signalDenoise(x::Vector{Float64}, 
    p;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)

    n = length(x)
    params = p[:params]
    
    @unpack d, p, alpha, beta, L, R, magnitudeString = params

    if isempty(L) 
        L = d[1]
    end
    
    if isempty(R)
        R = d[end]
    end

    xfull = vcat(L, x, R)
    

    f = (1/2n)*sum( (x-d).^2 )
    xdiff = diff(xfull)
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

d0 = df.y
downsampling_rate = 10
d = d0[1:downsampling_rate:length(d0)] 
magnitudeString = ""

# normalize!(x0)
normalize!(d)
magnitudeString = "Normalized"

x0 = Float64.(d)


p = 0.5
# p = 1
# p = 2

# alpha = 0
# alpha = 0.05
alpha = 0.5
# alpha = 1
# alpha = 2.0
# alpha = 5.0
# alpha = 10.0
# alpha = 100.0
# alpha = 1000.0

beta = minimum(abs.(x0))*1e-5

L = Float64[]
R = Float64[]

params = Dict(:d => d, :p => p, :alpha => alpha, :beta => beta, :L => L, :R => R, :magnitudeString => magnitudeString)
objective = signalDenoise;

pr = generate_pr(objective, x0, params=params)