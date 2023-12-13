using CSV
using DataFrames

include("objective.jl")

"""
    signalDenoise(x::Vector{Float64}, p; verbose::Bool=false, log::Bool=true, getGradientToo::Bool=true)

Denoises a signal vector `x` based on provided parameters and returns the denoised signal. Optionally, it also returns the gradient of the denoising function.

# Arguments
- `x::Vector{Float64}`: The input signal vector to be denoised.
- `p`: A dictionary containing parameters and data for the denoising process. 
- `verbose::Bool` (optional): If `true`, enables verbose logging. Defaults to `false`.
- `log::Bool` (optional): If `true`, logs certain information. Defaults to `true`.
- `getGradientToo::Bool` (optional): If `true`, the function also returns the gradient along with the denoised signal. Defaults to `true`.

# Returns
- If `getGradientToo` is `true`, returns a tuple `(f, g)` where `f` is the denoised signal and `g` is the gradient.
- If `getGradientToo` is `false`, returns only the denoised signal `f`.

# Examples
```julia
x = [1.0, 2.0, 3.0, 4.0]
p = Dict(:params => [1.5, 0.1, 0.1, 1.0, 2.0], :data => [1.0, 2.0, 3.0, 4.0])
denoised_signal = signalDenoise(x, p)
```
"""
function signalDenoise(x::Vector{Float64}, 
    p;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)

    n = length(x)
    # params = p.params
    params = p[:params]

    # d = p.data
    d = p[:data]
    
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

data0 = df.y
downsampling_rate = 10
data = data0[1:downsampling_rate:length(data0)] 
x0 = Float64.(data)

# p = 0.5
# p = 1
p = 2

# alpha = 0
# alpha = 0.5
# alpha = 1
# alpha = 2.0
# alpha = 5.0
alpha = 10.0
# alpha = 100

beta = minimum(abs.(x0))*1e-5

L = Float64[]
R = Float64[]

params = vcat([p, alpha, beta], L, R)
objective = signalDenoise;

pr = generate_pr(objective, x0, data=data, params=params)

# signalDenoise(pr.x0, pr.p)