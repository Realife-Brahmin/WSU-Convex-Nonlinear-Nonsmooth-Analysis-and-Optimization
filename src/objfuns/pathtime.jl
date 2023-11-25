using CSV
using DataFrames

include("objective.jl")

function pathtime(x::Vector{Float64}, 
    p;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)

    n = length(x)/2
    params = p[:params]

    w = x[1:n]
    z = x[n+1:2n]

    v = params[:v]
    my, mx = size(v)
    A = params[:A] # A is a tuple (x_a, y_a)
    B = params[:B] # B is a tuple (x_b, y_b)

    s = LinRange(0, 1, 1000)
    xx = (1 .- s)*A(1) + s*B(1)
    yy = (1 .- s)*A(2) + s*B(2)

    k = 1:n
    S = sin.(π*k*s)
    for k = 1:n
        S = sin(k*π*s)
        xx += w[k]*S
        yy += z[k]*S
    end
    
    xxm = 1 .+ xx*(mx-1)
    yym = 1 .+ yy*(my-1)

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
        g = zeros(n,)
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

x0 = Float64.(data)

params = Dict()
params = Dict(:v => v, :A=>A, :B=>B)

objective = pathtime;

pr = generate_pr(objective, x0, data=data, params=params)

# signalDenoise(pr.x0, pr.p)