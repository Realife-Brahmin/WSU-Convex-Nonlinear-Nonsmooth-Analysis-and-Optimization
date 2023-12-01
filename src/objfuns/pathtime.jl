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
    xxm = (xxm[2:end]+xxm[1:end-1])/2
    yym = (yym[2:end]+yym[1:end-1])/2
    xxm = max(min(xxm, mx), 1)
    yym = max(min(yym, my), 1)

    dist = sqrt(diff(xx).^2 + diff(yy).^2)
    vel = interp2(v, xxm, yym)
    f = sum(dist./vel)

    if getGradientToo
        del = sqrt(eps)
        g = zeros(2*n)
        for j = 1:2*n
            y = x
            y[j] += del
            df = pathtime(y, p)
            g[j] = (df-f)/del
        end
        return f, g
    else
        return f
    end

    @error "forbidden loc"

end

rawDataFolder = "rawData/"
filename = rawDataFolder * "SpeedData.csv"
df = CSV.File(filename) |> DataFrame

v = Matrix(df)
# rename!(df, [:x, :y])
n = 4
x0 = Float64.(0.1*randn(2*n))
A = (0.05, 0.05)
B = (0.95, 0.95)
params = Dict()
params = Dict(:v => v, :A=>A, :B=>B)

objective = pathtime;

pr = generate_pr(objective, x0, params=params)

# signalDenoise(pr.x0, pr.p)