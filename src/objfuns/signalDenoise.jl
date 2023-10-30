FuncParam = NamedTuple{(:params, :data), Tuple{Vector{Float64}, Matrix{Float64}}}

function signalDenoise(x::Vector{Float64}, 
    p::FuncParam;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)

    n = length(x)
    params = p.params
    dk = p.data
    if ~isempty(params) && length(params) == 5
        p, α, β, L, R = params
    elseif length(params) == 3
        L = dk[1]
        R = dk[end]
        p, α, β = params
    else
        @error "Couldn't define parameters correctly."
    end

    xfull = vcat(L, dk, R)
    

    f = (1/2)*sum( (x-dk).^2 )
    if p == 2
        f += (α/p)*diff

    if getGradientToo
        g = zeros(Float64, n)
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
data = vec(df.y)
x0 = data

objective = signalDenoise;

pr = generate_pr(objective, x0, data=data)