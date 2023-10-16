FuncParam = NamedTuple{(:params, :data), Tuple{Vector{Float64}, Matrix{Float64}}}


function drag(x::Vector{Float64}, 
    p::FuncParam;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)

    n = length(x)
    Δz = 1.0/(n+1)
    D0 = subDrag(x, p, 0)
    Dn = subDrag(x, p, n)
    D = zeros(Float64, n-1)

    for k = 1:n-1
        D[k] = subDrag(x, p, k)
    end

    f = D0 + sum(D) + Dn

    return f
end

function subDrag(x::Vector{Float64}, 
    p::FuncParam, k::Int64;
    verbose::Bool=false,
    log::Bool=true)

    n = length(x)
    if k > n || k < 0
        @error "Can't obtain D$(k) for a $(n)-vector x"
    end
    Δz = 1.0/(n+1)
    xk = (k == 0) ? 0.0 : x[k]
    # xkm1 = (k == 1) ? 0.0 : x[k-1]
    xkp1 = (k == n) ? 1.0 : x[k+1]
    
    Sk = (xkp1 - xk)/Δz
    Dk = (Sk^2/(1+Sk^2))*(xkp1^2 - xk^2)/2.0

    return Dk
end

n = 3;
x = collect(LinRange(0.0, 1.0, n+2)[2:n+1]);
p = empty_FuncParam();
D = drag(x, p)

