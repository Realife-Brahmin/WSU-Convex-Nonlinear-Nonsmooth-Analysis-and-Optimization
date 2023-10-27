FuncParam = NamedTuple{(:params, :data), Tuple{Vector{Float64}, Matrix{Float64}}}


function drag(x::Vector{Float64}, 
    p::FuncParam;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)

    n = length(x)
    δ = 1e-8
    Δz = 1.0/(n+1)
    D0 = subDrag(x, p, 0)
    Dn = subDrag(x, p, n)
    D = zeros(Float64, n-1)

    for k = 1:n-1
        D[k] = subDrag(x, p, k)
    end

    f = D0 + sum(D) + Dn

    if getGradientToo
        g = zeros(Float64, n)
        for k = 1:n
            xk_pert = copy(x)
            xk_pert[k] += δ
    
            Dkm1pert = subDrag(xk_pert, p, k-1)
            Dkm1 = (k == 1) ? D0 : D[k-1]

            Dkpert = subDrag(xk_pert, p, k)
            Dk = (k == n) ? Dn : D[k]
            
            g[k] = ((Dkm1pert - Dkm1) + (Dkpert - Dk)) / δ
        end
        return f, g
    else
        return f
    end

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
    xkp1 = (k == n) ? 1.0 : x[k+1]
    
    Sk = (xkp1 - xk)/Δz

    violation = max(0, xk-xkp1)
    Dk = (Sk^2/(1+Sk^2))*(xkp1^2 - xk^2)/2.0 + exp(violation) - 1

    return Dk
end

## NamedTuple pr (Problem) generation
objective = drag;
data = Matrix{Float64}(undef, 0, 0)
# n = 800
n = 300
# n = 2^10
# n = 2048
x0 = Float64.(collect(LinRange(0.0, 1.0, n+2)[2:n+1]))
params = Float64[]
p = (params=params, data=data)
println("Problem (pr::NamedTuple) generated")
pr = (p=p, x0=x0, objective=objective, alg=alg)