include("objective.jl")

function TestFunction2(x::Vector{Float64}, p; getGradientToo::Bool=true)
    
    p = p.params
    a, b, c = p
    den = 2 * length(x)
    f = sum(a* x.^4 + b* x.^2 + c*x) / den + 40

    if getGradientToo
        g = (4a* x.^3 + 2b* x + c*ones(length(x))) / den
        return f, g
    else
        return f
    end
end

objective = TestFunction2
x0 = -2 .+ 2 .* rand(15)
params = Float64.([1, -16, 5])

pr = generate_pr(objective, x0, params=params)