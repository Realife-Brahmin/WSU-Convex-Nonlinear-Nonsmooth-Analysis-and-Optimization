include("objective.jl")

function TestFunction1(x::Vector{Float64}, 
    pDict; 
    getGradientToo::Bool=true)

    # p = p.params
    # p = p[:params]
    p = pDict[:p]
    a = 4 - 2.1 * x[1]^2 + (1/3) * x[1]^4
    c = 4 * x[2]^2 - 4
    f = a * x[1]^2 + x[1] * x[2] + c * x[2]^2 + p[1]

    if getGradientToo
        g = zeros(Float64, 2)
        g[1] = 2*x[1]*a + x[1]^2 * ((4/3) * x[1]^3 - 4.2*x[1]) + x[2]
        g[2] = x[1] + 2*x[2]*c + 8*x[2]^3
        return f, g
    else
        return f
    end
end

objective = TestFunction1
x0 = [0.5, 0.5]
# params = Float64.([2])
params = Dict(:p => 2.0)
pr = generate_pr(objective, x0, params=params)
