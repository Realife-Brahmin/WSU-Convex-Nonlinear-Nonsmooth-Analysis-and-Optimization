include("objective.jl")

FuncParam = NamedTuple{(:params, :data), Tuple{Vector{Float64}, Matrix{Float64}}}

function sphere(x::Vector{Float64}, p::FuncParam;
    getGradientToo::Bool=true,
    verbose::Bool=false)

    f = dot(x, x)
    if getGradientToo
        g = 2x
        return f, g
    else
        return f
    end
end

objective = sphere;
x0 = collect(0.1:0.1:1)

pr = generate_pr(objective, x0)