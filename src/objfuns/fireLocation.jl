include("objective.jl")

function fireLocation(x::Vector{Float64}, 
    pDict;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)
    
    n = length(x)
    @unpack B, d = pDict
    r = Bx - d
    f = 1//2 *sum(r.^2)

    if !getGradientToo
        return f
    elseif getGradientToo
        g = sum(transpose(B)*r)
    else
        @error "floc"
    end
    
end

n = length(x0)
params = Dict();

objective = fireLocation;

pr = generate_pr(objective, x0, params=params)