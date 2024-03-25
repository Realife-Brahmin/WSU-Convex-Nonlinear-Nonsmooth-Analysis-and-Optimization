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

coords = [(10,12), (3,9), (2,3), (8,1)]
angles_deg = [95, 84, 75, 63]
m = tand.(90 .- angles_deg)
coords_x = [x for (x, y) ∈ coords]
coords_y = [y for (x, y) ∈ coords]
n = length(x0)
params = Dict();

objective = fireLocation;

pr = generate_pr(objective, x0, params=params)