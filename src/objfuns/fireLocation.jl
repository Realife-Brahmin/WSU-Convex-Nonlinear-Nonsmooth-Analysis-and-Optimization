include("objective.jl")

using Parameters

function fireLocation(x::Vector{Float64}, 
    pDict;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)
    
    n = length(x)
    @unpack B, d = pDict
    r = B*x - d
    f = 1//2 *sum(r.^2)

    if !getGradientToo
        return f
    elseif getGradientToo
        g = transpose(B)*r
        return f, g
    else
        @error "floc"
    end
    
end

angles_deg = [95, 84, 75, 63]
m = tand.(90 .- angles_deg)
coords = [(10, 12), (3, 9), (2, 3), (8, 1)]
coords_x, coords_y = first.(coords), last.(coords)
c = coords_y - m.*coords_x
B = hcat(-m, -ones(4))
d = c

n = length(x0)
x = x0
x0 = B\d
params = Dict(:B=>B, :d=>d);

objective = fireLocation;

pr = generate_pr(objective, x0, params=params, problemType="ECQP")

#test run
f0 = fireLocation(x0, params)