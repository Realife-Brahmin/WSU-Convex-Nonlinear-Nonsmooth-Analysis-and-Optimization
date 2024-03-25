include("objective.jl")

using Parameters

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

angles_deg = [95, 84, 75, 63]
m = tand.(90 .- angles_deg)
coords = [(10, 12), (3, 9), (2, 3), (8, 1)]
coords_x, coords_y = first.(coords), last.(coords)
c = coords_y - m.*coords_x
B = hcat(-m, -ones(4))
d = c
x0 = [0, 0]

n = length(x0)
x = x0
x0 = B\d
params = Dict();

objective = fireLocation;

pr = generate_pr(objective, x0, params=params)