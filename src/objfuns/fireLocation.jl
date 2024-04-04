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
m_line = tand.(90 .- angles_deg)
coords = [(10, 12), (3, 9), (2, 3), (8, 1)]
coords_x, coords_y = first.(coords), last.(coords)
c_line = coords_y - m_line.*coords_x
B = hcat(-m_line, -ones(4))
m, n = size(B)
d = c_line

x0 = B \ d
r0 = B * x0 - d
w0 = vcat(x0, r0)

params = Dict(:B=>B, :d=>d);

objective = fireLocation;

# pr = generate_pr(objective, x0, params=params, problemType="ECQP")

#test run
f0 = fireLocation(x0, params, getGradientToo=false)

# Henceforth, x, n, m lose their original meanings, in favour of w, n+m, n+m
G = vcat(hcat(zeros(n, n), zeros(n, m)), hcat(zeros(m, n), I(m)));
A = hcat(B, -I(m))
c = zeros(n+m)

pEQCP = Dict(:G=>G, :c=>c, :A=>A, :b=>b )

function equalityConstrainedQP(w::Vector{Float64},
    pDict;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=false)

    @unpack G, A, c, b = pDict
    f = 1//2*transpose(w)*G*w + transpose(w)*c

    if getGradientToo
        g = G*w + c
        return f, g
    else
        return f
    end
    
end

solState = SolStatePGCGType(w0, G, A, c, fk=f0)
