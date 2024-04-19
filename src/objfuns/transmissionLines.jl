include("objective.jl")
include("functionQP.jl")
include("../activeSetQP.jl")

using Parameters

function transmissionLines(w::Vector{Float64}, 
    pDict;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)
    
    n = length(x)
    # @unpack B, d = pDict[:params]
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

function convertPolygonToConstraints(vertices)
    # Number of vertices
    p = length(vertices)

    # Initialize A and b
    A = zeros(p, 2)
    b = zeros(p)
    Ae = zeros(0, 2)
    be = zeros(0)

    # Check if the polygon is a line segment
    if p == 2
        # Polygon is a line segment; set Ae and be instead of A and b
        Ae = zeros(1, 2)
        be = zeros(1)

        # Compute the line equation coefficients
        Ae[1, :] = [vertices[2][2] - vertices[1][2], vertices[1][1] - vertices[2][1]]
        be[1] = Ae[1, 1] * vertices[1][1] + Ae[1, 2] * vertices[1][2]

        # Clear the A and b for inequalities since it's a line segment
        A = zeros(0, 2)
        b = zeros(0)
    else
        # Compute inequalities for a polygon
        for k in 1:p
            # Compute indices for current and next vertex
            current_vertex = vertices[k]
            next_vertex = vertices[mod(k, p)+1] # Ensures that after the last vertex, it goes back to the first

            # Calculate the coefficients for the inequalities
            A[k, :] = [next_vertex[2] - current_vertex[2], current_vertex[1] - next_vertex[1]]
            b[k] = A[k, 1] * current_vertex[1] + A[k, 2] * current_vertex[2]
        end
    end

    return Ae, be, A, b
end

P0 = [(0, 2), (-4, 0), (-3, -2), (0, -2), (1, -1)]
xcentral = (0, 0) # should ideally check that xcentral actually lies inside P0
A0e, b0e, A0i, b0i = convertPolygonToConstraints(P0)

beta = 1.0
factor = 1.5
alpha = factor * beta
x01, x02 = xcentral
α, β = alpha, beta
G = [α+β 0 -α 0;
    0 α+β 0 -α;
    -α 0 α 0;
    0 -α 0 α]

c = [-β * x01, -β * x02, 0, 0]

c0 = 1 // 2 * β * norm(xcentral)^2
P1 = [(3, 1), (5, 3), (2, 5), (2, 3)]
P2 = [(1, -4), (1, -6), (4, -6), (4, -4)]
P3 = [(-2, -4), (-6, -3), (-3, -6)]
P4 = [(-4, 3), (-2, 4), (-2, 6), (-6, 3)]
P = [P0, P1, P2, P3, P4]

poly = 1 # Polyhedron Number, 0 reserved for the community with the power station

Ape, bpe, Api, bpi = convertPolygonToConstraints(P[poly])
Ae, be, A, b = Ape, bpe, Api, bpi
lb = myfill(c, -Inf)
ub = myfill(c, Inf)
mE, mI = length(be), length(b)

pQP = @packDict "{G, c, lb, ub, c0, mE, Ae, be, mI, A, b, poly}"

x0, f0 = computeFeasiblePointForLinearConstraints(pQP)

objective = QPObjectiveFunction
objectiveOriginal = transmissionLines
objectiveString = string(objectiveOriginal)
# params = pQP

pr = generate_pr(objective, x0, params=pQP, problemType="QP"; objectiveString=objectiveString)