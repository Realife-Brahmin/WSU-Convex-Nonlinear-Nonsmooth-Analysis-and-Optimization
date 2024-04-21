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

function polyhedronCentroid(vertices)
    # Extract x and y coordinates
    x_coords = [vertex[1] for vertex in vertices]
    y_coords = [vertex[2] for vertex in vertices]

    # Compute averages
    centroid_x = mean(x_coords)
    centroid_y = mean(y_coords)

    return (centroid_x, centroid_y)
end

function convertPolygonToConstraints(vertices)
    # Number of vertices
    p = length(vertices)

    # Check if the polygon is a line segment (only two vertices)
    if p == 2
        # Set up the equality constraint (line equation)
        Ae = zeros(1, 2)
        be = zeros(1)

        # Calculate the line equation coefficients using the given vertices
        x1, y1 = vertices[1]
        x2, y2 = vertices[2]
        Ae[1, :] = [y2 - y1, x1 - x2]
        be[1] = x1 * y2 - x2 * y1

        # No inequality constraints for a line segment
        A = zeros(0, 2)
        b = zeros(0)
    else
        # Inequality constraints for a polygon
        A = zeros(p, 2)
        b = zeros(p)

        for k in 1:p
            # Ensure vertices are cyclic by connecting the last to the first
            xk, yk = vertices[k]
            xkp1, ykp1 = vertices[mod(k, p)+1]

            # Calculate the coefficients for the inequalities
            A[k, :] = -[ykp1 - yk, xk - xkp1]
            b[k] = -(xk * ykp1 - yk * xkp1)
        end

        # No equality constraints for general polygons
        Ae = zeros(0, 2)
        be = zeros(0)
    end

    return Ae, be, A, b
end

function combineConstraints(A0e, b0e, Ape, bpe, A0i, b0i, Api, bpi)
    # Get the sizes of the problem
    nx = size(A0e, 2)  # Number of x variables
    ny = size(Ape, 2)  # Number of y variables
    mEx, mEy = length(b0e), length(bpe)
    mIx, mIy = length(b0i), length(bpi)

    Ae = vcat(hcat(A0e, zeros(mEx, ny)), hcat(zeros(mEy, nx), Ape))
    be = vcat(b0e, bpe)
    A = vcat(hcat(A0i, zeros(mIx, ny)), hcat(zeros(mIy, nx), Api))
    b = vcat(b0i, bpi)

    return Ae, be, A, b

end

function preparePolyProblem(;poly::Int=1,
    beta::Float64=1.0,
    factor::Float64=1.5)

    P0 = [(0, 2), (-4, 0), (-3, -2), (0, -2), (1, -1)]
    xcentral = (0, 0) # should ideally check that xcentral actually lies inside P0
    A0e, b0e, A0i, b0i = convertPolygonToConstraints(P0)
    P1 = [(3, 1), (5, 3), (2, 5), (2, 3)]
    P2 = [(1, -4), (1, -6), (4, -6), (4, -4)]
    P3 = [(-2, -4), (-6, -3), (-3, -6)]
    P4 = [(-4, 3), (-2, 4), (-2, 6), (-6, 3)]
    P = [P0, P1, P2, P3, P4]
    numOuterAreas = length(P)-1

    if poly < 1
        error("Polyhedron number must be a natural number not exceeding $(numOuterAreas)")
    end

    alpha = factor * beta

    if alpha < beta
        error("α must be greater than or equal to β")
    end

    x01, x02 = xcentral
    α, β = alpha, beta

    G = [α+β 0 -α 0;
        0 α+β 0 -α;
        -α 0 α 0;
        0 -α 0 α]

    c = [-β * x01, -β * x02, 0, 0]

    c0 = 1 // 2 * β * norm(xcentral)^2

    x0 = collect(polyhedronCentroid(P0))
    Ape, bpe, Api, bpi = convertPolygonToConstraints(P[poly+1])
    y0 = collect(polyhedronCentroid(P[poly+1]))

    Ae, be, A, b = combineConstraints(A0e, b0e, Ape, bpe, A0i, b0i, Api, bpi)
    lb = myfill(c, -Inf)
    ub = myfill(c, Inf)
    mE, mI = length(be), length(b)

    # @show G
    pQP = Dict(:G=>G, :c=>c, :lb=>lb, :ub=>ub, :c0=>c0, :mE=>mE, :Ae=>Ae, :be=>be, :mI=>mI, :A=>A, :b=>b, :poly=>poly, :numOuterAreas=>numOuterAreas, :P=>P, :xcentral=>xcentral, :alpha=>alpha, :beta=>beta)
    # pQP = @packDict "{G, c, lb, ub, c0, mE, Ae, be, mI, A, b, poly, numOuterAreas, P, xcentral, alpha, beta}" # does not work now that the variables are defined inside the function (local scope) instead of outside in Main (global scope)

    w0 = computeFeasiblePointForLinearConstraints(pQP)

    objective = QPObjectiveFunction
    objectiveOriginal = transmissionLines
    objectiveString = string(objectiveOriginal)

    pr = generate_pr(objective, w0, params=pQP, problemType="QP"; objectiveString=objectiveString)

    return pr

end

pr = preparePolyProblem()