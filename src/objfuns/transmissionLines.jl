include("objective.jl")
include("activeSetQuadraticProgramming.jl")

using Parameters

function transmissionLines(x::Vector{Float64}, 
    pDict;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)
    
    n = length(x)
    @unpack B, d = pDict[:params]
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

P0 = [(0, 2), (-4, 0), (-3, -2), (0, -2), (1, -1)]
P1 = [(3, 1), (5, 3), (2, 5), (2, 3)]
P2 = [(1, -4), (1, -6), (4, -6), (4, -4)]
P3 = [(-2, -4), (-6, -3), (-3, -6)]
P4 = [(-4, 3), (-2, 4), (-2, 6), (-6, 3)]

using Plots

function plot_convex_hull!(points, label)
    x_coords = [p[1] for p in points]
    y_coords = [p[2] for p in points]
    # Close the polygon by repeating the first point at the end
    plot!([x_coords; x_coords[1]], [y_coords; y_coords[1]], label=label)
    scatter!(x_coords, y_coords, label="")
end

# Create an empty plot
theme(:dao)
plot(title="All Convex Hulls", aspect_ratio=:equal)

# Plot each set of points
plot_convex_hull!(P0, "P0")
plot_convex_hull!(P1, "P1")
plot_convex_hull!(P2, "P2")
plot_convex_hull!(P3, "P3")
plot_convex_hull!(P4, "P4")

# Display the plot
display(plot!())

function convertPolygonToInequalities(vertices)
    # Number of vertices
    p = length(vertices)

    # Initialize A and b
    A = zeros(p, 2)
    b = zeros(p)

    for k in 1:p
        # Compute indices for current and next vertex
        current_vertex = vertices[k]
        next_vertex = vertices[mod(k, p)+1] # Ensures that after the last vertex, it goes back to the first

        # Calculate the coefficients for the inequalities
        A[k, :] = [next_vertex[2] - current_vertex[2], current_vertex[1] - next_vertex[1]]
        b[k] = A[k, 1] * current_vertex[1] + A[k, 2] * current_vertex[2]
    end

    return A, b
end

# Example usage:
# Define the vertices of the convex polygon P_j
# P_j = [(1, 2), (4, 3), (3, 7), (1, 5)]

# Convert to inequalities
# A0, b0 = convertPolygonToInequalities(P0)
# A1, b1 = convertPolygonToInequalities(P1)
(A0, b0), (A1, b1), (A2, b2), (A3, b3), (A4, b4) = convertPolygonToInequalities.([P0, P1, P2, P3, P4])
