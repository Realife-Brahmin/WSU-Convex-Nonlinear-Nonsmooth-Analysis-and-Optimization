using Plots
using Plots.PlotMeasures


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