using Plots
using Plots.PlotMeasures
using Colors
using ColorSchemes

"""
    plot_convex_hull!(points, label; color=:black)

Add a convex hull polygon and its vertices to an existing plot.

# Arguments
- `points::Vector{Tuple{Real,Real}}`: A vector of 2D points (tuples) that define the vertices of the convex hull.
- `label::String`: The label for the convex hull in the plot legend.

# Keyword Arguments
- `color::Symbol`: The color of the convex hull edges and vertices (default: `:black`).

# Description
This function takes a series of points defining the convex hull, connects them in order to form a closed polygon, and plots this polygon on the current plot. Additionally, it scatters the vertices with a specified color. The first point is repeated at the end to close the polygon. This function is intended to be used with an existing plot and will modify it in place.

# Examples
```julia
# Assume existing plot p
p = plot(aspect_ratio=:equal)

# Define points for convex hull
points = [(0, 0), (1, 2), (2, 1), (1, 0)]

# Plot convex hull on plot p
plot_convex_hull!(points, "My Hull", color=:blue)

# Display the plot
display(p)
```
"""
function plot_convex_hull!(points, label; color=:black)
    x_coords = [p[1] for p in points]
    y_coords = [p[2] for p in points]
    # Close the polygon by repeating the first point at the end
    plot!([x_coords; x_coords[1]], [y_coords; y_coords[1]], label=label, color=color, linewidth=2)
    scatter!(x_coords, y_coords, label="", color=color, markersize=2, alpha=0.5)
end

"""
    pickColor(arr, num::Int)

Select a color from a provided array based on an index number.

# Arguments
- `arr::Vector`: An array of colors or any elements from which to select.
- `num::Int`: The index number used to pick an element from the array.

# Returns
- An element from `arr` corresponding to the provided index `num`.

# Description
Given an array `arr` and an integer `num`, this function returns the element at the position in the array that corresponds to `num` modulo the length of `arr`. This is useful for cycling through a list of colors or elements when the number of items to color may exceed the length of the color array. It effectively makes the color array circular.

# Examples
```julia
color_array = [:red, :green, :blue]
picked_color = pickColor(color_array, 5)
# Since 5 modulo the length of color_array (3) is 2, picked_color will be :blue
```
"""
function pickColor(arr, num::Int)
    n = length(arr)
    numArr = mod(num, n)
    return arr[numArr+1]
end

"""
    calculatePlotLimits(P; buffer=1.0)

Calculate the x and y limits for a plot based on the provided set of convex hulls with an optional buffer.

# Arguments
- `P::Vector{Vector{Tuple{Real,Real}}}`: A vector of vectors, where each inner vector contains tuples representing points of a convex hull.

# Keyword Arguments
- `buffer::Float64`: An optional buffer space to add around the calculated limits (default: `1.0`).

# Returns
- `(xlims, ylims)`: A tuple containing the x and y limits for the plot.

# Description
This function iterates through a set of points defining multiple convex hulls, identifies the minimum and maximum x and y values, and calculates the limits for plotting these hulls. It adds an optional buffer to each limit to ensure that the hulls are not plotted on the very edge of the plot, enhancing visual appeal and readability.

# Examples
```julia
# Define a set of convex hulls
P = [
    [(0, 0), (1, 2), (2, 1)],
    [(1, 1), (2, 3), (3, 2)]
]

# Calculate plot limits with a buffer
xlims, ylims = calculatePlotLimits(P, buffer=0.5)

# Set the calculated limits in a plot
plot(xlims=xlims, ylims=ylims)
```
"""
function calculatePlotLimits(P; buffer=1.0)
    # Initialize min and max values
    x_min, x_max, y_min, y_max = Inf, -Inf, Inf, -Inf

    # Iterate through all hulls and update the min and max values
    for hull in P
        for (x, y) in hull
            x_min = min(x_min, x)
            x_max = max(x_max, x)
            y_min = min(y_min, y)
            y_max = max(y_max, y)
        end
    end

    # Apply the buffer to each limit
    xlims = (x_min - buffer, x_max + buffer)
    ylims = (y_min - buffer, y_max + buffer)

    return xlims, ylims
end

"""
    annotateWithOffset(x, y, txt; xoff=0.5, yoff=0.0)

Annotate a point on a plot with an offset to prevent overlapping with the point's marker.

# Arguments
- `x`: The x-coordinate of the point to annotate.
- `y`: The y-coordinate of the point to annotate.
- `txt`: The text of the annotation.

# Keyword Arguments
- `xoff::Float64`: The horizontal offset from the point (default: `0.5`).
- `yoff::Float64`: The vertical offset from the point (default: `0.0`).

# Description
`annotateWithOffset` is a utility function that adds a text annotation to a plot with specified offsets from a point. This is useful to avoid overlapping text with the point's marker or other plot elements. The offsets can be adjusted to place the text at an appropriate distance from the point.

# Examples
```julia
# Annotate a point at coordinates (3, 2) with the text "A1"
annotateWithOffset(3, 2, "A1")

# Annotate the same point with a larger horizontal offset and a vertical offset
annotateWithOffset(3, 2, "A1", xoff=1.0, yoff=0.5)
```
"""
function annotateWithOffset(x, y, txt; xoff=0.5, yoff=0.0)
    annotate!(x + xoff, y + yoff, text(txt, 12, :center, :bottom))
end

"""
    calculateCentroid(vertices)

Calculate the centroid (geometric center) of a set of vertices representing a polygon.

# Arguments
- `vertices::Vector{Tuple{Real, Real}}`: A vector of tuples where each tuple represents the x and y coordinates of a vertex.

# Returns
- `(centroid_x, centroid_y)`: A tuple representing the coordinates of the centroid.

# Description
Given a list of vertices that define a polygon, `calculateCentroid` computes the average of the x and y coordinates to find the centroid of the polygon. This function assumes that the input represents a simple, non-self-intersecting polygon.

# Examples
```julia
# Define vertices of a triangle
vertices = [(0, 0), (2, 0), (1, 2)]

# Calculate the centroid of the triangle
center = calculateCentroid(vertices)
println("Centroid of the polygon is at: ", center)

# Output
Centroid of the polygon is at: (1.0, 0.6666666666666666)
```
"""
function calculateCentroid(vertices)
    x_coords = [vertex[1] for vertex in vertices]
    y_coords = [vertex[2] for vertex in vertices]
    centroid_x = mean(x_coords)
    centroid_y = mean(y_coords)
    return (centroid_x, centroid_y)
end

"""
    plotTransmissionLines(resVec; savePlot=false)

Plot the optimal locations for transmission towers based on the results of an optimization process. The plot includes all convex hulls from a central point to the optimal points on surrounding polyhedrons.

# Arguments
- `resVec::Vector{ResultsType}`: A vector containing results from an optimization process. Each result should have a property `pr` with optimization parameters.

# Keyword Arguments
- `savePlot::Bool=false`: A flag to determine whether to save the plot to a file.

# Returns
- `Plot`: A plot object generated by the Plots.jl package.

# Description
The function takes a vector of result objects from an optimization process and plots the convex hulls of given points (P), lines connecting the central point (xcentral) to optimal points (x_i and y_i) on each polyhedron, and annotates the points and polyhedrons with their respective labels.

# Usage
Each element in `resVec` should follow a specific structure that includes a `pr` property. This property must contain at least `P`, `xcentral`, `alpha`, and `beta`. The function automatically calculates plot limits to fit all elements neatly within the plot area.

# Notes
- This function is a visualization tool to represent the output of an optimization process, specifically for the location of transmission towers.
- Colors are assigned from predefined palettes for distinction between different polyhedrons and their optimal points.
- The plot is displayed immediately and saved if `savePlot` is `true`.

# Examples
```julia
# Assume resVec is already populated with results from an optimization process
plotTransmissionLines(resVec)

# To save the plot to a file
plotTransmissionLines(resVec, savePlot=true)
```
"""
function plotTransmissionLines(resVec; savePlot=false)

    res1 = resVec[1]
    @unpack pr = res1
    @unpack P, xcentral, alpha, beta = pr[:p]  # Ensure correct extraction path if needed
    factor = alpha / beta

    xlims, ylims = calculatePlotLimits(P)

    numPoly = length(resVec)

    theme(:dao)

    p1 = plot(aspect_ratio=:equal)


    foptAll = 0.0

    color_palette = [
        colorant"Firebrick",
        colorant"CornflowerBlue",
        colorant"purple4",
        colorant"MidnightBlue",
        colorant"DarkOliveGreen",
        colorant"SaddleBrown",
    ]

    # Define a second color palette with different shades
    color_palette2 = [
        colorant"Firebrick",
        colorant"CornflowerBlue",
        colorant"purple4",
        colorant"MidnightBlue",
        colorant"DarkOliveGreen",
        colorant"SaddleBrown",
    ]


    polyPlotted = 0
    color = pickColor(color_palette, polyPlotted)
    color2 = pickColor(color_palette2, polyPlotted)

    # Plot central polyhedron and its central point
    plot_convex_hull!(P[1], "", color=color)
    scatter!([xcentral[1]], [xcentral[2]], label=L"x_0", markersize=7, color=color)
    annotateWithOffset(xcentral[1], xcentral[2], L"x_0", xoff=-0.5, yoff=-0.2)

    centroid = calculateCentroid(P[polyPlotted+1])
    polyLabel = LaTeXString("\$P_{$polyPlotted}\$")
    annotateWithOffset(centroid[1], centroid[2], polyLabel)

    polyPlotted += 1


    # Plot each set of points for other polyhedrons
    for polyNum in 1:numPoly
        # polyString = "P" * string(polyNum)
        polyString = ""
        xPolyStringLatex = LaTeXString("\$x_{$polyNum}\$")  # Using LaTeXStrings to create a dynamic label
        yPolySringLatex = LaTeXString("\$y_{$polyNum}\$")  # Using LaTeXStrings to create a dynamic label

        color = pickColor(color_palette, polyPlotted)
        color2 = pickColor(color_palette2, polyPlotted)

        plot_convex_hull!(P[polyNum+1], polyString, color=color)

        centroid = calculateCentroid(P[polyPlotted+1])
        polyLabel = LaTeXString("\$P_{$polyPlotted}\$")
        annotateWithOffset(centroid[1], centroid[2], polyLabel)

        resPoly = resVec[polyNum]
        woptPoly = resPoly.xvals[:, end]
        foptPoly = resPoly.fvals[end]

        # Extract and plot optimal points for this polyhedron
        xoptPoly, yoptPoly = woptPoly[1:2], woptPoly[3:4]

        scatter!([xoptPoly[1]], [xoptPoly[2]], label=xPolyStringLatex, color=color, markershape=:rect, markersize=5)
        
        scatter!([yoptPoly[1]], [yoptPoly[2]], label=yPolySringLatex, color=color2,
        markershape=:star5, markersize=9)

        # Annotate the points directly on the plot
        annotateWithOffset(xoptPoly[1], xoptPoly[2], xPolyStringLatex, xoff=-0.5, yoff=-0.1)

        annotateWithOffset(yoptPoly[1], yoptPoly[2], yPolySringLatex)


        plot!([xcentral[1], xoptPoly[1]], [xcentral[2], xoptPoly[2]], color=color, label="", linewidth=2, linestyle=:dash)

        plot!([xoptPoly[1], yoptPoly[1]], [xoptPoly[2], yoptPoly[2]], color=color2, label="", linewidth=2, linestyle=:dash)

        polyPlotted += 1
        foptAll += foptPoly
    end

    formatted_foptAll = @sprintf("%.4f", foptAll)

    plot!(p1, title="Optimal Transmission Tower Locations\n" * "α = $(alpha) and β = $(beta) | f = $(formatted_foptAll)",
    xlim=xlims,
    ylim=ylims,
    size=(600,600))
    # Display the plot
    display(p1)

    # Saving the plot if requested
    if savePlot
        folderName = string(dirname(dirname(@__DIR__))) * "/processedData/"
        filename = folderName * "transmissionLines_" * pr.alg[:method] * "_factor_" * string(factor) * ".png"

        if isfile(filename)
            rm(filename)  # Remove existing file to avoid conflict
        end

        savefig(p1, filename)
    end

    return p1
end


