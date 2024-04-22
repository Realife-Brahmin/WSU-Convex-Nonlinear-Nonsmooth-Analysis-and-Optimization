using Plots
using Plots.PlotMeasures
using Colors
using ColorSchemes

function plot_convex_hull!(points, label; color=:black)
    x_coords = [p[1] for p in points]
    y_coords = [p[2] for p in points]
    # Close the polygon by repeating the first point at the end
    plot!([x_coords; x_coords[1]], [y_coords; y_coords[1]], label=label, color=color, linewidth=2)
    scatter!(x_coords, y_coords, label="", color=color, markersize=2, alpha=0.5)
end

function pickColor(arr, num::Int)
    n = length(arr)
    numArr = mod(num, n)
    return arr[numArr+1]
end

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

function annotateWithOffset(x, y, txt; xoff=0.5, yoff=0.0)
    annotate!(x + xoff, y + yoff, text(txt, 12, :center, :bottom))
end

function calculateCentroid(vertices)
    x_coords = [vertex[1] for vertex in vertices]
    y_coords = [vertex[2] for vertex in vertices]
    centroid_x = mean(x_coords)
    centroid_y = mean(y_coords)
    return (centroid_x, centroid_y)
end

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


