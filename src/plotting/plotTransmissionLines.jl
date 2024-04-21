using Plots
using Plots.PlotMeasures

function plot_convex_hull!(points, label; color=:black)
    x_coords = [p[1] for p in points]
    y_coords = [p[2] for p in points]
    # Close the polygon by repeating the first point at the end
    plot!([x_coords; x_coords[1]], [y_coords; y_coords[1]], label=label, color=color)
    scatter!(x_coords, y_coords, label="", color=color)
end

function pickColor(arr::Vector, num::Int)
    n = length(arr)
    numArr = mod(num, n)
    return arr[numArr+1]
end

function plotTransmissionLines(resVec; savePlot=false)

    res1 = resVec[1]
    @unpack pr = res1
    @unpack P, xcentral, alpha, beta = pr[:p]  # Ensure correct extraction path if needed
    factor = alpha / beta
    numPoly = length(resVec)

    theme(:dao)

    # Define a color palette for the polyhedrons and their points   

    p1 = plot(aspect_ratio=:equal)


    foptAll = 0.0

    color_palette = [:red, :blue, :green, :orange, :purple, :yellow, :cyan, :magenta]
    polyPlotted = 0
    color = pickColor(color_palette, polyPlotted)

    # Plot central polyhedron and its central point
    # plot_convex_hull!(P[1], "P0", color=color)
    plot_convex_hull!(P[1], "", color=color)
    scatter!([xcentral[1]], [xcentral[2]], label=L"x_0", markersize=8, color=color)

    polyPlotted += 1

    # Plot each set of points for other polyhedrons
    for polyNum in 1:numPoly
        # polyString = "P" * string(polyNum)
        polyString = ""
        xPolyStringLatex = LaTeXString("\$x_{$polyNum}\$")  # Using LaTeXStrings to create a dynamic label
        yPolySringLatex = LaTeXString("\$y_{$polyNum}\$")  # Using LaTeXStrings to create a dynamic label

        color = pickColor(color_palette, polyPlotted)

        plot_convex_hull!(P[polyNum+1], polyString, color=color)
        resPoly = resVec[polyNum]
        woptPoly = resPoly.xvals[:, end]
        foptPoly = resPoly.fvals[end]

        # Extract and plot optimal points for this polyhedron
        xoptPoly, yoptPoly = woptPoly[1:2], woptPoly[3:4]

        scatter!([xoptPoly[1]], [xoptPoly[2]], label=xPolyStringLatex, color=color, markersize=6)
        scatter!([yoptPoly[1]], [yoptPoly[2]], label=yPolySringLatex, color=color, markersize=6)

        polyPlotted += 1
        foptAll += foptPoly
    end

    formatted_foptAll = @sprintf("%.4f", foptAll)

    plot!(p1, title="Optimal Transmission Tower Locations\n" * "α = $(alpha) and β = $(beta) | f = $(formatted_foptAll)")
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


