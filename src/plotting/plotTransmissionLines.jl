using Plots
using Plots.PlotMeasures

function plot_convex_hull!(points, label)
    x_coords = [p[1] for p in points]
    y_coords = [p[2] for p in points]
    # Close the polygon by repeating the first point at the end
    plot!([x_coords; x_coords[1]], [y_coords; y_coords[1]], label=label)
    scatter!(x_coords, y_coords, label="")
end

function plotTransmissionLines(resVec;
    savePlot=savePlot)


    res1 = resVec[1]
    @unpack pr = res1
    @unpack P = pr[:p]
    
    numPoly = length(resVec)
    # P0, P1, P2, P3, P4 = 
    theme(:dao)
    p1 = plot(title="All Convex Hulls", aspect_ratio=:equal)

    # Plot each set of points
    for polyNum ∈ 0:numPoly
        polyString = "P"*string(polyNum)
        plot_convex_hull!(P[polyNum+1], polyString)
    end

    # plot_convex_hull!(P0, "P0")
    # plot_convex_hull!(P1, "P1")
    # plot_convex_hull!(P2, "P2")
    # plot_convex_hull!(P3, "P3")
    # plot_convex_hull!(P4, "P4")

    # Display the plot
    display(plot!())

    folderName = string(dirname(dirname(@__DIR__)))*"/processedData/"
    filename = folderName*"transmissionLines_"*pr.alg[:method]*".png"

    if savePlot
        
        if isfile(filename)
            rm(filename)
        end
    
        Plots.savefig(p1, filename)
    end

    return p1

end

