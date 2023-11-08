function plotDragCurve(res;
    savePlot::Bool=true)

    pr = res.pr

    theme(:dao)
    f = res.fvals[end]
    r = res.xvals[:, end]
    n = length(r)
    z = collect(0.0:1.0/(n+1):1.0)[2:end-1]

    gr()
    p1 = Plots.plot(z, r, dpi = 7000,
        title = "Drag Functional Estimate: $(L"f_{min} =") $(round(f, digits=5)) \n using $(pr.alg.method) with $(n) points.",
        titlefont = font(12,"Computer Modern"),
        guidefont = font(15,"Computer Modern"),
        label = L"r(z)",
        ylabel = L"r(z)",
        xlabel = L"z",
        line = :stem,
        linewidth = 0.3,
        marker = :circle,
        markersize = 2,
        alpha = 1.0,
        aspect_ratio = :equal,
        xlims = (-0.01, 1.01),
        ylims = (-0.00, 1.01)
    )

    folderName = string(dirname(dirname(@__DIR__)))*"/processedData/"
    filename = folderName*"dragFunction_"*pr.alg.method*"_"*string(n)*".pdf"


    if savePlot
        
        if isfile(filename)
            rm(filename)
        end
    
        Plots.savefig(p1, filename)
    end

    return p1
end