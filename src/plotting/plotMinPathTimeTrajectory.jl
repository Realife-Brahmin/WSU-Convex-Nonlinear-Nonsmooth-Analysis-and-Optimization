using Plots

function plotMinPathTimeTrajectory(res;
    plotTrajectory=false,
    savePlot=savePlot)
    Plots.theme(:dao)
    pr = res.pr
    p = pr[:p]
    params = p[:params]
    v = params[:v]
    mx, my = size(v)
    A = params[:A]
    B = params[:B]
    Am = 1 .+ (mx-1, my-1).*A
    Bm = 1 .+ (mx-1, my-1).*B

    fval = round(res.fvals[end], digits=6)

    h1 = Plots.heatmap(v, 
    yflip = true, 
    color = palette(:vik25, rev=false), 
    clim = extrema(v), 
    aspect_ratio=:equal, 
    size=(600,600), 
    xlims=(1, mx), 
    ylims=(1, my),
    )

    # Add scatter points
    h11 = Plots.scatter(h1, [Am[1]], [Am[2]], label="Point A", markercolor=:antiquewhite2, markerstrokecolor=:black, markerstrokewidth=0.5, markersize=4)
    h12 = Plots.scatter(h11, [Bm[1]], [Bm[2]], label="Point B", markercolor=:antiquewhite4, markerstrokecolor=:black, markerstrokewidth=0.5, markersize=4)

    # legend!(:outertopright)
    xlabel!("X-axis")
    ylabel!("Y-axis")
    title!("Minimum Pathtime Trajectory\n"*"Optimal Time = $(fval) units\n"*"Method: $(pr.alg.method)\n"*"n = $(n)")

    xxm_int, yym_int = computeOptimalTrajectoryIndices(res)
    h2 = Plots.plot(h12,
    yflip = true, 
    xxm_int, 
    yym_int, 
    linecolor=:antiquewhite2, 
    linewidth=2, 
    label="Optimal Path",
    legend=:topright)

    folderName = string(dirname(dirname(@__DIR__)))*"/processedData/"
    filename = folderName*"minTimePath_"*pr.alg.method*"_n_"*string(n)*".png"

    if savePlot
        
        if isfile(filename)
            rm(filename)
        end
    
        Plots.savefig(h2, filename)
    end

    return h2

end

plotresults(res, savePlot=true)
