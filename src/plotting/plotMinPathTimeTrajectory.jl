using Plots

function plotMinPathTimeTrajectory(res;
    plotTrajectory=false,
    savePlot=savePlot)

    pr = res.pr
    p = pr[:p]
    params = p[:params]
    v = params[:v]
    mx, my = size(v)
    A = params[:A]
    B = params[:B]
    Am = 1 .+ (mx-1, my-1).*A
    Bm = 1 .+ (mx-1, my-1).*B

    h1 = Plots.heatmap(v, yflip = true, color = palette(:vik25, rev=false), clim = extrema(v))

    # Add scatter points
    h11 = Plots.scatter(h1, [Am[1]], [Am[2]], color=:white, label="Point A")
    h12 = Plots.scatter(h11, [Bm[1]], [Bm[2]], color=:white, label="Point B")

    xlabel!("X-axis")
    ylabel!("Y-axis")
    title!("Speed Heatmap")

    xxm_int, yym_int = computeOptimalTrajectoryIndices(res)
    h2 = Plots.plot(h12, xxm_int, yym_int, linecolor=:white, linewidth=2, label="Optimal Path", yflip = true)

end
