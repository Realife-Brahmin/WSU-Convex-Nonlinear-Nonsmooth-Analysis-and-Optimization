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

    Plots.heatmap(v, yflip = true, color = palette(:vik25, rev=false), clim = extrema(v))

    # Add scatter points
    scatter!([Am[1]], [Am[2]], color=:white, label="Point A")
    scatter!([Bm[1]], [Bm[2]], color=:white, label="Point B")

    xlabel!("X-axis")
    ylabel!("Y-axis")
    title!("Speed Heatmap")
end

plotresults(res, savePlot=true)