using Plots

function plotMinPathTimeTrajectory(res;
    plotTrajectory=false,
    savePlot=savePlot)

    pr = res.pr
    p = pr[:p]
    params = p[:params]
    v = params[:v]
    A = params[:A]
    B = params[:B]
    Plots.heatmap(v, yflip = true)
    xlabel!("X-axis")
    ylabel!("Y-axis")
    title!("Speed Heatmap")
end