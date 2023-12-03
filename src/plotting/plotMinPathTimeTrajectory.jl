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

function computeOptimalTrajectoryIndices(res)
    pr = res.pr
    wz_opt = res.xvals[:, end]
    
    n = Int(length(wz_opt)/2)

    p = pr[:p]
    params = p[:params]
    m = params[:m]
    w = wz_opt[1:n]
    z = wz_opt[n+1:2n]

    v = params[:v]
    my, mx = size(v)
    A = params[:A] # A is a tuple (x_a, y_a)
    B = params[:B] # B is a tuple (x_b, y_b)

    s = collect(LinRange(0, 1, m))
    pi_array = collect(1:n)*π
    xx = (1 .- s)*A[1] + s*B[1] + sin.(s * pi_array') * w
    yy = (1 .- s)*A[2] + s*B[2] + sin.(s * pi_array') * z
    xxm = 1 .+ xx*(mx-1)
    yym = 1 .+ yy*(my-1) 

    xxm_int = round.(Int, xxm)
    yym_int = round.(Int, yym)

    return xxm_int, yym_int
end
# plotresults(res, savePlot=true)