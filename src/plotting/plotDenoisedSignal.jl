using Plots
using Plots.PlotMeasures

function plotDenoisedSignal(res;
    savePlot::Bool=true)

    pr = res.pr
    method = pr.alg[:method]
    params = pr.p[:params]
    @unpack d, p, alpha, magnitudeString = params
    n = length(d)
    w = res.xvals[:, end]
    t = collect(1:1:n)
    minx, maxx = extrema(t)
    miny, maxy = extrema(d)

    theme(:dao)

    p1 = Plots.scatter(t, d,
        title = "Denoised Signal  $(magnitudeString) \n "*"using "*L"p = "*"$(p) and α = $(alpha) \n"* "Method: $(method)",
        label = "data-point",
        xlabel = "Time Instance [unit]",
        ylabel = "Signal value [unit]")

    if magnitudeString == ""
        size = (450, 400)
        aspect_ratio = :equal
    else
        size = (450, 400)
        aspect_ratio = :auto
    end

    p2 = Plots.plot!(p1, w,
        label = "fitted-function",
        linewidth = 3.0,
        linestyle = :solid,
        alpha = 0.5,
        legend = :topright,
        aspect_ratio = aspect_ratio,
        size = size,
        xlims=(minx*0.95, maxx*1.05), 
        ylims=(miny*0.95, maxy*1.05))

    folderName = string(dirname(dirname(@__DIR__)))*"/processedData/"
    filename = folderName*"denoisedSignal$(magnitudeString)_"*pr.alg[:method]*"_"*string(n)*"_p_$(replace(string(p), "." => "-"))_alpha_$(replace(string(alpha), "." => "-")).png"


    if savePlot
        
        if isfile(filename)
            rm(filename)
        end
    
        Plots.savefig(p2, filename)
    end
    
    return p2
end