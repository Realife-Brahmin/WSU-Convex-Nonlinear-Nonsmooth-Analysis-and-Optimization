using Plots

function plotDenoisedSignal(res;
    savePlot::Bool=true)

    pr = res.pr
    d = pr.x0
    params = pr.p.params
    p, alpha = params[1:2]
    n = length(d)
    w = res.xvals[:, end]
    t = collect(1:1:n)
    
    theme(:dao)

    p1 = Plots.scatter(t, d,
        title = "Denoised Signal \n"*"using "*L"p = "*"$(p) and α = $(alpha)",
        label = "data-point",
        xlabel = "Time Instance [unit]",
        ylabel = "Signal value [unit]")

    p2 = Plots.plot!(p1, w,
        label = "fitted-function",
        linewidth = 3.0,
        linestyle = :solid,
        alpha = 0.5)

    folderName = string(dirname(dirname(@__DIR__)))*"/processedData/"
    filename = folderName*"denoisedSignal_"*pr.alg.method*"_"*string(n)*".png"


    if savePlot
        
        if isfile(filename)
            rm(filename)
        end
    
        Plots.savefig(p2, filename)
    end
    
    return p2
end