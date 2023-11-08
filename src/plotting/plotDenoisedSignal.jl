using Plots

function plotDenoisedSignal(res;
    savePlot::Bool=true)

    pr = res.pr
    d = pr.x0
    n = length(d)
    w = res.xvals[:, end]
    t = collect(1:1:n)
    
    theme(:dao)

    p1 = Plots.scatter(t, d,
        title = "Denoised Signal",
        label = "insertTwoLabels",
        xlabel = "instance [unit]",
        ylabel = "signal value")

    p2 = Plots.plot!(p1, w)

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