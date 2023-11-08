include("plotDragCurve.jl")
include("plotDenoisedSignal.jl")

function plotresults(res;
    savePlot::Bool=true)
    pr = res.pr
    functionName = string(pr.objective)
    if functionName == "drag"
        plotDragCurve(res, savePlot=savePlot)
    elseif functionName == "signalDenoise"
        @error "Plotting under construction."
        plotDenoisedSignal(res, savePlot=savePlot)
    else
        println("Nothing to plot.")
    end
end