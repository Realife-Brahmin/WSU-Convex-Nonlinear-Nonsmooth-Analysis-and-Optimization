include("plotDragCurve.jl")

function plotresults(res;
    savePlot::Bool=true)
    pr = res.pr
    functionName = string(pr.objective)
    if functionName == "drag"
        plotDragCurve(res, savePlot=savePlot)
    elseif functionName == "signalDenoise"
        println("Plotting under construction.")
    else
        println("Nothing to plot.")
    end
end