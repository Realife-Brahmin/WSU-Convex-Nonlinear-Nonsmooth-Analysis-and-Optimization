include("plotDragCurve.jl")
include("plotDenoisedSignal.jl")
include("plotNeuralNetworkEvaluation.jl")

function plotresults(res;
    savePlot::Bool=true)
    
    pr = res.pr
    functionName = string(pr.objective)
    if functionName == "drag"
        plotDragCurve(res, savePlot=savePlot)
    elseif functionName == "signalDenoise"
        plotDenoisedSignal(res, savePlot=savePlot)
    elseif functionName == "nnloss"
        check_training_accuracy(res, savePlot=savePlot)
        check_test_accuracy(res, savePlot=savePlot)
    else
        println("Nothing to plot.")
    end
end