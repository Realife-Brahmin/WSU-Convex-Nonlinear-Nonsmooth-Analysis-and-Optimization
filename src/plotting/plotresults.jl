include("plot_fval_vs_iterations.jl")
include("plotDragCurve.jl")
include("plotDenoisedSignal.jl")
include("plotReceiverLocationPlot.jl")
include("plotNeuralNetworkEvaluation.jl")
include("plotMinPathTimeTrajectory.jl")
include("plotTransmissionLines.jl")

function plotresults(res;
    savePlot::Bool=true)
    
    plot_fval_vs_iterations(res, savePlot=savePlot)

    if typeof(res) <: Vector
        functionName = res[1].pr[:objectiveString]
    else
        functionName = res.pr[:objectiveString]
    end

    if functionName == "drag"
        plotDragCurve(res, savePlot=savePlot)
    # elseif functionName == "receiverLocation"
    #     plotReceiverLocationPlot(res, savePlot=savePlot)
    elseif functionName == "signalDenoise"
        plotDenoisedSignal(res, savePlot=savePlot)
    elseif functionName == "nnloss"
        check_training_accuracy(res, savePlot=savePlot)
        check_test_accuracy(res, savePlot=savePlot)
    elseif functionName == "pathtime"
        plotMinPathTimeTrajectory(res, savePlot=savePlot,
        plotTrajectory=false)
    elseif functionName == "transmissionLines"
        plotTransmissionLines(res, savePlot=savePlot)
    else
        println("Nothing to plot.")
    end
end