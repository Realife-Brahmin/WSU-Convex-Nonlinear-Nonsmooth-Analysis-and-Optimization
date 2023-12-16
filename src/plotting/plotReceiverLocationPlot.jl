using Plots

function plotReceiverLocationPlot(
    res; 
    savePlot::Bool = true)

    Plots.theme(:dao)
    pr = res.pr
    p = pr[:p]
    params = p[:params]
    
    @unpack pmean, P, μ, datasetName = params

    folderName = string(dirname(dirname(@__DIR__)))*"/processedData/"
    filename = folderName*"receiverLocation_"*datasetName*"_"*pr.alg.method*
    "_n_"*string(n)*"_m_"*string(m)*".png"

    if savePlot
        
        if isfile(filename)
            rm(filename)
        end
    
        Plots.savefig(h2, filename)
    end

    return h2

end

