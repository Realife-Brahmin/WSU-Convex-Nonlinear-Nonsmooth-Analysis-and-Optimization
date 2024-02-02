function plot_fval_vs_iterations(res;
    savePlot=true)
    
        pr = res.pr
        functionName = string(pr.objective)
        theme(:dao)
        fvals = res.fvals
        f = res.fvals[end]
        x = res.xvals[:, end]
        n = length(x)
        itr = length(fvals)
    
        gr()
        p1 = Plots.plot(1:itr, fvals,
            title = "Function Value vs Iteration: $(L"f_{min} =") $(round(f, digits=5)) \n using $(pr.alg[:method])",
            titlefont = font(12,"Computer Modern"),
            guidefont = font(15,"Computer Modern"),
            label = L"f",
            ylabel = L"f",
            xlabel = L"itr",
            linewidth = 3,
            markersize = 2,
            alpha = 1.0,
        )
    
        folderName = string(dirname(dirname(@__DIR__)))*"/processedData/"
        filename = folderName*functionName*"_fval_vs_itr_"*pr.alg[:method]*"_"*string(n)*".png"
    
    
        if savePlot
            
            if isfile(filename)
                rm(filename)
            end
        
            Plots.savefig(p1, filename)
        end
    
        return p1
end