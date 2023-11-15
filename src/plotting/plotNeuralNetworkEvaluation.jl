using Plots

function calculate_confusion_matrix(y_pred, y_true, threshold)

    @show length(y_pred)
    @show length(y_true)
    binary_predictions = Int.(y_pred .> threshold)
    binary_actual = Int.(y_true .> threshold)
    
    true_positives = sum((binary_predictions .== 1) .& (binary_actual .== 1))
    false_positives = sum((binary_predictions .== 1) .& (binary_actual .== 0))
    true_negatives = sum((binary_predictions .== 0) .& (binary_actual .== 0))
    false_negatives = sum((binary_predictions .== 0) .& (binary_actual .== 1))
    
    confusion_matrix = [true_negatives false_negatives;
                        false_positives true_positives]
    
    accuracy = (true_positives + true_negatives) / length(y_true)
    specificity = true_negatives / (true_negatives + false_positives)
    sensitivity = true_positives / (true_positives + false_negatives)
    
    return confusion_matrix, accuracy, specificity, sensitivity
end

function plot_and_save_heatmap(confusion_matrix, title, filename;
    savePlot=true,
    saveLocation="processedData/")

    heatmap_plot = Plots.heatmap(["Actual Negative", "Actual Positive"], 
        ["Predicted Negative", "Predicted Positive"], 
        confusion_matrix, color=:blues, aspect_ratio=1,
        title=title)

    if savePlot
        save_path = joinpath(saveLocation, filename)
        Plots.savefig(heatmap_plot, save_path)
    end

end

function check_training_accuracy(res;
    threshold=0.5,
    plotHeatMap=true, 
    savePlot=true,
    saveLocation="processedData/")

    pr = res.pr
    p = pr.p
    obj = pr.objective
    w_trained = res.xvals[:, end]
    classify = true

    params = p[:params]

    params[:classify] = classify

    y_pred_tr = vec(obj(w_trained, p))

    y_tr = params[:classData]

    confusion_matrix_tr, accuracy, specificity, sensitivity = calculate_confusion_matrix(y_pred_tr, y_tr, threshold)
    
    println("Training Data:
            Accuracy: $accuracy
            Specificity: $specificity
            Sensitivity: $sensitivity")

    if plotHeatMap
        plot_and_save_heatmap(confusion_matrix_tr, "Fit of FFNN on the Training Data",  "training_data_heatmap.png", savePlot=savePlot, saveLocation=saveLocation)
    end
end

function check_test_accuracy(res;
    threshold=0.5, 
    plotHeatMap=true, 
    savePlot=true,
    saveLocation="processedData/")

    pr = res.pr
    p = pr.p
    obj = pr.objective
    w_trained = res.xvals[:, end]
    classify = true

    params = p[:params]
    
    params[:classify] = classify
    testData = params[:testData]
    testClassData = params[:testClassData]
    params[:trainData] = testData
    params[:classData] = testClassData
    
    p[:params] = params

    y_pred_test = vec(obj(w_trained, p))
    y_test = params[:classData]

    confusion_matrix_test, accuracy_test, specificity_test, sensitivity_test = calculate_confusion_matrix(y_pred_test, y_test, threshold)
    
    println(
            "Testing Data:\n"*
            "Accuracy: $(accuracy_test)\n"*
            "Specificity: $(specificity_test)\n"
            *"Sensitivity: $(sensitivity_test)"
            )

    if plotHeatMap
        plot_and_save_heatmap(confusion_matrix_test, "Fit of FFNN on the Testing Data", "testing_data_heatmap.png", savePlot=savePlot, saveLocation=saveLocation)
    end
end