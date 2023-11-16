function preprocessLiverData()

    rawDataFolder = "rawData/";
    filename = rawDataFolder * "Indian Liver Patient Dataset (ILPD).csv";

    header = ["Age", "Gender", "TB", "DB", "Alkphos", "Sgpt", "Sgot", "TP", "ALB", "A/G Ratio", "Diagnosis"];

    df00 = CSV.File(filename, header=header) |> DataFrame;
    df0 = dropmissing(df00);
    df = deepcopy(df0);
    # Convert "Gender" to 0 and 1
    df.Gender = map(gender -> gender == "Male" ? 0 : 1, df.Gender);
    # Convert "Diagnosis" from 1 to 1/3 and 2 to 2/3
    df.Diagnosis = map(diagnosis -> diagnosis == 1 ? 1/3 : 2/3, df.Diagnosis);

    # Normalize the first 10 columns
    for col in names(df)[1:10]
        min_val, max_val = minimum(df[!, col]), maximum(df[!, col]);
        df[!, col] = map(x -> (x - min_val) / (max_val - min_val), df[!, col]);
    end

    select!(df, Not(["ALB", "A/G Ratio"]))

    return df
    
end