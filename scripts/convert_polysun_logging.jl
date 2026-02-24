using CSV
using DataFrames
# using TimeSeries
# using TimeSeriesResampler: resample, mean, ohlc, sum, TimeFrame

filepath = "profiles/Timesteps_variation_1.csv"
dt = 3600

df = CSV.read(filepath, DataFrame, delim='\t', decimal=',')
df = df[!, Not(end)]
df[!, 9:end] = df[!, "Zeitschritt (s)"] .* df[!, 9:end]
relevant_cols = names(df[:, 9:end])
# number_of_relevant_cols = size(relevant_cols)[1]

df_dt = DataFrame(timestep = [0:dt:31536000;])
for col_name in relevant_cols
    df_dt[!, col_name] = missings(Float64, nrow(df_dt))
end
mat_df = Matrix(df)
for (idx, timestep) in enumerate([0:dt:31536000-dt;])
    if idx % 1000 == 0
        println(idx)
    end

    t1 = timestep
    t2 = df_dt[idx+1, :timestep]

    upper_h_idx = mat_df[:, 3] .== div(t1,3600)*3600
    df_dt[idx, 2:7] = Matrix(df[upper_h_idx, 9:14] ./ df[upper_h_idx, 2])

    mask = mat_df[:, 3] .>= t1 .&& mat_df[:, 3] .< t2
    df_dt[idx, 8:end] = sum.(eachcol(df[mask, 15:end])) ./ dt
end


# df = df[!, 8:end]
new_filepath = filepath[1:end-4]*"_"*string(dt)*".csv"
try
    CSV.write(new_filepath, df_dt[1:end-1, :], delim=";", decimal=',')
catch
    CSV.write(new_filepath[1:end-4]*"new.csv", delim=";", decimal=',')
end
