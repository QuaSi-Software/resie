using CSV
using DataFrames
filepath = raw"c:\Users\jenter\Documents\Resie\profiles\Profile Rohdaten\Verbrauchsprofile\2025-11-28_MB_Germersheim_Alle-Gebaude_berta-rudi_power.csv" # TODO change
timestep_col_name = "timestep" # TODO change
profile_start_date = "01.01.2024 00:00" # TODO change

df = CSV.read(filepath, DataFrame, delim=";", decimal=',')
# allowmissing!(df)
# df[:, end-3:end] .= ifelse.(df[!, end-3:end] .< 0, missing, df[!, end-3:end])
# df = df[:, 1:end-2]
# df = Impute.interpolate(df)
timestep_size = df[5, timestep_col_name] - df[6, timestep_col_name]
for col_name in names(df)
    if col_name != timestep_col_name 
    
        target_path = raw"c:\Users\jenter\Documents\Resie\profiles\MA\\" * col_name * ".prf" # TODO change

        header_variables = ["# data_type:", "# time_definition:", "# profile_start_date:", 
                            "# profile_start_date_format:", "# timestamp_format:", "# interpolation_type:", 
                            "# time_shift_seconds:"]
        if occursin("irr", col_name) # search in column name for certain word to filter
            header_values = ["intensive", "startdate_timestamp", profile_start_date, 
                            "dd.mm.yyyy HH:MM", "seconds", "linear_solar_radiation", 0] # TODO change
        else
            header_values = ["intensive", "startdate_timestamp", profile_start_date, 
                             "dd.mm.yyyy HH:MM", "seconds", "linear_classic", 0] # TODO change
        end

        header = DataFrame(timestep = header_variables, a = header_values)
        rename!(header, :a => col_name) 

        CSV.write(target_path, header, header=false, delim='\t')

        column_df = df[:, ["timestep",col_name]]

        CSV.write(target_path, column_df, header=false, append=true, delim=';', decimal='.')
        print(target_path*"\n")
    end
end