using CSV
using DataFrames
filepath = raw"c:\Users\jenter\Documents\resie\profiles\MA\csv\cbmp_up_max_15min.csv" # TODO change
timestep_col_name = "timestep" # TODO change
profile_start_date = "01.01.2024 00:00" # TODO change

df = CSV.read(filepath, DataFrame, delim=";", decimal=',')
# allowmissing!(df)
# df[:, end-3:end] .= ifelse.(df[!, end-3:end] .< 0, missing, df[!, end-3:end])
# df = df[:, 1:end-2]
# df = Impute.interpolate(df)
# df[end:end+20, :] = df[1:20, :]
# timestep_size = df[5, timestep_col_name] - df[6, timestep_col_name]
for col_name in names(df)
    if col_name != timestep_col_name 
    
        target_path = raw"c:\Users\jenter\Documents\resie\profiles\MA\\" * col_name * ".prf" # TODO change

        header_variables = ["# data_type:", "# time_definition:", "# timestamp_format:" 
                            # "# profile_start_date_format:", "# timestamp_format:", "# interpolation_type:"
                            ]
        if occursin("irr", col_name) # search in column name for certain word to filter
            header_values = ["intensive", "datestamp", "dd.mm.yyyy HH:MM" 
                             # "dd.mm.yyyy HH:MM", "dd.mm.yyyy HH:MM", "stepwise"
                             ] # TODO change
        else
            header_values = ["intensive", "datestamp", "dd.mm.yyyy HH:MM" 
                             # "dd.mm.yyyy HH:MM", "dd.mm.yyyy HH:MM", "stepwise"
                             ] # TODO change
        end

        header = DataFrame(timestep = header_variables, a = header_values)
        rename!(header, :a => col_name) 

        CSV.write(target_path, header, header=false, delim='\t')

        column_df = df[:, ["timestep",col_name]]

        CSV.write(target_path, column_df, header=false, append=true, delim=';', decimal='.')
        print(target_path*"\n")
    end
end