using CSV
using DataFrames

function convert_to_prf(filepath="profiles/M5BAT/Batt1.csv", target_folder="profiles/M5BAT/Batt1", timestep_col_name="DateAndTime")
    df = CSV.read(filepath, DataFrame, delim=";", decimal=',')
    df = df[!, ["DateAndTime", "P_AC", "Q_AC", "SOC", "I_DC_Batt", "U_DC_Batt"]]
    df.SOC = df.SOC * 0.1
    df.I_DC_Batt = df.I_DC_Batt .* 0.1
    df.U_DC_Batt = df.U_DC_Batt .* 0.1
    df.P_DC_Batt_out = ifelse.(df.I_DC_Batt .> 0, df.I_DC_Batt .* df.U_DC_Batt, 0)
    df.P_DC_Batt_in = ifelse.(df.I_DC_Batt .< 0, abs.(df.I_DC_Batt) .* df.U_DC_Batt, 0)
    mkpath(target_folder)
    # allowmissing!(df)
    # df[:, end-3:end] .= ifelse.(df[!, end-3:end] .< 0, missing, df[!, end-3:end])
    # df = df[:, 1:end-2]
    # df = Impute.interpolate(df)
    # timestep_size = df[5, timestep_col_name] - df[6, timestep_col_name]
    for col_name in names(df)
        if col_name != timestep_col_name 
        
            target_path = target_folder * "/" * col_name * ".prf" # TODO change

            header_variables = ["# data_type:", "# time_definition:",  
                                "# timestamp_format:", "# interpolation_type:", 
                                "# time_shift_seconds:"]
            if occursin("irr", col_name)
                header_values = ["extensive", "datestamp",
                                "yyyy-mm-dd HH:MM:SS", "linear_solar_radiation",
                                0] # TODO change
            else
                header_values = ["intensive", "datestamp",
                                "yyyy-mm-dd HH:MM:SS", "linear_classic",
                                0] # TODO change
            end

            header = DataFrame(timestep = header_variables, a = header_values)
            rename!(header, :a => col_name) 

            CSV.write(target_path, header, header=false, delim='\t')

            column_df = df[:, [timestep_col_name, col_name]]

            CSV.write(target_path, column_df, header=false, append=true, delim=';', decimal='.')
            print(target_path*"\n")
        end
    end
end

for (path, _, files) in walkdir("profiles/M5BAT/")
    for file in files
        if file != "BESS.csv"
            convert_to_prf(joinpath(path, file), joinpath(path, file)[1:end-4])
        end
    end
end