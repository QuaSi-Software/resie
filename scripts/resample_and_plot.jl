    # t1 = 0.0
    # t1 = @elapsed begin
    # end 
    # append!(unit.time_array_1, t1*10^6)
    # if sim_params["time"] > 24*365*4*900-7200
    #     unit.time_array_1 = unit.time_array_1[2:end]
    #     unit.time_array_2 = unit.time_array_2[2:end]
    #     unit.time_array_3 = unit.time_array_3[2:end]
    #     @info unit.time_array_1
    #     @info unit.time_array_2
    #     @info unit.time_array_3
    #     @info "Mean:$(round(mean(unit.time_array_1), digits=3)) us; std: $(round(std(unit.time_array_1), digits=3)) us; max:$(round(maximum(unit.time_array_1), digits=3)) us"
    #     @info "Mean:$(round(mean(unit.time_array_2), digits=3)) us; std: $(round(std(unit.time_array_2), digits=3)) us; max:$(round(maximum(unit.time_array_2), digits=3)) us"
    #     @info "Mean:$(round(mean(unit.time_array_3), digits=3)) us; std: $(round(std(unit.time_array_3), digits=3)) us; max:$(round(maximum(unit.time_array_3), digits=3)) us"
    # end
using CSV
using DataFrames
using Dates
using Statistics
using PlotlyJS
# plotlyjs()
using Infiltrator

function resample_hour_to_15min(df_hour::DataFrame; ref15=nothing, mode::String="repeat", timecol::Union{String,Nothing}=nothing)
    parse_dt(x) = begin
        if x === missing || x === nothing
            return nothing
        elseif x isa DateTime
            return x
        elseif x isa Date
            return DateTime(x)
        elseif x isa AbstractString
            for fmt in (dateformat"yyyy-mm-ddTHH:MM:SS", dateformat"yyyy-mm-dd HH:MM:SS", dateformat"yyyy/mm/dd HH:MM:SS", dateformat"dd.mm.yyyy HH:MM:SS")
                try
                    return DateTime(x, fmt)
                catch
                end
            end
            try
                return DateTime(x)
            catch
                return nothing
            end
        else
            return nothing
        end
    end

    n_hour = nrow(df_hour)
    if ref15 !== nothing
        n_ref = nrow(ref15)
        if n_hour * 4 != n_ref
            @warn "ref15 length does not equal 4 * nrow(df_hour): $(n_ref) != $(4*n_hour)"
        end
    end

    colnames = names(df_hour)

    # prepare output columns as Any arrays
    outcols = Dict{String, Vector{Any}}()
    for c in colnames
        outcols[c] = Vector{Any}()
    end

    for r in eachrow(df_hour)
        # prepare values (apply multiply if requested)
        vals = Dict{String,Any}()
        for c in colnames
            v = r[c]
            if mode == "multiply" && v isa Number && v !== missing
                v = v * 4
            end
            vals[c] = v
        end

        for k in 0:3
            for c in colnames
                v = vals[c]
                if timecol !== nothing && c == timecol
                    dt = parse_dt(v)
                    if dt !== nothing
                        push!(outcols[c], dt + Minute(15*k))
                    else
                        push!(outcols[c], v) # fallback: repeat raw value
                    end
                else
                    push!(outcols[c], v)
                end
            end
        end
    end

    return DataFrame(outcols)
end

function prepare_df(filepath_measure="C:/Users/demartin/Software/resie_battery/resie/profiles/M5BAT/Batt9_shift.csv", 
                    filepath_resie="C:/Users/demartin/Software/resie_battery/resie/output/out_LI_LFP_9_shift.csv",
                    filepath_sam="C:/Users/demartin/Software/SAM Nrel/results_m5bat_Li_LFP_9_5min_shift.csv")
    original_header_measure = ["DateAndTime", "P_AC_Set", "P_AC_Set", "P_AC", "Q_AC", "SOC", "I_DC_Batt", "U_DC_Batt",
                               "Mode_PQ", "Mode_Stop", "Mode_Silent", "Mode_Wait", "interpolated"]
    new_header_measure = ["date_raw", "P_AC_Set", "P_react_AC_Set", "P_AC", "P_react_AC", "SOC", "I_bat",
                          "V_bat", "Mode_PQ", "Mode_Stop", "Mode_Silent", "Mode_Wait", "interpolated"]
    original_header_resie = ["Time [dd.mm.yyyy HH:MM:SS]", "BYD_HVM_8.3 Capacity", "BYD_HVM_8.3 CellVoltage", "BYD_HVM_8.3 charge_efficiency",
                             "BYD_HVM_8.3 Cycles", "BYD_HVM_8.3 discharge_efficiency", "BYD_HVM_8.3 ExtractedCharge", "BYD_HVM_8.3 Load",
                             "BYD_HVM_8.3 Load%", "BYD_HVM_8.3 LossesGains", "BYD_HVM_8.3 m_e_ac_230v IN", "BYD_HVM_8.3 m_e_ac_230v OUT",
                             "BYD_HVM_8.3 SOC", "BYD_HVM_8.3 Temperature", "TST_BUS_EL_01 Balance", "TST_CH_01 m_e_ac_230v OUT",
                             "TST_CH_01 Supply", "TST_DEM_01 Demand", "TST_DEM_01 m_e_ac_230v IN", "TST_GRI_EL_01_out Input_sum",
                             "TST_GRI_EL_01_out m_e_ac_230v IN", "TST_GRI_EL_01 m_e_ac_230v OUT", "TST_GRI_EL_01 Output_sum"]
    new_header_resie = ["date", "E_max", "V_cell", "eta_ch", "n_efc", "eta_dis", "1-Q", "E", 
                        "E_%", "E_loss", "E_in_Wh", "E_out_Wh", "SOC", "Temp", 
                        "Balance", "E_ch", "E_ch_2", "E_dis", "E_dis_2", "E_grid_out_sum", "E_grid_out", "E_grid_in", "E_grid_in_sum"]
    original_header_sam = ["Time stamp", "Battery total charge | (Ah)", "Electricity to/from battery AC | (kW)", "Battery state of charge | (%)", 
                           "Battery number of cycles", "Battery maximum charge with degradation | (Ah)", "Battery maximum charge at temperature | (Ah)", 
                           "Battery cell voltage | (V)", "Battery current | (A)", "Battery voltage | (V)", "Battery relative capacity to nameplate | (%)", 
                           "Battery relative capacity to nameplate (calendar) | (%)", "Battery relative capacity to nameplate (cycling) | (%)",
                           "Battery temperature | (C)", "Battery average cycle DOD", "Battery cycle depth of discharge | (%)"]
    new_header_sam = ["timestamp", "Q", "P_kW", "SOC", "n_avg", "Q_max", "Q_max_temp", 
                      "V_cell", "I_bat", "V_bat", "Q_max_%", "Q_max_cal_%", "Q_max_cycl_%", "Temp",
                      "DOD_avg", "DOD_cycle"]
    df_measure = CSV.read(filepath_measure, DataFrame, header=new_header_measure, skipto=2, delim=';')
    df_resie = CSV.read(filepath_resie, DataFrame, header=new_header_resie, skipto=2, delim=';', decimal=',')

    df_measure.date_raw = DateTime.(df_measure.date_raw, "yyyy-mm-dd HH:MM:SS")
    df_resie.date = DateTime.(df_resie.date, "dd.mm.yyyy HH:MM:SS")

    timestep = Second(df_resie.date[2] - df_resie.date[1]).value

    df_resie.eta_ch = ifelse.(isnan.(df_resie.eta_ch), missing, df_resie.eta_ch)
    df_resie.eta_dis = ifelse.(isnan.(df_resie.eta_dis), missing, df_resie.eta_dis)
    df_resie.SOC_adj = (df_resie.SOC .- 0) ./ 1.0
    df_resie.V_bat = df_resie.V_cell .* 240

    df_measure.SOC = df_measure.SOC .* 0.1
    df_measure.SOC_adj = df_measure.SOC
    df_measure.V_bat = df_measure.V_bat .* 0.1
    df_measure.V_cell = df_measure.V_bat ./ 240
    # filter values with small I_bat to get more consistent data
    I_min = 0
    df_measure.I_bat = ifelse.(abs.(df_measure.I_bat) .> I_min .* 10, df_measure.I_bat .* 0.1, 0)
    df_measure.P_AC = ifelse.(abs.(df_measure.I_bat) .> I_min, df_measure.P_AC .* 1000, 0)
    df_measure.P_react_AC = ifelse.(abs.(df_measure.I_bat) .> I_min, df_measure.P_react_AC .* 1000, 0)
    
    df_measure.P_AC_Set = df_measure.P_AC_Set .* 1000
    df_measure.P_react_AC_Set = df_measure.P_react_AC_Set .* 1000
    df_measure.P_bat = df_measure.I_bat .* df_measure.V_bat
    df_measure.S_AC = sqrt.(df_measure.P_AC .^ 2 + df_measure.P_react_AC .^ 2)
    df_measure.E_in_Wh = ifelse.(df_measure.P_bat .< 0, .-df_measure.P_bat .* (1/3600) , 0)
    df_measure.E_out_Wh = ifelse.(df_measure.P_bat .> 0, df_measure.P_bat .* (1/3600) , 0)
    
    df_measure.n_efc = (df_measure.E_in_Wh .+ df_measure.E_out_Wh) ./ (df_resie.E_max[1] * 2)
    df_measure.eta_ch_inv = ifelse.(df_measure.P_bat .< 0 .&& df_measure.S_AC .> 0, .-df_measure.P_bat ./ df_measure.S_AC, missing)
    df_measure.eta_dis_inv = ifelse.(df_measure.P_bat .> 0 .&& df_measure.S_AC .> 0, df_measure.S_AC ./ df_measure.P_bat, missing)
   
    # aggregate measurements to timestep of resie
    filter_cols = ["SOC_adj", "SOC", "V_bat", "V_cell"]
    mean_cols = ["P_AC_Set", "P_react_AC_Set", "P_AC", "P_react_AC", "I_bat", "P_bat", "eta_ch_inv", "eta_dis_inv", "S_AC"]
    sum_cols = ["E_in_Wh", "E_out_Wh", "n_efc"]

    df_measure = df_measure[(df_measure.date_raw .>= df_resie.date[1]) .& (df_measure.date_raw .<= df_resie.date[end]), :]

    # create time bins according to df_resie
    df_measure.date = floor.(df_measure.date_raw, Second(timestep))
    # create groups with each group containing all rows for one time bin
    grouped = groupby(df_measure, :date)
    # combine returns one row per group by using the pattern:
    # source_columns .=> function .=> new_column_names
    df_measure_new = combine(grouped,
                             mean_cols .=> mean ∘ skipmissing,
                             sum_cols  .=> mean ∘ skipmissing,
                             filter_cols .=> last,
                             renamecols=false
                            )

    # measure_filter = Array{Union{Float64, DateTime}}(undef, size(df_resie, 1), size(filter_cols, 1))
    # measure_mean = Array{Union{Float64, Missing}}(undef, size(df_resie, 1), size(mean_cols, 1))
    # measure_sum = Array{Union{Float64, Missing}}(undef, size(df_resie, 1), size(sum_cols, 1))

    # for i in 1:nrow(df_resie)
    #     if i == size(df_resie, 1)
    #         filtered_rows = df_measure[(df_measure.date .>= df_resie.date[i]) .& (df_measure.date .< df_resie.date[i] + Second(timestep)), :]
    #     else
    #         filtered_rows = df_measure[(df_measure.date .>= df_resie.date[i]) .& (df_measure.date .< df_resie.date[i+1]), :]
    #     end
    #     means = map(col -> mean(skipmissing(filtered_rows[!, col])), mean_cols)
    #     measure_mean[i, :] = means'
    #     sums = map(col -> sum(skipmissing(filtered_rows[!, col])), sum_cols)
    #     measure_sum[i, :] = sums'
    #     measure_filter[i, :] = Vector(filtered_rows[end, filter_cols])     
    # end
    # df_measure_new = DataFrame(hcat(df_resie.date, measure_filter, measure_mean, measure_sum), ["date"; filter_cols; mean_cols; sum_cols])

    df_measure_new.n_efc = cumsum(df_measure_new.n_efc)

    energy_diff_start_end = abs(df_measure.SOC[1] - df_measure.SOC[end]) / 100 * df_resie.E_max[1]
    round_trip_efficiency = (sum(df_measure_new.E_out_Wh) + energy_diff_start_end/2) / (sum(df_measure_new.E_in_Wh) + energy_diff_start_end/2)

    E_AC_in_Wh = sum(ifelse.(df_measure_new.P_bat .< 0 .&& df_measure_new.S_AC .> 0, df_measure_new.P_bat .* -1000 .* (timestep/3600), 0))
    E_AC_out_Wh = sum(ifelse.(df_measure_new.P_bat .> 0 .&& df_measure_new.S_AC .> 0, df_measure_new.P_bat .* 1000 .* (timestep/3600), 0))
    
    df_resie.P_bat = (df_resie.E_out_Wh .- df_resie.E_in_Wh) ./ (timestep/3600)

    if filepath_sam != nothing
        df_sam = CSV.read(filepath_sam, DataFrame, header=new_header_sam, skipto=2)

        year_resie = Year(df_resie.date[1]).value
        start_date = DateTime(year_resie, 01, 01, 0, 0, 0)
        end_date = start_date + Second(timestep*(size(df_sam, 1)-1))
        df_sam.date = collect(start_date:Second(timestep):end_date)
        df_sam = df_sam[df_resie.date[1] .<= df_sam.date .<= df_resie.date[end], :]
        df_sam = df_sam[(df_sam.date .== df_resie.date), :]

        df_sam.P_bat = df_sam.P_kW .* 1000
        df_sam.SOC_adj = df_sam.SOC
        df_sam.E_in_Wh = ifelse.(df_sam.P_bat .< 0, .-df_sam.P_bat .* (timestep/3600) , 0)
        df_sam.E_out_Wh = ifelse.(df_sam.P_bat .> 0, df_sam.P_bat .* (timestep/3600) , 0)
        df_sam.n_efc = df_sam.n_avg .* df_sam.DOD_avg ./ 100

        df_resie.Q = df_sam.Q_max .- df_resie[!, "1-Q"] .* 100

    else
        df_sam = nothing
    end

    return df_measure_new, df_resie, df_sam
end

function plot_df(df_measure::Union{DataFrame, Nothing}, 
                 df_resie::Union{DataFrame, Nothing}, 
                 df_sam::Union{DataFrame, Nothing}, 
                 col_name_1, col_name_2::Union{String, Nothing}=nothing)
    traces = GenericTrace[]
    if !isnothing(df_measure)
        trace_m_1 = scatter(x=df_measure.date, y=df_measure[:, col_name_1], 
                        mode="lines", name="$col_name_1 Measurement", line=attr(width=1))
        push!(traces, trace_m_1)
    end
    if !isnothing(df_resie)
        trace_r_1 = scatter(x=df_resie.date, y=df_resie[:, col_name_1], 
                        mode="lines", name="$col_name_1 ReSiE", line=attr(width=1))
        push!(traces, trace_r_1)
    end
    if !isnothing(df_sam)
        trace_s_1 = scatter(x=df_sam.date, y=df_sam[:, col_name_1], 
                        mode="lines", name="$col_name_1 SAM", line=attr(width=1))
        push!(traces, trace_s_1)
    end
    if !isnothing(col_name_2)
        if !isnothing(df_measure)
            trace_m_2 = scatter(x=df_measure.date, y=df_measure[:, col_name_2], 
                        mode="lines", name="$col_name_2 Measurement", line=attr(width=1), yaxis="y2")
            push!(traces, trace_m_2)
        end
        if !isnothing(df_resie)
            trace_r_2 = scatter(x=df_resie.date, y=df_resie[:, col_name_2], 
                            mode="lines", name="$col_name_2 ReSiE", line=attr(width=1), yaxis="y2")
            push!(traces, trace_r_2)
        end
        if !isnothing(df_sam)
            trace_s_2 = scatter(x=df_sam.date, y=df_sam[:, col_name_2], 
                            mode="lines", name="$col_name_2 SAM", line=attr(width=1), yaxis="y2")
            push!(traces, trace_s_2)
        end
    end
    layout = Layout(template=:plotly_white, 
                    yaxis_title_text=col_name_1, 
                    yaxis2=attr(side="right", overlaying="y", title=col_name_2))
    plot(traces, layout, config=PlotConfig(scrollZoom=true))
end
