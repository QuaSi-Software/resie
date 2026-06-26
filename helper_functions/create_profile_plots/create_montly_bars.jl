using CSV
using DataFrames
using Dates
using StatsPlots

csv_path_100 = "./output/csv_price_optim/out_Hafner_Puffer_gross-initial_load-0.5_Hafner_Puffer_gross-volume-100_Hafner_Puffer_gross-h_to_r-2.25.csv"
csv_path_200 = "./output/csv_price_optim/out_Hafner_Puffer_gross-initial_load-0.25_Hafner_Puffer_gross-volume-200_Hafner_Puffer_gross-h_to_r-2.25.csv"
csv_path_500 = "./output/csv_price_optim/out_Hafner_Puffer_gross-initial_load-0.1_Hafner_Puffer_gross-volume-500_Hafner_Puffer_gross-h_to_r-5.06.csv"
csv_path_800 = "./output/csv_price_optim/out_Hafner_Puffer_gross-initial_load-0.0625_Hafner_Puffer_gross-volume-800_Hafner_Puffer_gross-h_to_r-4.0.csv"
csv_path_1000 = "./output/csv_price_optim/out_Hafner_Puffer_gross-initial_load-0.05_Hafner_Puffer_gross-volume-1000_Hafner_Puffer_gross-h_to_r-3.58.csv"

csv_paths = [csv_path_100, csv_path_200, csv_path_500, csv_path_800, csv_path_1000]

cols_electrictiy = ["Hafner_PV_Freiflaeche m_e_ac_230v OUT", "m_e_ac_230v EnergyFlow Hafner_PV_Freiflaeche->Hafner_WP", "m_e_ac_230v EnergyFlow Hafner_PV_Freiflaeche->Hafner_Stromnetz_OUT", "m_e_ac_230v EnergyFlow Hafner_Stromnetz_IN->Hafner_WP"]
cols_heat = ["m_heat_WP EnergyFlow Hafner_WP->Hafner_DEM_Waerme", "secondary_m_heat_WP EnergyFlow secondary_Hafner_WP->Hafner_DEM_Waerme", "m_heat_WP EnergyFlow Hafner_Puffer_gross->Hafner_DEM_Waerme"]
cols_cost = ["energy_cost_total", "energy_cost_grid", "energy_cost_pv"]
cols_losses = ["Hafner_Puffer_gross LossesGains"]

value_cols = vcat(cols_electrictiy, cols_heat, cols_cost, cols_losses)
# value_cols = cols_cost

for csv_path in csv_paths
    df = CSV.read(csv_path, DataFrame; delim=';', decimal=',')
    df.timestamp = DateTime.(df[!, 1], dateformat"dd.mm.yyyy HH:MM:SS")
    df.month = month.(df.timestamp)
    df.energy_cost_grid = df[!, "m_e_ac_230v EnergyFlow Hafner_Stromnetz_IN->Hafner_WP"] .* df[!, "Grid_Price_IN Temperature"] ./ 10^6
    df.energy_cost_pv = df[!, "m_e_ac_230v EnergyFlow Hafner_PV_Freiflaeche->Hafner_WP"] .* 0.08 ./ 10^3
    df.energy_cost_total = df.energy_cost_grid .+ df.energy_cost_pv
    # --- Aggregate monthly sums ---
    monthly = combine(groupby(df, [:month]),
                                [col => sum => col for col in value_cols]...)
    volume = split(split(csv_path, "volume-")[2], "_Hafner")[1]
    monthly.volume .= volume
    new_cols_electricity = ["Energie Strom Erzeugung PV / Wh", "Energie Strom PV->WP / Wh", "Energie Strom PV->Stromnetz / Wh", "Energie Strom Stromnetz->WP / Wh"]
    new_cols_heat = ["Energie Wärme WP_PV->Bedarf / Wh", "Energie Wärme WP_Stromnetz->Bedarf / Wh", "Energie Wärme Speicher->Bedarf / Wh"]
    new_cols_cost = ["Kosten Strom WP gesamt / €", "Kosten Strom Stromnetz->WP / €", "Kosten Strom PV->WP / €"]
    new_cols_losses = ["Energie Wärmeverluste Speicher"]
    old_new_cols = ["Speichervolumen / m³"]
    new_cols = vcat(["Monat"], new_cols_electricity, new_cols_heat, new_cols_cost, new_cols_losses, old_new_cols)
    rename!(monthly, new_cols)
    CSV.write("./output/csv_price_optim/" * "monthly_sums_" * volume * ".csv", monthly)
end 

name = split(csv_paths[1], '/')[3]
Y = []
df_all = DataFrame()
for csv_path in csv_paths
    df = CSV.read(csv_path, DataFrame; delim=';', decimal=',')

    # --- Parse time column ---
    df.timestamp = DateTime.(df[!, 1], dateformat"dd.mm.yyyy HH:MM:SS")
    df.month = month.(df.timestamp)

    # calculate cost columns
    df.energy_cost_grid = df[!, "m_e_ac_230v EnergyFlow Hafner_Stromnetz_IN->Hafner_WP"] .* df[!, "Grid_Price_IN Temperature"] ./ 10^6
    df.energy_cost_pv = df[!, "m_e_ac_230v EnergyFlow Hafner_PV_Freiflaeche->Hafner_WP"] .* 0.08 ./ 10^3
    df.energy_cost_total = df.energy_cost_grid .+ df.energy_cost_pv

    # --- Aggregate monthly sums ---
    monthly = combine(groupby(df, [:month]),
                                [col => sum => Symbol(col * "_sum") for col in value_cols]...)
    monthly.volume .= split(split(csv_path, "volume-")[2], "_Hafner")[1]
    
    sort!(monthly, :month)
    monthly.month_label = Dates.format.(Date.(2025, monthly.month, 1), "u")
    if isempty(df_all)
        global df_all = monthly
    else
        global df_all = vcat(df_all, monthly)
    end
    if isempty(Y)
        global Y = monthly[!, value_cols[1] * "_sum"]
    else
        global Y = hcat(Y, monthly[!, value_cols[1] * "_sum"])
    end
end

names = hcat(unique(df_all.volume)...)
# --- Plot (grouped via offsets) ---
x = 1:size(Y, 1)
n = size(Y, 2)                       # number of series (5)
bw = min(0.8 / n, 0.28)              # bar width per series
offsets = ((1:n) .- (n + 1) / 2) .* (bw * 1.15)  # symmetric offsets

# first series
p = bar(x .+ offsets[1], Y[:, 1];
        bar_width = bw, label = names[1])

for j in 2:n
    bar!(x .+ offsets[j], Y[:, j];
         bar_width = bw, label = names[j])
end

plot!(p;
      xticks=(x, unique(df_all.month_label)), xrotation=0,
      legend=:top,
      xlabel="Month", ylabel=value_cols[1])

savefig(p, "./output/bar_" * name * "_" * value_cols[1] * ".png")

# sums_matrix = reduce(hcat, [monthly[!, Symbol(col * "_sum")] for col in value_cols])

# p = groupedbar(
#     monthly.month_label,
#     sums_matrix;
#     label = hcat(value_cols...),
#     rotation = 45,
#     bar_position = :stack,
#     legend = :top
# )

# @df df_all groupedbar(
#         :month;
#         bar_position = :dodge,
#         group = :volume,
#         legend = :top
#     )

