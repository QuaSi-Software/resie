using CSV
using DataFrames
using PlotlyJS

file = "./output/results_price_optim.csv"
name = split(file, '/')[end][1:end-4]
volume_col = "components Hafner_Puffer_gross volume"

df = CSV.read(file, DataFrame; delim=';', decimal=',')

ei_idx = findfirst(==("error"), names(df))
cols = names(df)[ei_idx+1:end]

df_sorted = sort(df, volume_col)

plot_df = df_sorted[:, vcat(volume_col, cols)]
plot_df = select(plot_df, Not(["hours_to_full", "hours_to_empty", "Hafner_BUS_Waerme Balance", "objective"]))
new_cols = ["Speichervolumen / m³", "Energie Beladung Speicher", "Energie Entladung Speicher", "Speicher Verluste",
            "Eigennutzungsgrad", "Eigenversorgungsgrad", "Kosten Strom Stromnetz->WP / €", "Kosten Strom PV->WP / €", 
            "Kosten Strom WP gesamt / €"]
rename!(plot_df, new_cols)
# p = plot(plot_df.Volume, Matrix(plot_df[!, 2:end]); names=permutedims(names(plot_df)))
traces = [
        scatter(
            x = plot_df[!, "Speichervolumen / m³"],
            y = plot_df[!, v],
            mode = "lines",
            name = v
        )
        for v in names(plot_df)[2:end]
    ]
layout = Layout(
        title = name,
        xaxis_title = "Speichervolumen / m³",
        yaxis_title = "Value"
    )
p = plot(traces, layout)

# Save to HTML
savefig(p, "./output/plot_" * name * ".html")