using CSV
using DataFrames
using Dates
using Statistics

# --- CONFIG ---
infile = raw"c:/Users/jenter/Documents/Backup Masterarbeit/CBMP/Auswertungen/aFRR_prices_2024_year_no_leapday.csv"
outdir = raw"c:/Users/jenter/Documents/Backup Masterarbeit/CBMP/Auswertungen"

year_dt0      = DateTime(2024, 1, 1, 0, 0, 0)
sec_per_bin   = 900
bins_per_year = 365 * 24 * 4  # 35040

out_up_max   = joinpath(outdir, "cbmp_up_max_15min_2024.csv")
out_up_mean  = joinpath(outdir, "cbmp_up_mean_15min_2024.csv")
out_dn_ext   = joinpath(outdir, "cbmp_down_extreme_15min_2024.csv")
out_dn_mean  = joinpath(outdir, "cbmp_down_mean_15min_2024.csv")
out_flags    = joinpath(outdir, "cbmp_activation_flags_15min_2024.csv")

# --- helper: find columns robustly ---
canonical(s) = lowercase(replace(String(s), r"_+$" => ""))

function findcol(df::DataFrame, candidate::String)
    nms = names(df)
    cn  = canonical.(nms)
    idx = findfirst(==(canonical(candidate)), cn)
    idx === nothing && error("Column '$candidate' not found. Found: $(String.(nms))")
    nm = nms[idx]
    return nm isa Symbol ? nm : Symbol(nm)
end

function main()
    # Read CSV; empty strings -> missing
    df = DataFrame(CSV.File(infile; missingstring="", normalizenames=true))
    rename!(df, Symbol.(names(df)))  # ensure all names are Symbols

    col_ts = findcol(df, "timestamp")
    col_up = findcol(df, "Price_Up_EUR_MWh_")
    col_dn = findcol(df, "Price_Down_EUR_MWh_")

    # Ensure numeric timestamp (seconds since start)
    df[!, col_ts] = Float64.(df[!, col_ts])

    # 15-min bin index (0..35039)
    df[!, :bin] = Int.(floor.(df[!, col_ts] ./ sec_per_bin))
    df = df[(df.bin .>= 0) .& (df.bin .< bins_per_year), :]

    upv = df[!, col_up]
    dnv = df[!, col_dn]

    # "no call" ONLY if both present AND exactly equal
    no_call = (.!ismissing.(upv)) .& (.!ismissing.(dnv)) .& (upv .== dnv)

    # call masks per direction (counts also rows where both exist but differ)
    is_up_call = (.!no_call) .& (.!ismissing.(upv))
    is_dn_call = (.!no_call) .& (.!ismissing.(dnv))

    # Build call-only frames
    df_up = df[is_up_call, [:bin, col_up]]
    df_dn = df[is_dn_call, [:bin, col_dn]]

    # Aggregations for bins where calls exist
    up_agg = combine(groupby(df_up, :bin),
        col_up => (x -> maximum(skipmissing(x))) => :cbmp_up_max,
        col_up => (x -> mean(collect(skipmissing(x)))) => :cbmp_up_mean
    )

    # DOWN "extreme": most negative value (minimum) within the 15-min bin
    dn_agg = combine(groupby(df_dn, :bin),
        col_dn => (x -> minimum(skipmissing(x))) => :cbmp_down_extreme,
        col_dn => (x -> mean(collect(skipmissing(x)))) => :cbmp_down_mean
    )

    # Activation flags per bin (1 if any call exists in that bin)
    up_flag = combine(groupby(df_up, :bin), nrow => :n_up)
    up_flag[!, :activation_up] = Int.(up_flag.n_up .> 0)
    select!(up_flag, [:bin, :activation_up])

    dn_flag = combine(groupby(df_dn, :bin), nrow => :n_dn)
    dn_flag[!, :activation_down] = Int.(dn_flag.n_dn .> 0)
    select!(dn_flag, [:bin, :activation_down])

    # Full year timeline (35040 bins)
    bins = collect(0:(bins_per_year-1))
    dt15 = year_dt0 .+ Minute.(15 .* bins)
    ts_str = Dates.format.(dt15, dateformat"dd/mm/yyyy HH:MM")

    full = DataFrame(bin=bins, timestamp=ts_str)

    full = leftjoin(full, up_agg, on=:bin)
    full = leftjoin(full, dn_agg, on=:bin)
    full = leftjoin(full, up_flag, on=:bin)
    full = leftjoin(full, dn_flag, on=:bin)

    # Fill missing prices with baseline = 0
    full[!, :cbmp_up_max]       = coalesce.(full.cbmp_up_max, 0.0)
    full[!, :cbmp_up_mean]      = coalesce.(full.cbmp_up_mean, 0.0)
    full[!, :cbmp_down_extreme] = coalesce.(full.cbmp_down_extreme, 0.0)
    full[!, :cbmp_down_mean]    = coalesce.(full.cbmp_down_mean, 0.0)

    # Fill missing flags with 0
    full[!, :activation_up]   = coalesce.(full.activation_up, 0)
    full[!, :activation_down] = coalesce.(full.activation_down, 0)

    # Write 4 price series
    CSV.write(out_up_max,  select(full, :timestep, :cbmp_up_max))
    CSV.write(out_up_mean, select(full, :timestep, :cbmp_up_mean))
    CSV.write(out_dn_ext,  select(full, :timestep, :cbmp_down_extreme))
    CSV.write(out_dn_mean, select(full, :timestep, :cbmp_down_mean))

    # Write flags (5th file)
    CSV.write(out_flags, select(full, :timestep, :activation_up, :activation_down))

    println("Done. Files written to: ", outdir)
    println("Rows per output (should be 35040): ", nrow(full))
    println("Activated UP bins: ", sum(full.activation_up))
    println("Activated DOWN bins: ", sum(full.activation_down))
end

main()
