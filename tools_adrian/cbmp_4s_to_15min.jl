using CSV
using DataFrames
using Dates
using Statistics

# --- CONFIG ---
infile = raw"c:/Users/jenter/Documents/Backup Masterarbeit/CBMP/Auswertungen/aFRR_prices_2024_year_no_leapday.csv"
outdir = raw"c:/Users/jenter/Documents/Backup Masterarbeit/CBMP/Auswertungen"

year_dt0      = DateTime(2024, 1, 1, 0, 0, 0)
sec_per_bin   = 900                       # 15 min
bins_per_year = 365 * 24 * 4              # 35040 (no leap day)

out_up_max  = joinpath(outdir, "cbmp_up_max_15min_2024.csv")
out_up_mean = joinpath(outdir, "cbmp_up_mean_15min_2024.csv")
out_dn_ext  = joinpath(outdir, "cbmp_down_extreme_15min_2024.csv")
out_dn_mean = joinpath(outdir, "cbmp_down_mean_15min_2024.csv")

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


# --- MAIN ---
function main()
    # Read CSV; empty strings -> missing
    df = DataFrame(CSV.File(infile; missingstring="", normalizenames=true))
    rename!(df, Symbol.(names(df)))  # ensure all column names are Symbols


    col_isp = findcol(df, "ISP_CET_CEST")
    col_ts  = findcol(df, "timestamp")
    col_up  = findcol(df, "Price_Up_EUR_MWh")
    col_dn  = findcol(df, "Price_Down_EUR_MWh")

    # Ensure numeric timestamp
    df[!, col_ts] = Float64.(df[!, col_ts])

    # Bin index (0..35039)
    df[!, :bin] = Int.(floor.(df[!, col_ts] ./ sec_per_bin))
    df = df[(df.bin .>= 0) .& (df.bin .< bins_per_year), :]

    # Identify "no call" rows:
    # rule: no call ONLY if both values are present AND exactly equal
    upv = df[!, col_up]
    dnv = df[!, col_dn]
    no_call = (.!ismissing.(upv)) .& (.!ismissing.(dnv)) .& (upv .== dnv)

    # Build "call" masks for each direction:
    # A row counts as UP-call value if:
    #   - it is NOT a no_call row
    #   - and UP is not missing
    is_up_call = (.!no_call) .& (.!ismissing.(upv))

    # A row counts as DOWN-call value if:
    #   - it is NOT a no_call row
    #   - and DOWN is not missing
    is_dn_call = (.!no_call) .& (.!ismissing.(dnv))

    df_up = df[is_up_call, [:bin, col_up]]
    df_dn = df[is_dn_call, [:bin, col_dn]]

    # Aggregate UP: max and mean
    up_agg = combine(groupby(df_up, :bin),
        col_up => (x -> maximum(skipmissing(x))) => :cbmp_up_max,
        col_up => (x -> mean(collect(skipmissing(x)))) => :cbmp_up_mean
    )

    # Aggregate DOWN:
    # "extreme" (your "highest by magnitude") for DOWN is the MOST NEGATIVE value => minimum
    dn_agg = combine(groupby(df_dn, :bin),
        col_dn => (x -> minimum(skipmissing(x))) => :cbmp_down_extreme,
        col_dn => (x -> mean(collect(skipmissing(x)))) => :cbmp_down_mean
    )

    # Full 15-min timeline (35040 rows)
    bins = collect(0:(bins_per_year-1))
    dt15 = year_dt0 .+ Minute.(15 .* bins)
    ts_str = Dates.format.(dt15, dateformat"dd/mm/yyyy HH:MM")
    base = DataFrame(bin=bins, timestamp=ts_str)

    full = leftjoin(base, up_agg, on=:bin)
    full = leftjoin(full, dn_agg, on=:bin)

    # Write outputs (empty cell for missing = no call in that direction in that 15min)
    CSV.write(out_up_max,  select(full, :timestamp, :cbmp_up_max))
    CSV.write(out_up_mean, select(full, :timestamp, :cbmp_up_mean))
    CSV.write(out_dn_ext,  select(full, :timestamp, :cbmp_down_extreme))
    CSV.write(out_dn_mean, select(full, :timestamp, :cbmp_down_mean))

    println("Done. Files written to: ", outdir)
    println("Rows per output (should be 35040): ", nrow(full))
    println("15-min bins with UP calls: ", count(.!ismissing.(full.cbmp_up_mean)))
    println("15-min bins with DOWN calls: ", count(.!ismissing.(full.cbmp_down_mean)))
    println("4s rows classified as NO CALL: ", count(no_call))
    println("4s rows classified as UP CALL (incl. diff UP/DOWN rows): ", nrow(df_up))
    println("4s rows classified as DOWN CALL (incl. diff UP/DOWN rows): ", nrow(df_dn))
end

main()
