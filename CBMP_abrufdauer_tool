using CSV
using DataFrames
using Dates
using Statistics

println("== Start Monatsauswertung (nur share_f10..share_f100 + callblock_prob + callshare_mean) ==")

# ------------------------------------------------------------
# Einstellungen
# ------------------------------------------------------------
in_file  = raw"C:/Users/jenter/Documents/Backup Masterarbeit/CBMP/aFRR_prices_2024_year_no_leapday.csv"  # <-- anpassen
out_dir  = dirname(in_file)

out_up_month = joinpath(out_dir, "cbmp_month_UP_compact.csv")
out_dn_month = joinpath(out_dir, "cbmp_month_DOWN_compact.csv")

col_isp = Symbol("ISP (CET/CEST)")
col_ts  = :timestamp
col_up  = Symbol("Price Up (EUR/MWh)")
col_dn  = Symbol("Price Down (EUR/MWh)")

# Faktoren 0.1 .. 1.0
factors = collect(0.1:0.1:1.0)

# 15 Minuten / 4 Sekunden
BLOCK_S = 900
NSTEP   = Int(BLOCK_S ÷ 4)  # 225

# ISP-Parsing (nur Blockstart)
isp_dt_format = dateformat"dd/mm/yyyy HH:MM:SS"
function parse_isp_start(isp::AbstractString)::DateTime
    m = match(r"^\s*(\d{2}/\d{2}/\d{4} \d{2}:\d{2}:\d{2})", isp)
    m === nothing && error("Konnte Startzeit nicht aus ISP extrahieren: $(repr(isp))")
    return DateTime(m.captures[1], isp_dt_format)
end

println("Lese Input: ", in_file)

# ------------------------------------------------------------
# Einlesen
# ------------------------------------------------------------
df = CSV.read(in_file, DataFrame;
    delim = ',',
    quotechar = '"',
    pool = false,
    missingstring = ["", "NA", "NaN", "N/A"],
    types = Dict(
        col_isp => Union{Missing,String},
        col_ts  => Int,
        col_up  => Union{Missing,Float64},
        col_dn  => Union{Missing,Float64},
    ),
)

dropmissing!(df, col_isp)
sort!(df, col_ts)

# Block-ID aus timestamp (ab 0)
df.block_id = div.(df[!, col_ts], BLOCK_S)

gdf = groupby(df, :block_id)
println("Anzahl 15-Min-Blocks: ", length(gdf))

# ------------------------------------------------------------
# Pro-Block Sammeln: UP und DOWN separat
# ------------------------------------------------------------
months = Int[]
callshare_up = Float64[]
callblock_up = Int[]
shares_up = Vector{Vector{Float64}}()

callshare_dn = Float64[]
callblock_dn = Int[]
shares_dn = Vector{Vector{Float64}}()

for sub in gdf
    dt0 = parse_isp_start(sub[1, col_isp])
    mth = month(dt0)
    push!(months, mth)

    up = sub[!, col_up]
    dn = sub[!, col_dn]

    up_call = .!ismissing.(up) .&  ismissing.(dn)
    dn_call = .!ismissing.(dn) .&  ismissing.(up)

    # preisunabhängig
    push!(callshare_up, count(up_call) / NSTEP)
    push!(callblock_up, any(up_call) ? 1 : 0)

    push!(callshare_dn, count(dn_call) / NSTEP)
    push!(callblock_dn, any(dn_call) ? 1 : 0)

    # Grenzpreise über Abruf-4s-Steps
    max_up = any(up_call) ? maximum(skipmissing(up[up_call])) : missing
    min_dn = any(dn_call) ? minimum(skipmissing(dn[dn_call])) : missing

    s_up = zeros(Float64, length(factors))
    s_dn = zeros(Float64, length(factors))

    for (k, f) in pairs(factors)
        # UP: Grenzpreis = max, Abruf wenn CBMP(t) >= bid (= f*max)
        if max_up !== missing
            bid = f * max_up
            s_up[k] = count(up_call .& (up .>= bid)) / NSTEP
        end

        # DOWN: Grenzpreis = min, Abruf wenn CBMP(t) <= bid (= f*min)
        if min_dn !== missing
            bid = f * min_dn
            s_dn[k] = count(dn_call .& (dn .<= bid)) / NSTEP
        end
    end

    push!(shares_up, s_up)
    push!(shares_dn, s_dn)
end

# ------------------------------------------------------------
# Monatsmittel bilden
# ------------------------------------------------------------
function month_mean_compact(months::Vector{Int}, callshare::Vector{Float64}, callblock::Vector{Int},
                            shares::Vector{Vector{Float64}}, factors::Vector{Float64})
    dfB = DataFrame(month = months, callshare_mean = callshare, callblock = callblock)

    # share_f10..share_f100 als Spalten
    for (k, f) in pairs(factors)
        dfB[!, Symbol("share_f$(Int(round(f*100)))")] = [shares[i][k] for i in eachindex(shares)]
    end

    share_cols = [Symbol("share_f$(k)") for k in 10:10:100]

    g = groupby(dfB, :month)
    out = combine(g) do sdf
        row = DataFrame(month = [sdf[1, :month]])
        row.callshare_mean = [mean(sdf.callshare_mean)]
        row.callblock_prob = [mean(sdf.callblock)]
        for c in share_cols
            row[!, c] = [mean(sdf[!, c])]
        end
        row
    end

    sort!(out, :month)
    return out
end

up_month = month_mean_compact(months, callshare_up, callblock_up, shares_up, factors)
dn_month = month_mean_compact(months, callshare_dn, callblock_dn, shares_dn, factors)

# Output exakt gewünschte Spaltenreihenfolge
share_cols = [Symbol("share_f$(k)") for k in 10:10:100]
select!(up_month, vcat([:month, :callblock_prob, :callshare_mean], share_cols))
select!(dn_month, vcat([:month, :callblock_prob, :callshare_mean], share_cols))

# ------------------------------------------------------------
# Export
# ------------------------------------------------------------
CSV.write(out_up_month, up_month; delim=',', quotechar='"')
CSV.write(out_dn_month, dn_month; delim=',', quotechar='"')

println("== Fertig ==")
println("UP   -> ", out_up_month)
println("DOWN -> ", out_dn_month)
