using CSV
using DataFrames
using Dates
using Statistics

println("== Start 15-Minuten-Auswertung (share_f10..share_f100 + callblock_prob + callshare_mean) ==")

# ------------------------------------------------------------
# Einstellungen
# ------------------------------------------------------------
in_file  = raw"C:\Users\jenter\Documents\Backup Masterarbeit\CBMP\aFRR_prices_2024_year_no_leapday.csv"  # <-- anpassen
out_dir  = dirname(in_file)

out_up_15m = joinpath(out_dir, "cbmp_15min_UP_compact.csv")
out_dn_15m = joinpath(out_dir, "cbmp_15min_DOWN_compact.csv")

col_isp = Symbol("ISP (CET/CEST)")
col_ts  = :timestamp
col_up  = Symbol("Price Up (EUR/MWh)")
col_dn  = Symbol("Price Down (EUR/MWh)")

# Faktoren 0.1 .. 1.0
factors = collect(0.1:0.1:1.0)

# 15 Minuten / 4 Sekunden
BLOCK_S = 900
NSTEP   = Int(BLOCK_S ÷ 4)  # 225

println("Lese Input: ", in_file)

# ------------------------------------------------------------
# ISP-Parsing (nur Blockstart, optional für Debug)
# ------------------------------------------------------------
isp_dt_format = dateformat"dd/mm/yyyy HH:MM:SS"
function parse_isp_start(isp::AbstractString)::DateTime
    m = match(r"^\s*(\d{2}/\d{2}/\d{4} \d{2}:\d{2}:\d{2})", isp)
    m === nothing && error("Konnte Startzeit nicht aus ISP extrahieren: $(repr(isp))")
    return DateTime(m.captures[1], isp_dt_format)
end

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
println("Anzahl 15-Min-Blocks (vor Auffüllen): ", length(gdf))

# ------------------------------------------------------------
# Output-Tabellen
# ------------------------------------------------------------
function empty_out()
    DataFrame(
        timestamp = Int[],
        callblock_prob = Float64[],   # 0/1 pro Block
        callshare_mean = Float64[],   # Anteil Abruf-4s-Steps im Block
        share_f10 = Float64[], share_f20 = Float64[], share_f30 = Float64[], share_f40 = Float64[],
        share_f50 = Float64[], share_f60 = Float64[], share_f70 = Float64[], share_f80 = Float64[],
        share_f90 = Float64[], share_f100 = Float64[]
    )
end

out_up = empty_out()
out_dn = empty_out()

# ------------------------------------------------------------
# Pro 15-Minuten-Block berechnen (nur aus vorhandenen Blocks)
# ------------------------------------------------------------
for sub in gdf
    # Blockstart-Timestamp in Sekunden (15-min Raster)
    t0 = Int(sub[1, col_ts] - (sub[1, col_ts] % BLOCK_S))

    up = sub[!, col_up]
    dn = sub[!, col_dn]

    # Abruf-Indikatoren je 4s (wie vorher definiert)
    up_call = .!ismissing.(up) .&  ismissing.(dn)
    dn_call = .!ismissing.(dn) .&  ismissing.(up)

    # ===== UP =====
    cb_up = any(up_call) ? 1.0 : 0.0
    cs_up = count(up_call) / NSTEP

    max_up = any(up_call) ? maximum(skipmissing(up[up_call])) : missing
    s_up = zeros(Float64, length(factors))
    if max_up !== missing
        for (k, f) in pairs(factors)
            bid = f * max_up
            s_up[k] = count(up_call .& (up .>= bid)) / NSTEP
        end
    end

    row_up = (; timestamp = t0, callblock_prob = cb_up, callshare_mean = cs_up,
              share_f10 = s_up[1], share_f20 = s_up[2], share_f30 = s_up[3], share_f40 = s_up[4],
              share_f50 = s_up[5], share_f60 = s_up[6], share_f70 = s_up[7], share_f80 = s_up[8],
              share_f90 = s_up[9], share_f100 = s_up[10])
    push!(out_up, row_up)

    # ===== DOWN =====
    cb_dn = any(dn_call) ? 1.0 : 0.0
    cs_dn = count(dn_call) / NSTEP

    min_dn = any(dn_call) ? minimum(skipmissing(dn[dn_call])) : missing
    s_dn = zeros(Float64, length(factors))
    if min_dn !== missing
        for (k, f) in pairs(factors)
            bid = f * min_dn
            s_dn[k] = count(dn_call .& (dn .<= bid)) / NSTEP
        end
    end

    row_dn = (; timestamp = t0, callblock_prob = cb_dn, callshare_mean = cs_dn,
              share_f10 = s_dn[1], share_f20 = s_dn[2], share_f30 = s_dn[3], share_f40 = s_dn[4],
              share_f50 = s_dn[5], share_f60 = s_dn[6], share_f70 = s_dn[7], share_f80 = s_dn[8],
              share_f90 = s_dn[9], share_f100 = s_dn[10])
    push!(out_dn, row_dn)
end

sort!(out_up, :timestamp)
sort!(out_dn, :timestamp)

# ------------------------------------------------------------
# CHECK + Auffüllen fehlender 15-Minuten-Blöcke mit 0-Abruf
# ------------------------------------------------------------
expected_blocks = Int(365 * 24 * 3600 ÷ BLOCK_S)  # 35040
full_block_ids  = collect(0:(expected_blocks - 1))

function zero_block_row(block_id::Int)
    t = block_id * BLOCK_S
    return (; timestamp = t, callblock_prob = 0.0, callshare_mean = 0.0,
            share_f10 = 0.0, share_f20 = 0.0, share_f30 = 0.0, share_f40 = 0.0,
            share_f50 = 0.0, share_f60 = 0.0, share_f70 = 0.0, share_f80 = 0.0,
            share_f90 = 0.0, share_f100 = 0.0)
end

# UP
existing_blocks_up = Set(div.(out_up.timestamp, BLOCK_S))
missing_blocks_up  = setdiff(full_block_ids, collect(existing_blocks_up))
println("UP – fehlende Blöcke: ", length(missing_blocks_up))
if !isempty(missing_blocks_up)
    println("UP – fehlende block_id(s): ", missing_blocks_up)
    for b in missing_blocks_up
        push!(out_up, zero_block_row(b))
    end
    sort!(out_up, :timestamp)
end

# DOWN
existing_blocks_dn = Set(div.(out_dn.timestamp, BLOCK_S))
missing_blocks_dn  = setdiff(full_block_ids, collect(existing_blocks_dn))
println("DOWN – fehlende Blöcke: ", length(missing_blocks_dn))
if !isempty(missing_blocks_dn)
    println("DOWN – fehlende block_id(s): ", missing_blocks_dn)
    for b in missing_blocks_dn
        push!(out_dn, zero_block_row(b))
    end
    sort!(out_dn, :timestamp)
end

# Sanity checks
@assert nrow(out_up) == expected_blocks "UP hat nicht $(expected_blocks) Zeilen!"
@assert nrow(out_dn) == expected_blocks "DOWN hat nicht $(expected_blocks) Zeilen!"
@assert out_up.timestamp[1] == 0
@assert out_up.timestamp[end] == (expected_blocks - 1) * BLOCK_S

println("Anzahl 15-Min-Blocks (nach Auffüllen): ", expected_blocks)

# ------------------------------------------------------------
# Export
# ------------------------------------------------------------
CSV.write(out_up_15m, out_up; delim=',', quotechar='"')
CSV.write(out_dn_15m, out_dn; delim=',', quotechar='"')

println("== Fertig ==")
println("UP   -> ", out_up_15m, " | rows: ", nrow(out_up))
println("DOWN -> ", out_dn_15m, " | rows: ", nrow(out_dn))
println("\nHinweis:")
println("- callblock_prob ist pro 15-Min-Block 0 oder 1 (Mittelung => Wahrscheinlichkeit).")
println("- callshare_mean ist pro 15-Min-Block der Zeitanteil (0..1) an 4s-Steps mit Abruf in der Richtung.")
