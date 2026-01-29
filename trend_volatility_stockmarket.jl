using CSV                      # CSV-Dateien lesen und schreiben
using DataFrames               # Tabellenverarbeitung
using Statistics               # mean, std

# ------------------------------------------------------------
# Hilfsfunktion: robuste Konvertierung nach Float64
# - akzeptiert Number, String, Missing
# - ersetzt Dezimalkomma
# - leere Strings → missing
# ------------------------------------------------------------
function tofloat(x)
    if x isa Missing
        return missing
    elseif x isa Number
        return Float64(x)
    else
        s = strip(String(x))                 # Whitespace entfernen
        s = replace(s, "," => ".")           # Dezimalkomma → Punkt
        return s == "" ? missing : tryparse(Float64, s)
    end
end

# ------------------------------------------------------------
# Kernfunktion: Trend & Volatilität berechnen
# ------------------------------------------------------------
function compute_trend_volatility(
    t::Vector{Float64},                     # timestep in Sekunden
    p::Vector{Float64};                     # Preis in €/MWh
    horizon_s::Int = 3 * 3600                # Vorwärtshorizont (3 h)
)

    n = length(t)                            # Anzahl Zeitschritte

    trend = Vector{Union{Missing, Float64}}(missing, n)  # €/MWh pro Stunde
    vol   = Vector{Union{Missing, Float64}}(missing, n)  # €/MWh

    j_end = 1                                # Fenster-Endindex (Two-Pointer)

    for i in 1:n                             # Über alle Startzeitpunkte
        while j_end < n && t[j_end + 1] <= t[i] + horizon_s
            j_end += 1                       # Fenster bis +3h erweitern
        end

        if j_end - i + 1 < 2                 # Mindestens 2 Punkte nötig
            continue
        end

        idx = i:j_end                        # Indexbereich des Fensters
        x = (t[idx] .- t[i]) ./ 3600.0       # Zeit relativ zu i in Stunden
        y = p[idx]                           # Preise im Fenster

        vol[i] = std(y)                      # Volatilität = Std-Abw.

        x̄ = mean(x)                          # Mittelwert Zeit
        ȳ = mean(y)                          # Mittelwert Preis

        denom = sum((x .- x̄).^2)             # Nenner Regression

        if denom != 0.0
            trend[i] = sum((x .- x̄) .* (y .- ȳ)) / denom
            # → Steigung in €/MWh pro Stunde
        end
    end

    return trend, vol
end

# ------------------------------------------------------------
# Einlesen der CSV-Datei
# ------------------------------------------------------------
input_path  = "C:/Users/jenter/Documents/resie/profiles/MA/csv/DE_Grosshandelspreise_2024_Viertelstunde.csv"
output_path = "C:/Users/jenter/Documents/resie/profiles/MA/csv/prices_with_trend_vol_3h.csv"

df_raw = CSV.read(input_path, DataFrame; normalizenames=true)

# Sicherstellen, dass die benötigten Spalten existieren
@assert :timestep in propertynames(df_raw)
@assert :price_power_EUR_MWh in propertynames(df_raw)

# Robuste Typkonvertierung
df = DataFrame(
    timestep            = [tofloat(v) for v in df_raw.timestep],
    price_power_EUR_MWh = [tofloat(v) for v in df_raw.price_power_EUR_MWh]
)

# Zeilen mit fehlenden Werten entfernen
df = dropmissing(df)

# In Float64-Vektoren überführen
t = Vector{Float64}(df.timestep)
p = Vector{Float64}(df.price_power_EUR_MWh)

# ------------------------------------------------------------
# Trend & Volatilität berechnen
# ------------------------------------------------------------
trend, volatility = compute_trend_volatility(t, p)

# ------------------------------------------------------------
# Output-DataFrame erstellen
# ------------------------------------------------------------
df_out = DataFrame(
    timestep              = df.timestep,
    price_power_EUR_MWh   = df.price_power_EUR_MWh,
    trend_EUR_MWh_per_h   = trend,
    volatility_EUR_MWh    = volatility
)

# ------------------------------------------------------------
# CSV schreiben
# ------------------------------------------------------------
CSV.write(output_path, df_out)

println("✔ Fertig. Datei geschrieben: ", output_path)
