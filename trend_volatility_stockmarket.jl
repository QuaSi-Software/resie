using CSV                          # CSV-Dateien lesen/schreiben
using DataFrames                   # Tabellarische Daten verarbeiten
using Statistics                   # Mittelwert/Standardabweichung
using Printf                       # (optional) saubere Formatierung

input_path  = "C:/Users/jenter/Documents/resie/profiles/MA/csv/DE_Grosshandelspreise_2024_Viertelstunde.csv"   # Pfad zur Eingabe-CSV
output_path = "C:/Users/jenter/Documents/resie/profiles/MA/csv"                  # Pfad zur Ausgabe-CSV

horizon_s = 3 * 60 * 60           # Trend-/Volatilitäts-Horizont: 3 Stunden in Sekunden (10800 s)

df = CSV.read(input_path, DataFrame)                                    # CSV einlesen (erwartet Spalten: timestep, price_power_EUR_MWh)

t = Vector{Float64}(df.timestep)                                        # timestep als Float64-Vektor (Sekunden)
p = Vector{Float64}(df.price_power_EUR_MWh)                             # Preis als Float64-Vektor (€/MWh)

n = nrow(df)                                                            # Anzahl Zeilen / Zeitpunkte

trend_slope_EUR_MWh_per_h = Vector{Union{Missing, Float64}}(missing, n) # Output-Vektor: Trend (Steigung) in €/MWh pro Stunde
volatility_EUR_MWh        = Vector{Union{Missing, Float64}}(missing, n) # Output-Vektor: Volatilität (Std-Abw.) in €/MWh

j_end = 1                                                               # Laufender Index für das Fensterende (Two-Pointer-Technik, effizient)

for i in 1:n                                                            # Über alle Zeitpunkte laufen
    while j_end < n && t[j_end + 1] <= t[i] + horizon_s                 # Fensterende so weit schieben, bis "nächste 3 Stunden" abgedeckt sind
        j_end += 1                                                      # Fensterende inkrementieren
    end                                                                  # Ende der while-Schleife

    if j_end - i + 1 < 2                                                # Wenn weniger als 2 Punkte im Fenster sind (Regression nicht möglich)
        trend_slope_EUR_MWh_per_h[i] = missing                          # Trend nicht berechenbar
        volatility_EUR_MWh[i]        = missing                          # Volatilität nicht berechenbar
        continue                                                        # Nächster Zeitschritt
    end                                                                  # Ende der Mindestanzahl-Prüfung

    idx = i:j_end                                                       # Indexbereich für das Vorwärtsfenster (inkl. aktueller Punkt bis +3h)
    x = (t[idx] .- t[i]) ./ 3600.0                                      # Zeitachse relativ zum aktuellen Zeitpunkt i, in Stunden (0..3h)
    y = p[idx]                                                          # Preise im Fenster

    volatility_EUR_MWh[i] = std(y)                                      # Volatilität als Standardabweichung der Preise im Fenster

    x_mean = mean(x)                                                    # Mittelwert der Zeitwerte
    y_mean = mean(y)                                                    # Mittelwert der Preiswerte

    denom = sum((x .- x_mean) .^ 2)                                     # Nenner der Steigung (Varianz von x * (n-1) in Summenform)
    if denom == 0.0                                                     # Falls alle x gleich wären (sollte praktisch nur bei Fehlerdaten passieren)
        trend_slope_EUR_MWh_per_h[i] = missing                          # Trend nicht definierbar
    else                                                                 # Normalfall: Regression möglich
        numer = sum((x .- x_mean) .* (y .- y_mean))                     # Zähler der Steigung (Kovarianz-Summe)
        trend_slope_EUR_MWh_per_h[i] = numer / denom                    # Steigung a der linearen Regression y = a*x + b (€/MWh pro Stunde)
    end                                                                  # Ende denom-Prüfung
end                                                                      # Ende for-Schleife

df_out = DataFrame(                                                     # Neue Ausgabetabelle erstellen
    timestep              = df.timestep,                                # timestep unverändert übernehmen (Sekunden)
    price_power_EUR_MWh   = df.price_power_EUR_MWh,                     # Preis unverändert übernehmen (€/MWh)
    trend_EUR_MWh_per_h   = trend_slope_EUR_MWh_per_h,                  # Trend (Steigung) über die nächsten 3h
    volatility_EUR_MWh    = volatility_EUR_MWh                          # Volatilität (Std-Abw.) über die nächsten 3h
)

CSV.write(output_path, df_out)                                          # Ausgabe-CSV schreiben

println("Fertig. Datei geschrieben: ", output_path)                     # Kurzes Status-Print
