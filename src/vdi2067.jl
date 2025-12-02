"""
    vdi2067_annuity(sim, components; params=VDIParams())

Berechne nach VDI-Regelwerk (vereinfachte Version)
- Kapitalgebundene jährliche Annuität (inkl. Ersatzbeschaffungen & Restwert)
- Betriebs-, bedarfsabhängige Kosten, Sonstige Kosten
- Erlöse (Einspeisung, Regelreservevermarktung)
- Gesamtannuität und Wärmepreis (EUR/kWh)

Inputs:
- sim::Dict: Simulations-Dictionary mit Zeitreihen (Wh) und Preisprofile (EUR/MWh)
  Erwartete Keys (mindestens):
    "Grid_IN" :: Vector{Float64}  # Netzbezug (Wh)
    "Grid_Out" :: Vector{Float64} # Einspeisung (Wh)
    "Heat"    :: Vector{Float64}  # Wärmebedarf (Wh)
    optional: "PV_selfconsumption" :: Vector{Float64}
    "Grid_price" either scalar or Vector{Float64} (EUR/MWh)
    "FeedIn_price" either scalar or Vector{Float64} (EUR/MWh)
    "Regelenergie_price" Vector{Float64} (EUR/MWh)
    "Regelleistung_price" Vector{Float64} (EUR/MW)  # Bei Nutzung mit Vorhalteleistung in MW multiplizieren
    "Reserved_power" scalar or Vector (W) (Vorgehaltene Leistung in W)
- components::Vector{Component}:
    Strukturen zur Beschreibung der Investkomponenten (A0, Lebensdauer, Name)
- params::VDIParams: Zins- / Preisänderungsfaktoren und Default-Faktoren

Outputs:
- Dict mit Annuitäten-Komponenten und Wärmepreis (€/kWh)
"""
module VDI2067

export VDIParams, Component, vdi2067_annuity

struct VDIParams
    T::Int
    i_cap::Float64
    r_cap::Float64
    i_op::Float64
    r_op::Float64
    i_energy::Float64
    r_energy::Float64
    i_revenue::Float64
    r_revenue::Float64
    f_inst::Float64
    f_winsp::Float64
    spec_inspection_hours::Float64
    grid_price_addon::Float64
    include_replacements_at_horizon::Bool  # Steuerung: Ersatzbeschaffungen, die genau auf T liegen, berücksichtigen
end

VDIParams(; T=20,
          i_cap=0.03, r_cap=0.00,
          i_op=0.03, r_op=0.01,
          i_energy=0.03, r_energy=0.02,
          i_revenue=0.03, r_revenue=0.00,
          f_inst=0.01, f_winsp=0.015, spec_inspection_hours=5.0,
          grid_price_addon=0.0,
          include_replacements_at_horizon=true) = VDIParams(T, i_cap, r_cap, i_op, r_op, i_energy, r_energy, i_revenue, r_revenue, f_inst, f_winsp, spec_inspection_hours, grid_price_addon, include_replacements_at_horizon)

# Einfache Komponentendefinition für kapitalbezogene Kalkulationen
struct Component
    name::String
    A0::Float64    # Investkosten zum Zeitpunkt 0 in EUR
    TN::Int        # Nutzungsdauer in Jahren
    share_peripheral::Float64 # optionaler Anteil, falls nur Teil der Anlage zugehörig
end

"""
Hilfsfunktionen
""" 
# Annuitätsfaktor (VDI-Form): a = q^T * (q - 1) / (q^T - 1), mit q = 1 + i. Grenzfall i == 0 -> 1/T
function annuity_factor(i::Float64, T::Int)
    if T == 0
        return 0.0
    end
    q = 1.0 + i
    # Grenzfall q == 1 => i == 0 -> Limit a = 1 / T
    if isapprox(q, 1.0; atol=1e-12)
        return 1.0 / T
    end
    qT = q^T
    return qT * (q - 1.0) / (qT - 1.0)
end

# Preisänderungsfaktor (transformiert den Betrag des 1. Jahres in Barwert über T),
# gemäß VDI: b = (1 - (q_v / q)^T) / (q - q_v), wobei q = 1 + i, q_v = 1 + r
function price_change_factor(r, i, T)
    """
    Berechnet den VDI-Preisänderungsfaktor:
      b = (1 - (q_v/q)^T) / (q - q_v)
    mit q = 1 + i, q_v = 1 + r.
    Grenzfall q == q_v -> b = T / q.
    """
    if T == 0
        return 1.0
    end
    q = 1.0 + i
    qv = 1.0 + r
    denom = q - qv
    if isapprox(denom, 0.0; atol=1e-12)
        return T / q
    end
    return (1.0 - (qv / q)^T) / denom
end

# Wandle einen skalaren Preis in einen Vektor um, falls erforderlich
function vecize_price(price, N)
    if isa(price, Number)
        return fill(price, N)
    elseif length(price) == N
        return price
    else
        error("Preis-Array-Längenfehler: erwartet Länge $N, erhalten $(length(price))")
    end
end

# Barwert geplanter Ersatzbeschaffungen (ohne die anfängliche Investition A0)
function pv_replacements(A0::Float64, TN::Int, T::Int, i::Float64, r::Float64)
    pv = 0.0
    if TN <= 0
        return pv
    end
    t = TN
    while t < T - 1e-9 # Ersatztermine bei TN, 2*TN, ... strikt kleiner als T
        price_at_purchase = A0 * (1 + r)^(t)
        pv += price_at_purchase / (1 + i)^t
        t += TN
    end
    return pv
end

# Barwert des Restwerts zum Zeitpunkt T und Abzinsung auf t=0 (VDI-konform)
function pv_restwert(A0::Float64, TN::Int, T::Int, i::Float64, r::Float64; include_repl_at_T::Bool=true)
    # VDI-konform:
    # - Ermittle die zuletzt getätigte Beschaffung (vor oder auf T)
    # - Falls Ersatz genau auf T erfolgt, ergibt das den vollen Restwert = A0*(1+r)^T
    # - andernfalls linear anteiliger Restwert nach nicht-genutzter Restnutzungsdauer
    if TN <= 0 || T <= 0
        return 0.0
    end

    # Behandlung eines Ersatzes, der exakt auf T fällt:
    n = floor(Int, T / TN)                # Anzahl vollständiger Lebenszyklen
    last_purchase = n * TN

    if !include_repl_at_T && last_purchase == T && n > 0
        # Option: Ersatz am T als nach dem Betrachtungszeitraum interpretieren
        last_purchase = (n - 1) * TN
    end

    # last_purchase <= T (ggf. gleich T)
    if last_purchase == T
        # Ersatz fällt genau auf T -> voller Restwert (Preisfortschreibung bis T)
        rest_at_T = A0 * (1.0 + r)^T
    else
        # letzte Beschaffung liegt vor T -> linearer Restwert
        age_at_T = T - last_purchase
        remaining_fraction = max(0.0, (TN - age_at_T) / TN)
        price_at_last_purchase = A0 * (1.0 + r)^(last_purchase)
        rest_at_T = price_at_last_purchase * remaining_fraction
    end

    pv = rest_at_T / (1.0 + i)^T
    return pv
end

# Barwert aller kapitalgebundenen Kosten (Anfangsinvestition + Ersatzbeschaffung - Restwert)
function pv_total_capital(A0::Float64, TN::Int, T::Int, i::Float64, r::Float64)
    pv_init = A0  # Zeitpunkt 0
    pv_repl = pv_replacements(A0, TN, T, i, r)
    pv_rest = pv_restwert(A0, TN, T, i, r)
    return pv_init + pv_repl - pv_rest
end

# Berechne die jährliche Kapitalannuität für eine Komponente
function annuity_capital_component(comp::Component, params::VDIParams)
    pv = pv_total_capital(comp.A0, comp.TN, params.T, params.i_cap, params.r_cap)
    a = annuity_factor(params.i_cap, params.T)
    return pv * a
end

# Energiekosten im ersten Jahr aus den Zeitreihen (EUR)
function first_year_energy_cost(sim::Dict, params::VDIParams)
    # Erwartet: alle Energiedaten in Wh, Preise in EUR/MWh
    Grid_IN = sim["Grid_IN"]
    N = length(Grid_IN)
    grid_price = vecize_price(sim["Grid_price"], N)
    if haskey(sim, "PV_selfconsumption")
        pv_self = sim["PV_selfconsumption"]
    else
        pv_self = zeros(N)
    end
    # Aufschlag für Steuern/Abgaben (absoluter Aufschlag in EUR/MWh)
    grid_price = [p + params.grid_price_addon for p in grid_price]
    # Für Eigenerzeugung falls ein Preisprofil angegeben ist (z. B. Grenzkosten)
    feedin_price = haskey(sim, "FeedIn_price") ? vecize_price(sim["FeedIn_price"], N) : zeros(N)
    selfgen_price = haskey(sim, "SelfGen_price") ? vecize_price(sim["SelfGen_price"], N) : fill(0.0, N)

    # Kosten für Netzbezug:
    cost_grid = sum( Grid_IN .* grid_price ) * 1e-6   # EUR
    # Kosten für selbstverbrauchte Eigenerzeugung (falls Grenzkosten vorhanden)
    cost_self = sum( pv_self .* selfgen_price ) * 1e-6
    return cost_grid + cost_self
end

# Annuitiere Energiekosten
function annuity_energy(sim::Dict, params::VDIParams)
    E1 = first_year_energy_cost(sim, params)
    b = price_change_factor(params.r_energy, params.i_energy, params.T)
    a = annuity_factor(params.i_energy, params.T)
    return E1 * b * a
end

# Erlöse des ersten Jahres
function first_year_revenues(sim::Dict)
    N = length(sim["Grid_IN"])
    # Regelenergie-Erträge: Produkt aus Netzbezug und Regelenergie-Preisprofil
    re_price = haskey(sim, "Regelenergie_price") ? vecize_price(sim["Regelenergie_price"], N) : zeros(N)
    revenue_re_energie = sum(sim["Grid_IN"] .* re_price) * 1e-6
    # Regelleistungs-Erträge: Preis * vorgehaltene Leistung (Konvertierung W -> MW)
    regelleistung_price = haskey(sim, "Regelleistung_price") ? vecize_price(sim["Regelleistung_price"], N) : zeros(N)
    reserved_power_W = haskey(sim, "Reserved_power") ? sim["Reserved_power"] : 0.0
    reserved_power_MW = reserved_power_W * 1e-6
    revenue_regelleistung = sum(regelleistung_price .* reserved_power_MW) # Preis in EUR/MW -> EUR
    # Einspeisevergütung Eigenerzeugung:
    grid_out = haskey(sim, "Grid_Out") ? sim["Grid_Out"] : zeros(N)
    feed_in_price = haskey(sim, "FeedIn_price") ? vecize_price(sim["FeedIn_price"], N) : zeros(N)
    revenue_feedin = sum(grid_out .* feed_in_price) * 1e-6

    # Summe der Erlöse
    return revenue_re_energie + revenue_regelleistung + revenue_feedin
end

function annuity_revenues(sim::Dict, params::VDIParams)
    E1 = first_year_revenues(sim)
    b = price_change_factor(params.r_revenue, params.i_revenue, params.T)
    a = annuity_factor(params.i_revenue, params.T)
    return E1 * b * a
end

# Betriebskosten (Wartung etc.) basieren auf Anteil der Investkosten oder Stundenaufwand * Stundensatz
# Optional können 'op_bases' im Simulations-Dict bereitgestellt werden; Standard ist Prozent der Investition
function first_year_operation_costs(components::Vector{Component}, params::VDIParams)
    # Standard: prozentualer Anteil an den Investitionskosten
    total_cap = sum(c.A0 for c in components)
    # Verwende die Parameter-Faktoren (f_inst & f_winsp)
    inst_cost = total_cap * params.f_inst
    winsp_cost = total_cap * params.f_winsp
    # zusätzliche Basis (kann im sim-Dict angegeben werden)
    return inst_cost + winsp_cost
end

function annuity_operation(components::Vector{Component}, params::VDIParams)
    # Basis des ersten Jahres
    A1 = first_year_operation_costs(components, params)
    a = annuity_factor(params.i_op, params.T)
    b = price_change_factor(params.r_op, params.i_op, params.T)
    return A1 * a * b
end

# Sonstige Kosten: Prozentsatz auf Investitionen (Versicherungen, Verwaltung, o.ä.)
function annuity_sonstige(components::Vector{Component}, params::VDIParams; pct=0.02)
    total_cap = sum(c.A0 for c in components)
    A1 = total_cap * pct
    a = annuity_factor(params.i_op, params.T)
    b = price_change_factor(params.r_op, params.i_op, params.T)
    return A1 * a * b
end

# Hauptfunktion zur Berechnung der Annuitäten nach VDI
function vdi2067_annuity(sim::Dict, components::Vector{Component}; params=VDIParams())
    # Kapitalannuitäten für alle Komponenten berechnen
    capital_annuity_per_comp = Dict{String, Float64}()
    for c in components
        capital_annuity_per_comp[c.name] = annuity_capital_component(c, params)
    end
    A_cap_total = sum(values(capital_annuity_per_comp))

    # Betriebskosten
    A_betrieb = annuity_operation(components, params)

    # Bedarfsgebundene Energiekosten
    A_bedarf = annuity_energy(sim, params)

    # Sonstige Kosten
    A_sonst = annuity_sonstige(components, params)

    # Erlöse
    A_erloese = annuity_revenues(sim, params)

    # Gesamtannuität
    A_ges = A_cap_total + A_betrieb + A_bedarf + A_sonst - A_erloese

    # Wärmepreis: Anteil pro kWh Wärme
    heat_profile = sim["Heat"]
    annual_heat_kWh = sum(heat_profile) / 1000.0
    heat_price_eur_per_kwh = annual_heat_kWh == 0 ? Inf : A_ges / annual_heat_kWh

    return Dict(
        "A_cap_total" => A_cap_total,
        "A_betrieb" => A_betrieb,
        "A_bedarf" => A_bedarf,
        "A_sonst" => A_sonst,
        "A_erloese" => A_erloese,
        "A_ges" => A_ges,
        "heat_price_eur_per_kwh" => heat_price_eur_per_kwh,
        "annual_heat_kwh" => annual_heat_kWh,
        "capital_annuity_breakdown" => capital_annuity_per_comp
    )
end

end # module VDI2067