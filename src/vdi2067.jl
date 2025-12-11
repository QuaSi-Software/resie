module VDI2067
using OrderedCollections: OrderedDict

export VDIParams, VDIComponent, vdi2067_annuity,
       VDI_SCENARIO_NO, VDI_SCENARIO_MOD, VDI_SCENARIO_PRO,
       heatpump_component, boiler_component, buffertank_component, battery_component

############################################################
#  PARAMETER STRUCTURE – VDI 2067 COMPATIBLE
############################################################

struct VDIParams
    T::Int            # Evaluation period (years)
    i_cap::Float64    # Discount rate for capital-related annuity
    r_cap::Float64    # Escalation of capital costs (investments, replacements)
    r_energy::Float64 # Escalation of energy prices (demand-dependent costs)
    r_op::Float64     # Escalation of maintenance & inspection
    r_inst::Float64   # Escalation of repair / reinstatement
    r_rev::Float64    # Escalation of revenues (feed-in, control energy)
    grid_price_addon_eur_mwh::Float64   # taxes, levies, grid fees
end


############################################################
#  COMPONENT DEFINITION — COST GROUP 420/430
############################################################
# A0        : initial investment at t = 0 (EUR)
# TN        : technical lifetime (years)
# f_wartung : maintenance + inspection factor (share of A0 in year 1)
# f_instand : repair / reinstatement factor (share of A0 in year 1)
############################################################

struct VDIComponent
    name::String
    A0::Float64
    TN::Int
    f_wartung::Float64
    f_instand::Float64
end


############################################################
#  CONVENIENCE CONSTRUCTORS FOR YOUR FOUR COMPONENT TYPES
############################################################
# Only A0 and TN must be provided by the calling program.
# The maintenance/repair factors are hard-coded according
# to your specification / VDI assumptions.
############################################################

"Heat pump: maintenance 1.5 %, repair 1.0 % of A0 in first year."
heatpump_component(A0::Real, TN::Int) =
    VDIComponent("HeatPump", float(A0), TN, 0.015, 0.01)

"Boiler / electric heater: maintenance 2.0 %, repair 1.0 %."
boiler_component(A0::Real, TN::Int) =
    VDIComponent("Boiler", float(A0), TN, 0.02, 0.01)

"Buffer tank: maintenance 0.5 %, repair 0.5 %."
buffertank_component(A0::Real, TN::Int) =
    VDIComponent("BufferTank", float(A0), TN, 0.005, 0.005)

"Battery: maintenance 1.0 %, repair 2.0 %."
battery_component(A0::Real, TN::Int) =
    VDIComponent("Battery", float(A0), TN, 0.01, 0.02)


############################################################
#  SCENARIO DEFINITIONS (NO / MODERATE / PROGRESSIVE)
############################################################

# Scenario 1 — no escalation
VDI_SCENARIO_NO = VDIParams(
    20,       # Evaluation period (years)
    0.03,     # discount rate +3.0%
    0.000,    # no capital cost 
    0.000,    # no energy cost 
    0.000,    # no maintenance cost 
    0.000,    # no repair cost 
    0.000,    # no revenues 
    120.0     # grid price addon €/MWh
)

# Scenario 2 — moderate escalation
VDI_SCENARIO_MOD = VDIParams(
    20,       # Evaluation period (years)
    0.03,     # discount rate +3.0%
    0.012,    # capital cost +1.2%
    0.020,    # energy cost +1.0%
    0.005,    # maintenance cost +0.5%
    0.005,    # repair cost +0.5%
    -0.01,    # revenues -10%
    120.0     # grid price addon €/MWh
)

# Scenario 3 — progressive escalation
VDI_SCENARIO_PRO = VDIParams(
    20,       # Evaluation period (years)
    0.03,     # discount rate +3.0%
    0.018,    # capital +1.8%
    0.040,    # energy cost +2.0%
    0.010,    # maintenance +1.0%
    0.010,    # repair +1.0%
    -0.02,    # revenues -2.0%
    120.0     # grid price addon €/MWh
)


############################################################
#  VDI 2067 MATHEMATICAL FUNCTIONS
############################################################

# Annuity factor a = q^T * (q − 1) / (q^T − 1), q = 1 + i
function annuity_factor(i, T)
    q = 1 + i
    T == 0 && return 0.0
    isapprox(q, 1.0) && return 1/T
    return (q^T * (q - 1)) / (q^T - 1)
end

# Price change factor b = (1 − (qv / q)^T) / (q − qv)
function price_change_factor(r, i, T)
    q = 1 + i
    qv = 1 + r
    isapprox(q, qv) && return T / q
    return (1 - (qv/q)^T) / (q - qv)
end

# Accept scalar or vector price input
vecize_price(p, N) =
    isa(p, Number) ? fill(p, N) :
    (length(p) == N ? p : error("Price vector length mismatch."))


############################################################
#  CAPITAL COSTS — REPLACEMENTS AND RESIDUAL VALUE
############################################################

function npv_replacements(A0, TN, T, i, r)
    npv = 0.0
    t = TN
    while t < T
        cost_t = A0 * (1+r)^t
        npv += cost_t / (1+i)^t
        t += TN
    end
    return npv
end

function npv_restwert(A0, TN, T, i, r)
    n = div(T, TN)
    last_purchase = n * TN

    if last_purchase == T
        rest_T = A0 * (1+r)^T
    else
        age = T - last_purchase
        frac = max(0.0, (TN - age) / TN)
        rest_T = A0 * (1+r)^last_purchase * frac
    end

    return rest_T / (1+i)^T
end

function capital_annuity(comp::VDIComponent, p::VDIParams)
    pv_total = comp.A0 +
        npv_replacements(comp.A0, comp.TN, p.T, p.i_cap, p.r_cap) -
        npv_restwert(comp.A0, comp.TN, p.T, p.i_cap, p.r_cap)

    return pv_total * annuity_factor(p.i_cap, p.T)
end


############################################################
#  OPERATING COSTS: MAINTENANCE (B) + REPAIR (IN)
############################################################

function op_annuity(components::Vector{VDIComponent}, p::VDIParams)
    # First-year maintenance & inspection cost
    A_B1  = sum(c.A0 * c.f_wartung for c in components)
    # First-year repair / reinstatement cost
    A_IN1 = sum(c.A0 * c.f_instand for c in components)

    a   = annuity_factor(p.i_cap, p.T)
    b_B  = price_change_factor(p.r_op,   p.i_cap, p.T)
    b_IN = price_change_factor(p.r_inst, p.i_cap, p.T)

    return A_B1 * a * b_B + A_IN1 * a * b_IN
end


############################################################
#  MISCELLANEOUS COSTS (INSURANCE, ADMIN,…)
############################################################

function misc_annuity(components::Vector{VDIComponent}, p::VDIParams)
    total_cap = sum(c.A0 for c in components)
    A1 = 0.02 * total_cap        # rule-of-thumb: 2% of investment

    a = annuity_factor(p.i_cap, p.T)
    b = price_change_factor(p.r_op, p.i_cap, p.T)

    return A1 * a * b
end


############################################################
#  ENERGY COSTS — TIME-RESOLVED → ANNUITY
############################################################

function energy_annuity(sim::Union{Dict, OrderedDict}, p::VDIParams)
    IN = sim["Grid_IN"]                                         # Wh time series
    base_price = vecize_price(sim["Grid_price"], length(IN))    # €/MWh (market price)

    # Add fixed components (taxes, grid fees, levies)
    # This is a linear addition: total electricity price = market + fixed components
    price_eur_mwh = base_price .+ p.grid_price_addon_eur_mwh

    # A_V1: energy costs of first year [EUR]
    # Convert Wh * €/MWh → EUR
    A1 = sum(IN .* price_eur_mwh) * 1e-6

    a = annuity_factor(p.i_cap, p.T)
    b = price_change_factor(p.r_energy, p.i_cap, p.T)

    return A1 * a * b
end


############################################################
#  REVENUES — CONTROL ENERGY AND FEED-IN
############################################################

function revenue_control(sim::Union{Dict, OrderedDict}, p::VDIParams)
    haskey(sim, "Regelenergie_Wh") || return 0.0

    RE = sim["Regelenergie_Wh"]
    price = vecize_price(sim["RegelenergiePreis_EUR_MWh"], length(RE))

    A1 = sum(RE .* price) * 1e-6

    a = annuity_factor(p.i_cap, p.T)
    b = price_change_factor(p.r_rev, p.i_cap, p.T)

    return A1 * a * b
end

function revenue_feedin(sim::Union{Dict, OrderedDict}, p::VDIParams)
    OUT = sim["Grid_Out"]
    price = vecize_price(sim["FeedIn_price"], length(OUT))

    A1 = sum(OUT .* price) * 1e-6

    a = annuity_factor(p.i_cap, p.T)
    b = price_change_factor(p.r_rev, p.i_cap, p.T)

    return A1 * a * b
end


############################################################
#  MAIN FUNCTION — TOTAL ANNUITY AND HEAT PRICE
############################################################

function vdi2067_annuity(sim::Union{Dict, OrderedDict}, components::Vector{VDIComponent}, p::VDIParams)

    A_cap   = sum(capital_annuity(c, p) for c in components)
    A_op    = op_annuity(components, p)
    A_misc  = misc_annuity(components, p)
    A_energy = energy_annuity(sim, p)

    A_rev_control = revenue_control(sim, p)
    A_rev_feed    = revenue_feedin(sim, p)

    A_total = A_cap + A_op + A_misc + A_energy -
              (A_rev_control + A_rev_feed)

    Q_heat_kWh = sum(sim["Heat"]) / 1000.0
    heat_price = A_total / Q_heat_kWh

    return Dict(
        "A_cap" => A_cap,
        "A_op" => A_op,
        "A_misc" => A_misc,
        "A_energy" => A_energy,
        "A_rev_control" => A_rev_control,
        "A_rev_feed" => A_rev_feed,
        "A_total" => A_total,
        "heat_price_eur_per_kwh" => heat_price
    )
end

end # module
