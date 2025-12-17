module VDI2067
using OrderedCollections: OrderedDict

export VDIParams, VDIComponent, vdi2067_annuity,
       VDI_SCENARIO_NONE, VDI_SCENARIO_MOD, VDI_SCENARIO_PRO,
       heatpump_component, boiler_component, buffertank_component, battery_component

############################################################
#  PARAMETER STRUCTURE – VDI 2067 COMPATIBLE
############################################################

struct VDIParams
    T::Int                      # Evaluation period (years)
    i_cap::Float64              # interest factor for capital-related annuity
    r_cap::Float64              # Escalation of capital costs (investments, replacements)
    r_energy::Float64           # Escalation of energy prices (demand-dependent costs)
    r_op::Float64               # Escalation of maintenance & inspection
    r_inst::Float64             # Escalation of repair / reinstatement
    r_misc::Float64             # Escalation of miscellaneous costs (insurance, admin,…)
    r_rev::Float64              # Escalation of revenues (feed-in, control energy)
    grid_price_addon::Float64   # taxes, levies, grid fees
    AW_PV::Float64              # "Anzusetzender Wert" bei Direktvermarktung"
    AW_Wind::Float64            # "Anzusetzender Wert" bei Direktvermarktung"
end


############################################################
#  COMPONENT DEFINITION — COST GROUP 420/430
############################################################

struct VDIComponent
    name::String            # name of component (i.e. "HeatPump")
    A0::Float64             # initial investment at t=0 (EUR)
    TN::Int                 # technical lifetime (years)
    f_wartung::Float64      # maintenance + inspection factor (share of A0 in year 1)
    f_instand::Float64      # repair / reinstatement factor (share of A0 in year 1)
end


############################################################
#  COMPONENT CONSTRUCTORS FOR THE FOUR CONSIDERED COMPONENT TYPES
############################################################
# Only initial investment cost must be provided by parameterstudy.jl
# The technical component lifetime and maintenance/repair factors are hard-coded according to VDI specification and own assumptions
############################################################

"Heat pump: maintenance 1.5 %, repair 1.0 % of A0 in first year."
heatpump_component(A0::Real) =
    VDIComponent("HeatPump", float(A0), 20, 0.015, 0.01)

"Boiler / electric heater: maintenance 2.0 %, repair 1.0 %."
boiler_component(A0::Real) =
    VDIComponent("Boiler", float(A0), 15, 0.02, 0.01)

"Buffer tank: maintenance 0.5 %, repair 0.5 %."
buffertank_component(A0::Real) =
    VDIComponent("BufferTank", float(A0), 15, 0.005, 0.005)

"Battery: maintenance 1.0 %, repair 2.0 %."
battery_component(A0::Real) =
    VDIComponent("Battery", float(A0), 12, 0.01, 0.02)


############################################################
#  SCENARIO DEFINITIONS for different cost escalation (NONE / MODERATE / PROGRESSIVE)
############################################################

# Scenario 1 — no escalation
VDI_SCENARIO_NONE = VDIParams(
    20,       # Evaluation period (years)
    0.03,     # discount rate +3.0%
    0.000,    # no capital cost (Excel)
    0.000,    # no energy cost (Excel)
    0.000,    # no maintenance cost (Excel)
    0.000,    # no repair cost (Excel)
    0.000,    # no revenues
    0.000,    # no miscellaneous (Excel)
    #TODO grid fees based on peak power demand  
    120.0,     # grid price addon €/MWh
    79.90,     # AW_PV €/MWh
    7.35*1.42 # AW_Wind €/MWh
    )

# Scenario 2 — moderate escalation
VDI_SCENARIO_MOD = VDIParams(
    20,       # Evaluation period (years)
    0.03,     # discount rate +3.0%
    0.010,    # capital cost +1.2% (Excel)
    0.010,    # energy cost +1.0% (Excel)
    0.010,    # maintenance cost +0.5% (Excel)
    0.010,    # repair cost +0.5% (Excel)
    -0.01,    # revenues -10%
    0.010,    # miscellaneous +0.5% (Excel)
    #TODO grid fees based on peak power demand 
    120.0,     # grid price addon €/MWh
    # TODO escalations for AW_PV and AW_Wind?
    79.90,    # AW_PV €/MWh
    7.35*1.42 # AW_Wind €/MWh
    )

# Scenario 3 — progressive escalation
VDI_SCENARIO_PRO = VDIParams(
    20,       # Evaluation period (years)
    0.03,     # discount rate +3.0%
    0.020,    # capital +1.8% (Excel)
    0.020,    # energy cost +2.0% (Excel)
    0.020,    # maintenance +1.0% (Excel)
    0.020,    # repair +1.0% (Excel)
    -0.02,    # revenues -2.0%
    0.020,    # miscellaneous +1.0% (Excel)
    #TODO grid fees based on peak power demand 
    120.0,     # grid price addon €/MWh
    # TODO escalations for AW_PV and AW_Wind?
    79.90,     # AW_PV €/MWh
    7.35*1.42 # AW_Wind €/MWh
    )


############################################################
#  VDI 2067 MATHEMATICAL FUNCTIONS
############################################################

# Annuity factor 
# a = q^T * (q − 1) / (q^T − 1)
function annuity_factor(i, T)
    q = 1 + i
    T == 0 && return 0.0
    isapprox(q, 1.0) && return 1/T
    return (q^T * (q - 1)) / (q^T - 1)
end

# Price change factor 
# b = (1 − (qv / q)^T) / (q − qv)
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

function npv_replacements(A0, TN, T, i, r)      # net present value of component replacements
    npv = 0.0
    t = TN
    while t < T
        cost_t = A0 * (1+r)^t
        npv += cost_t / (1+i)^t
        t += TN
    end
    return npv
end

function residual_value(A0, TN, T, i, r)    # residual value at end of evaluation period
    q = 1 + i
    qv = 1 + r

    n = div(T, TN)                  # amount of component replacements in evaluation period
    t_last_repl = n * TN            # year of the last component replacement

    if t_last_repl == T
        return 0.0                       # replacement at the end of the evaluation period → no residual value
    end
    
    
    A_last = A0 * qv^t_last_repl    # cost of last replacement in the evaluation period

    write_off = ((n + 1) * TN - T) / TN     # function for linear write off / depreciation

    RW_T = A_last * write_off       # residual value of last replacement at the end of evaluation period
    
    return RW_T / q^T   # Discounting residual value to t=0
end

function capital_annuity(comp::VDIComponent, p::VDIParams)      # capital cost-related annuity
    return  (comp.A0 +
            npv_replacements(comp.A0, comp.TN, p.T, p.i_cap, p.r_cap) -
            residual_value(comp.A0, comp.TN, p.T, p.i_cap, p.r_cap)) * 
            annuity_factor(p.i_cap, p.T)
end


############################################################
#  OPERATING COSTS: OPERATION (B) + MAINTENANCE (IN)
############################################################

function op_annuity(components::Vector{VDIComponent}, p::VDIParams)     # operation cost-related annuity
    # First-year operating costs
    # TODO keep no operating costs (h/a and EUR/h) ?
    A_B1  = 0.0
    # First-year inspection, maintenance, repair and reinstatement cost
    A_IN1 = sum(c.A0 * (c.f_instand + c.f_wartung) for c in components)

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
    b = price_change_factor(p.r_misc, p.i_cap, p.T)

    return A1 * a * b
end


############################################################
#  ENERGY COSTS — TIME-RESOLVED → ANNUITY
############################################################

function energy_annuity(sim::Union{Dict, OrderedDict}, p::VDIParams)
    # TODO not only Grid_IN but also control_reserve?
    IN = sim["Grid_IN"]                                         # Wh time series
    base_price = vecize_price(sim["Grid_price"], length(IN))    # €/MWh (market price)

    # Add fixed components (taxes, grid fees, levies)
    # This is a linear addition: total electricity price = market + fixed components
    # TODO make more complex (e.g. peak-based grid fees)
    price_eur_mwh = base_price .+ p.grid_price_addon

    # A_V1: energy costs of first year [EUR]
    # Convert Wh * €/MWh → EUR
    A1 = sum(IN .* price_eur_mwh) * 1e-6

    a = annuity_factor(p.i_cap, p.T)
    b = price_change_factor(p.r_energy, p.i_cap, p.T)

    return A1 * a * b
end


############################################################
#  REVENUE — CONTROL ENERGY + CONTROL CAPACITY
############################################################
"""
    revenue_control(sim::Dict, p::VDIParams)

Calculates annualized revenues from control reserve provision according to VDI 2067.

Assumptions (model-consistent):
- Time resolution is fixed to 15 minutes (Δt = 0.25 h)
- All price profiles are provided at 15-minute resolution
- Control energy prices are given in EUR/MWh
- Control capacity prices are given in EUR/MW
- Offered control capacity is implicitly derived from the simulated
  reserve energy (Grid input for control reserve)
- Capacity prices remain constant over 4-hour blocks, but are provided
  as 15-minute profiles (already expanded)

Required entries in `sim`:
- "Control_energy"                :: Vector{Float64}  (Wh per timestep)
- "Control_energy_price"      :: Vector{Float64}  (EUR/MWh)
- "Control_power_price"      :: Vector{Float64}  (EUR/MW)

Returns:
- Annualized control reserve revenue [EUR/a]
"""
function revenue_control(sim::Dict, p::VDIParams)

    # Fixed model time step
    Δt = 0.25  # hours (15 minutes)

    # 1) CONTROL ENERGY (WORK) REVENUE
    A_E_energy = 0.0

    if haskey(sim, "Control_energy") && haskey(sim, "Control_energy_price")
        E = sim["Control_energy"]     # Wh
        price_E = vecize_price(sim["Control_energy_price"], length(E))

        # Wh * EUR/MWh → EUR
        A_E_energy = sum(E .* price_E) * 1e-6
    end

    # 2) CONTROL CAPACITY (POWER) REVENUE
    A_E_capacity = 0.0

    if haskey(sim, "Control_energy") && haskey(sim, "Control_power_price")
        E = sim["Control_energy"]     # Wh
        price_P = vecize_price(sim["Control_power_price"], length(E))

        # Offered capacity inferred from reserve energy:
        # P_t = E_t / Δt  → Wh / h = W → MW
        P_MW = (E ./ Δt) .* 1e-6

        # EUR/MW * MW * h → EUR
        A_E_capacity = sum(P_MW .* price_P .* Δt)
    end

    # 3) TOTAL FIRST-YEAR CONTROL RESERVE REVENUE
    A_E1 = A_E_energy + A_E_capacity

    # 4) ANNUITIZATION
    a = annuity_factor(p.i_cap, p.T)
    b = price_change_factor(p.r_rev, p.i_cap, p.T)

    return A_E1 * a * b
end

############################################################
#  REVENUES — FEED-IN
############################################################

function revenue_feedin(sim::Union{Dict, OrderedDict}, p::VDIParams)

    # 1) PHOTOVOLTAIC
    A_E_PV = 0.0

    if haskey(sim, "Grid_Out_PV")
        E_PV = sim["Grid_Out_PV"]                                  
        market_value_PV = vecize_price(sim["Market_price_PV"], length(E_PV))

        # Effective remuneration per timestep: max(MW, AW)
        price_eff_PV = max.(market_value_PV, p.AW_PV)

        # Wh * €/MWh → EUR
        A_E_PV = sum(E_PV .* price_eff_PV) * 1e-6
    end

    # 2) WIND
    A_E_Wind = 0.0

    if haskey(sim, "Grid_Out_Wind_Wh")
        E_Wind = sim["Grid_Out_Wind_Wh"]
        market_value_Wind = vecize_price(sim["Market_price_Wind"], length(E_Wind))

        # Effective remuneration per timestep: max(MW, AW)
        price_eff_Wind = max.(market_value_Wind, p.AW_Wind)

        # Wh * €/MWh → EUR
        A_E_Wind = sum(E_Wind .* price_eff_Wind) * 1e-6
    end

    # 3) TOTAL FIRST-YEAR FEED-IN REVENUE
    A_E2 = A_E_PV + A_E_Wind

    # 4) VDI 2067 ANNUITIZATION
    a = annuity_factor(p.i_cap, p.T)
    b = price_change_factor(p.r_rev, p.i_cap, p.T)

    return A_E2 * a * b
end


############################################################
#  MAIN FUNCTION — TOTAL ANNUITY AND HEAT PRICE
############################################################

function vdi2067_annuity(sim::Union{Dict, OrderedDict}, components::Vector{VDIComponent}, p::VDIParams)

    sim["Grid_IN"] = sim["m_power Grid_In->HeatPump"] .+ sim["m_power Grid_In->Boiler"] .+
                     sim["m_power Grid_In->Demand_Power"] .+ sim["m_power Grid_In->Battery"]
    sim["Grid_Out_PV"] = sim["m_power Photovoltaic->Grid_OUT"]
    sim["Grid_Out_Wind"] = sim["m_power WindFarm->Grid_OUT"]
    sim["Control_energy"] = sim["ControlReserve m_power OUT"]
    sim["Control_power_price"] =
    sim["Control_energy_price"] =
    sim["Grid_price"] =
    
    A_cap   = round(sum(capital_annuity(c, p) for c in components), digits=2)
    A_op    = round(op_annuity(components, p), digits=2)
    A_misc  = round(misc_annuity(components, p), digits=2)
    A_energy = round(energy_annuity(sim, p), digits=2)

    A_rev_control = round(revenue_control(sim, p), digits=2)
    A_rev_feed    = round(revenue_feedin(sim, p), digits=2)

    A_total = round(A_cap + A_op + A_misc + A_energy -
              (A_rev_control + A_rev_feed), digits=2)

    # calculating a heat price out of the total annual costs and the produced heat does not make sense
    # Q_heat_kWh = sum(sim["Heat"]) / 1000.0
    # heat_price = round(A_total / Q_heat_kWh, digits=2)   # EUR/kWh


    return Dict(
        "A_cap" => A_cap,
        "A_op" => A_op,
        "A_misc" => A_misc,
        "A_energy" => A_energy,
        "A_rev_control" => A_rev_control,
        "A_rev_feed" => A_rev_feed,
        "A_total" => A_total,
        # "heat_price_eur_per_kwh" => heat_price
    )
end

end # module
