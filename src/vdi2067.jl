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
    grid_price_addon::Float64   # taxes, levies, grid fees  #TODO is this needed?
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
    f_bedien::Float64       # operation factor (h/a in year 1)
    incentive_p::Float64    # incentive percentage to be multiplied with A0
    incentive_max::Float64  # maximum incentive for single investment
end

# COMPONENT CONSTRUCTORS FOR THE FOUR CONSIDERED COMPONENT TYPES
# Only initial investment cost must be provided by parameterstudy.jl
# The technical component lifetime and maintenance/repair/operating factors are hard-coded according to VDI specification and own assumptions

"Heat pump: maintenance 1.5 %, repair 1.0 % of A0 in first year."
heatpump_component(A0::Real) =
    VDIComponent("HeatPump", float(A0), 20.0, 0.015, 0.01, 5.0, 0.2, 0.0)

"Boiler / electric heater: maintenance 2.0 %, repair 1.0 %."
boiler_component(A0::Real) =
    VDIComponent("Boiler", float(A0), 15.0, 0.02, 0.01, 5.0, 0.2, 0.0)

"Buffer tank: maintenance 0.5 %, repair 0.5 %."
buffertank_component(A0::Real) =
    VDIComponent("BufferTank", float(A0), 15.0, 0.005, 0.005, 0.0, 0.2, 0.0)

"Battery: maintenance 1.0 %, repair 2.0 %."
battery_component(A0::Real) =
    VDIComponent("Battery", float(A0), 12.0, 0.01, 0.02, 0.0, 0.0, 0.0)


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
    120.0,    # grid price addon €/MWh
    106.8,    # AW_PV €/MWh
    83.496    # AW_Wind €/MWh, 58.8*1.42
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
    #TODO grid fees based on historical average data
    120.0,     # grid price addon €/MWh
    # TODO escalations for AW_PV and AW_Wind?
    106.8,    # AW_PV €/MWh
    83.496    # AW_Wind €/MWh, 58.8*1.42
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
    #TODO grid fees based on historical average data
    120.0,     # grid price addon €/MWh
    # TODO escalations for AW_PV and AW_Wind?
    106.8,     # AW_PV €/MWh
    83.496    # AW_Wind €/MWh, 58.8*1.42
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

function capital_annuity_incentive(comp::VDIComponent, p::VDIParams)      # capital cost-related annuity with incentive considered
    A0_inc = max(comp.A0 * comp.incentive_p, comp.incentive_max)
    return  (A0_inc +
            npv_replacements(A0_inc, comp.TN, p.T, p.i_cap, p.r_cap) -
            residual_value(A0_inc, comp.TN, p.T, p.i_cap, p.r_cap)) * 
            annuity_factor(p.i_cap, p.T)
end


############################################################
#  OPERATING COSTS: OPERATION (B) + MAINTENANCE (IN)
############################################################

function op_annuity(components::Vector{VDIComponent}, p::VDIParams)     # operation cost-related annuity

    # first year operating costs (f_bedien in h/a, hourly labor cost -> 30 EUR/h)
    A_B1  = sum((c.f_bedien * 30.0) for c in components if c.A0 > 0)  # only use f_bedien if component is actually used -> p_th > 0

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

function energy_annuity(sim::Dict, p::VDIParams)
    IN = sim["Grid_IN"] .* 1e-6        # convert Wh time series in MWh
    base_price = 214.0        # vecize_price(sim["Grid_price"], length(IN))    # €/MWh (market price)

    # A_V1: energy costs of first year [EUR]
    # MWh * EUR/MWh → EUR
    A1 = sum(IN .* base_price)  

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
        E = sim["Control_energy"] .* 1e-6     # convert Wh in MWh
        price_E = vecize_price(sim["Control_energy_price"], length(E))

        # MWh * EUR/MWh → EUR
        A_E_energy = sum(E .* price_E)
    end

    # 2) CONTROL CAPACITY (POWER) REVENUE
    A_E_capacity = 0.0

    if haskey(sim, "Control_energy") && haskey(sim, "Control_power_price")
        E = sim["Control_energy"]     # Wh
        price_P = vecize_price(sim["Control_power_price"], length(E))

        # Offered capacity inferred from reserve energy:
        # P_t = E_t / Δt  → Wh / h = W → MW
        P_MW = (E ./ Δt) .* 1e-6

        # MW * EUR/MW → EUR
        A_E_capacity = sum(P_MW .* price_P)
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

"""
    revenue_feedin(sim::Dict, p::VDIParams)

Calculates annualized feed-in revenues for PV and wind plantsaccording to EEG.

Assumptions:
- Time resolution: 15 minutes (900 s)
- Feed-in is remunerated via market value + market premium
- Market premium is defined as: max(0, AW - market price)
- PV and wind are treated separately with individual market values
  and "anzulegender Wert" (AW)

Required entries in `sim`:
- "Grid_Out_PV"                :: Vector{Float64}
- "Grid_Out_Wind"              :: Vector{Float64}
- "Market_Value_PV"            :: Vector{Float64}
- "Market_Value_Wind"          :: Vector{Float64}

Required entries in `p`:
- p.AW_PV     :: €/MWh
- p.AW_Wind   :: €/MWh

Returns:
- Annualized feed-in revenue [EUR/a]
"""

function revenue_feedin(sim::Dict, p::VDIParams)

    # 1) PHOTOVOLTAIC
    A_E_PV = 0.0

    if haskey(sim, "Grid_Out_PV")
        E_PV = sim["Grid_Out_PV"] .* 1e-6      # Wh -> MWh                                  
        market_value_PV = vecize_price(sim["Market_Value_PV"], length(E_PV))

        # Effective remuneration per timestep: max(MW, AW)
        # the bigger number (market value or AW) per timestep is used
        price_eff_PV = max.(market_value_PV, p.AW_PV)

        # MWh * €/MWh → EUR
        A_E_PV = sum(E_PV .* price_eff_PV)
    end

    # 2) WIND
    A_E_Wind = 0.0

    if haskey(sim, "Grid_Out_Wind")
        E_Wind = sim["Grid_Out_Wind"] .* 1e-6     # Wh -> MWh  
        market_value_Wind = vecize_price(sim["Market_Value_Wind"], length(E_Wind))

        # Effective remuneration per timestep: max(MW, AW)
        # the bigger number (market value or AW) per timestep is used
        price_eff_Wind = max.(market_value_Wind, p.AW_Wind)

        # MWh * €/MWh → EUR
        A_E_Wind = sum(E_Wind .* price_eff_Wind)
    end

    # 3) TOTAL FIRST-YEAR FEED-IN REVENUE
    A_E2 = A_E_PV + A_E_Wind

    # 4) VDI 2067 ANNUITIZATION
    a = annuity_factor(p.i_cap, p.T)
    b = price_change_factor(p.r_rev, p.i_cap, p.T)

    return A_E2 * a * b
end

function co2_yearly(sim::Dict)
    IN = sim["Grid_IN"] + sim["Control_energy"] .* 1e-3        # convert Wh time series in kWh
    co2_intensity = vecize_price(sim["CO2_Grid"], length(IN))    # g/kWh 

    # what about PV and Wind?

    return sum(IN .* co2_intensity)  

end

############################################################
#  MAIN FUNCTION — TOTAL ANNUITY AND HEAT PRICE
############################################################

function vdi2067_annuity(sim::Union{Dict,OrderedDict}, components::Vector{VDIComponent}, p::VDIParams)
    sim_new = Dict()
    sim_new["Grid_IN"] = sim["m_power EnergyFlow Grid_IN->HeatPump"] .+ sim["m_power EnergyFlow Grid_IN->Boiler"] .+
                         sim["m_power EnergyFlow Grid_IN->Demand_Power"] .+ sim["m_power EnergyFlow Grid_IN->Battery"]
    sim_new["Grid_Out_PV"] = sim["m_power EnergyFlow Photovoltaic->Grid_OUT"]
    sim_new["Grid_Out_Wind"] = sim["m_power EnergyFlow WindFarm->Grid_OUT"]
    sim_new["Control_energy"] = sim["ControlReserve m_power OUT"]
    sim_new["Control_power_price"] = sim["Reserve_Power_Price"]
    sim_new["Control_energy_price"] = sim["Reserve_Energy_Price"]
    sim_new["Grid_price"] = sim["Stock_Price"]
    sim_new["Market_Value_PV"] = sim["Market_Price_PV"]
    sim_new["Market_Value_Wind"] = sim["Market_Price_Wind"]
    sim_new["CO2_Grid"] = sim["CO2_Grid"]

    sim = OrderedDict()

    A_cap   = sum(capital_annuity(c, p) for c in components)
    A_op    = op_annuity(components, p)
    A_misc  = misc_annuity(components, p)
    A_energy = energy_annuity(sim_new, p)
    A_cap_incentive = sum(capital_annuity_incentive(c, p) for c in components)
    A_rev_control = revenue_control(sim_new, p)
    A_rev_feed    = revenue_feedin(sim_new, p)

    A_total = A_cap + A_op + A_misc + A_energy -
              (A_rev_control + A_rev_feed)
              
    A_total_incentive = A_cap_incentive + A_op + A_misc + A_energy -
                        (A_rev_control + A_rev_feed)
    
    # CO2_yearly = co2_yearly(sim_new)


    # calculating a heat price out of the total annual costs and the produced heat does not make sense
        # Q_heat_kWh = sum(sim_new["Heat"]) / 1000.0
        # heat_price = round(A_total / Q_heat_kWh, digits=2)   # EUR/kWh


    return OrderedDict(
        "A_cap" => A_cap,
        "A_cap_incentive" => A_cap_incentive,
        "A_op" => A_op,
        "A_misc" => A_misc,
        "A_energy" => A_energy,
        "A_rev_control" => A_rev_control,
        "A_rev_feed" => A_rev_feed,
        "A_total" => A_total,
        "A_total_incentive" => A_total_incentive,
        # "CO2_yearly" => CO2_yearly
        # "heat_price_eur_per_kwh" => heat_price
    )
end

end # module
