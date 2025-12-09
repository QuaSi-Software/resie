"""
    vdi2067_annuity(sim, components; params=VDIParams())

Calculate according to VDI standard (simplified version):
- Capital-bound annual annuity (incl. replacements & residual value)
- Operating, demand-dependent costs, miscellaneous costs
- Revenues (grid feed-in, control energy commercialization)
- Total annuity and heat price (EUR/kWh)

Inputs:
- sim::Dict: Simulation dictionary with time series (Wh) and price profiles (EUR/MWh)
  Expected keys (at minimum):
    "Grid_IN" :: Vector{Float64}  # Grid purchase (Wh)
    "Grid_Out" :: Vector{Float64} # Grid feed-in (Wh)
    "Heat"    :: Vector{Float64}  # Heat demand (Wh)
    optional: "PV_selfconsumption" :: Vector{Float64}
    "Grid_price" either scalar or Vector{Float64} (EUR/MWh)
    "FeedIn_price" either scalar or Vector{Float64} (EUR/MWh)
    "ControlEnergy_price" Vector{Float64} (EUR/MWh)
    "ControlPower_price" Vector{Float64} (EUR/MW)  # When used with reserved power in MW multiply
    "Reserved_power" scalar or Vector (W) (Reserved power in W)
- components::Vector{VDIComponents}:
    Structures describing investment components (A0, lifetime, name)
- params::VDIParams: Interest / price change factors and default factors

Outputs:
- Dict with annuity components and heat price (€/kWh)
"""
module VDI2067

export VDIParams, VDIComponents, vdi2067_annuity

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
    include_replacements_at_horizon::Bool  # Control: whether to consider replacement purchases that occur exactly at T
end

VDIParams(; T=20,
          i_cap=0.03, r_cap=0.00,
          i_op=0.03, r_op=0.01,
          i_energy=0.03, r_energy=0.02,
          i_revenue=0.03, r_revenue=0.00,
          f_inst=0.01, f_winsp=0.015, spec_inspection_hours=5.0,
          grid_price_addon=0.0,
          include_replacements_at_horizon=false) = VDIParams(T, i_cap, r_cap, i_op, r_op, i_energy, r_energy, i_revenue, r_revenue, f_inst, f_winsp, spec_inspection_hours, grid_price_addon, include_replacements_at_horizon)

# Simple component definition for capital-related calculations
# TODO Parameters will probably be integrated into ReSiE Components in the future and can then be adjusted via input.json
struct VDIComponents
    name::String
    A0::Float64    # Investment costs at time 0 in EUR
    TN::Int        # Useful life in years
end

"""
Helper functions
""" 
# Annuity factor (VDI form): a = q^T * (q - 1) / (q^T - 1), with q = 1 + i. Limiting case i == 0 -> 1/T
function annuity_factor(i::Float64, T::Int)
    if T == 0
        return 0.0
    end
    q = 1.0 + i
    # Limiting case q == 1 => i == 0 -> Limit a = 1 / T
    if isapprox(q, 1.0; atol=1e-12)
        return 1.0 / T
    end
    return q^T * (q - 1.0) / (q^T - 1.0)
end

# Price change factor (transforms the first year amount into present value over T),
# according to VDI: b = (1 - (q_v / q)^T) / (q - q_v), where q = 1 + i, q_v = 1 + r
function price_change_factor(r, i, T)
    """
    Calculates the VDI price change factor:
      b = (1 - (q_v/q)^T) / (q - q_v)
    with q = 1 + i, q_v = 1 + r.
    Limiting case q == q_v -> b = T / q.
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

# Convert a scalar price to a vector if required
function vecize_price(price, N)
    if isa(price, Number)
        return fill(price, N)
    elseif length(price) == N
        return price
    else
        error("Price array length error: expected length $N, received $(length(price))")
    end
end

# Net present value of planned replacement purchases (excluding initial investment A0)
function npv_replacements(A0::Float64, TN::Int, T::Int, i::Float64, r::Float64)
     npv = 0.0
     if TN <= 0
         return npv
     end
     t = TN
     while t < T - 1e-9 # Replacement times at TN, 2*TN, ... strictly less than T
         price_at_purchase = A0 * (1 + r)^(t)
         npv += price_at_purchase / (1 + i)^t
         t += TN
     end
     return npv
 end

# Net present value of residual value at time T and discounting to t=0 (VDI-compliant)
function npv_restwert(A0::Float64, TN::Int, T::Int, i::Float64, r::Float64; include_repl_at_T::Bool=true)
     # VDI-compliant:
     # - Determine the last purchase made (before or at T)
     # - If replacement occurs exactly at T, this yields the full residual value = A0*(1+r)^T
     # - otherwise linear pro-rata residual value according to unused remaining useful life
     if TN <= 0 || T <= 0
         return 0.0
     end

     # Handle a replacement occurring exactly at T:
     n = floor(Int, T / TN)                # Number of complete life cycles
     last_purchase = n * TN

     if !include_repl_at_T && last_purchase == T && n > 0
         # Option: interpret replacement at T as after the observation period
         last_purchase = (n - 1) * TN
     end

     # last_purchase <= T (possibly equal to T)
     if last_purchase == T
         # Replacement occurs exactly at T -> full residual value (price escalation to T)
         rest_at_T = A0 * (1.0 + r)^T
     else
         # Last purchase occurs before T -> linear residual value
         age_at_T = T - last_purchase
         remaining_fraction = max(0.0, (TN - age_at_T) / TN)
         price_at_last_purchase = A0 * (1.0 + r)^(last_purchase)
         rest_at_T = price_at_last_purchase * remaining_fraction
     end

     npv = rest_at_T / (1.0 + i)^T
     return npv
 end

# Net present value of all capital-related costs (initial investment + replacement purchases - residual value)
function npv_total_capital(A0::Float64, TN::Int, T::Int, i::Float64, r::Float64)
    npv_repl = npv_replacements(A0, TN, T, i, r)
    npv_rest = npv_restwert(A0, TN, T, i, r)
    return A0 + npv_repl - npv_rest
 end

 # Calculate the annual capital annuity for a component
 function annuity_capital_component(comp::VDIComponents, params::VDIParams)
     npv = npv_total_capital(comp.A0, comp.TN, params.T, params.i_cap, params.r_cap)
     a = annuity_factor(params.i_cap, params.T)
     return npv * a
 end

 # Energy costs in the first year from the time series (EUR)
 function first_year_energy_cost(sim::Dict, params::VDIParams)
     # Expects: all energy data in Wh, prices in EUR/MWh
     Grid_IN = sim["Grid_IN"]
     N = length(Grid_IN)
     grid_price = vecize_price(sim["Grid_price"], N)
     if haskey(sim, "PV_selfconsumption")
         pv_self_consumption = sim["PV_selfconsumption"]
     else
         pv_self_consumption = zeros(N)
     end
     # Surcharge for taxes/levies (absolute surcharge in EUR/MWh)
     grid_price = [p + params.grid_price_addon for p in grid_price]
     # For own generation if a price profile is provided (e.g. marginal costs)
     # TODO not used further; either remove or deduct from costs; also used in first_year_revenues
     feedin_price = haskey(sim, "FeedIn_price") ? vecize_price(sim["FeedIn_price"], N) : zeros(N) 
     selfgen_price = haskey(sim, "SelfGen_price") ? vecize_price(sim["SelfGen_price"], N) : fill(0.0, N)

     # Costs for grid purchase:
     cost_grid = sum( Grid_IN .* grid_price ) * 1e-6   # EUR
     # Costs for self-consumed own generation (if marginal costs available)
     # TODO does self-consumption generate costs or are these already included in the PV system annuity or is this something else?
     cost_self = sum( pv_self_consumption .* selfgen_price ) * 1e-6
     return cost_grid + cost_self
 end

# Annualize energy costs
function annuity_energy(sim::Dict, params::VDIParams)
    E1 = first_year_energy_cost(sim, params)
    b = price_change_factor(params.r_energy, params.i_energy, params.T)
    a = annuity_factor(params.i_energy, params.T)
    return E1 * b * a
end

# First year revenues
function first_year_revenues(sim::Dict)
    N = length(sim["Grid_IN"])
    # Control energy revenues: product of grid purchase and control energy price profile
    control_energy_price = haskey(sim, "ControlEnergy_price") ? vecize_price(sim["ControlEnergy_price"], N) : zeros(N)
    revenue_control_energy = sum(sim["Grid_IN"] .* control_energy_price) * 1e-6
    # Control power revenues: price * reserved power (conversion W -> MW)
    control_power_price = haskey(sim, "ControlPower_price") ? vecize_price(sim["ControlPower_price"], N) : zeros(N)
    reserved_power_W = haskey(sim, "Reserved_power") ? sim["Reserved_power"] : 0.0
    reserved_power_MW = reserved_power_W * 1e-6
    revenue_control_power = sum(control_power_price .* reserved_power_MW) # Price in EUR/MW -> EUR
    # Feed-in tariff own generation:
    grid_out = haskey(sim, "Grid_Out") ? sim["Grid_Out"] : zeros(N)
    feed_in_price = haskey(sim, "FeedIn_price") ? vecize_price(sim["FeedIn_price"], N) : zeros(N)
    revenue_feedin = sum(grid_out .* feed_in_price) * 1e-6

    # Sum of revenues
    return revenue_control_energy + revenue_control_power + revenue_feedin
end

function annuity_revenues(sim::Dict, params::VDIParams)
    E1 = first_year_revenues(sim)
    b = price_change_factor(params.r_revenue, params.i_revenue, params.T)
    a = annuity_factor(params.i_revenue, params.T)
    return E1 * b * a
end

# Operating costs (maintenance etc.) based on percentage of investment costs or effort hours * hourly rate
# Optional 'op_bases' can be provided in the simulation dictionary; standard is percentage of investment
function first_year_operation_costs(components::Vector{VDIComponents}, params::VDIParams)
     # Standard: percentage share of investment costs
     total_cap = sum(c.A0 for c in components)
     # Use the parameter factors (f_inst & f_winsp)
     inst_cost = total_cap * params.f_inst
     winsp_cost = total_cap * params.f_winsp
     # additional basis (can be specified in sim dictionary)
     return inst_cost + winsp_cost
 end

function annuity_operation(components::Vector{VDIComponents}, params::VDIParams)
     # Basis of first year
     A1 = first_year_operation_costs(components, params)
     a = annuity_factor(params.i_op, params.T)
     b = price_change_factor(params.r_op, params.i_op, params.T)
     return A1 * a * b
 end

 # Miscellaneous costs: percentage of investment costs (insurance, administration, etc.)
function annuity_sonstige(components::Vector{VDIComponents}, params::VDIParams; pct=0.02)
     total_cap = sum(c.A0 for c in components)
     A1 = total_cap * pct
     a = annuity_factor(params.i_op, params.T)
     b = price_change_factor(params.r_op, params.i_op, params.T)
     return A1 * a * b
 end

 # Main function to calculate annuities according to VDI
function vdi2067_annuity(sim::Dict, components::Vector{VDIComponents}; params=VDIParams())
     # Calculate capital annuities for all components
     capital_annuity_per_comp = Dict{String, Float64}()
     for c in components
         capital_annuity_per_comp[c.name] = annuity_capital_component(c, params)
     end
     A_cap_total = sum(values(capital_annuity_per_comp))

     # Operating costs
     A_betrieb = annuity_operation(components, params)

     # Demand-based energy costs
     A_bedarf = annuity_energy(sim, params)

     # Miscellaneous costs
     A_sonst = annuity_sonstige(components, params)

     # Revenues
     A_erloese = annuity_revenues(sim, params)

     # Total annuity
     A_ges = A_cap_total + A_betrieb + A_bedarf + A_sonst - A_erloese

     # Heat price: share per kWh heat
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