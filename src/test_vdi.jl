###############################
#  VDI 2067 TEST PROGRAM
#  Run this directly in VS Code terminal
###############################

# Load your module (adjust the path if necessary)
include("VDI2067.jl")
using .VDI2067   # If module is defined inside included file

println("\n=== VDI 2067 TEST START ===\n")

############################################################
# 1. CREATE SIMPLE, COMPLETE INPUTS FOR TESTING
############################################################

# Example: 1 year, hourly resolution
N = 8760

# --- Fixed electricity price (€/MWh) ---
fixed_price_eur_mwh = 200.0   # = 0.20 €/kWh

# For testing: constant electrical consumption
grid_in_wh  = fill(1000.0, N)   # 1 kWh per hour
grid_out_wh = zeros(N)          # no export

# Heat demand (just for denominator)
heat_wh = fill(5000.0, N)       # 5 kWh heat per hour

# Required price vectors
feedin_price_eur_mwh = fill(0.0, N)
regelenergie_wh = zeros(N)
regelenergie_preis = zeros(N)

# Build simulation dictionary expected by vdi2067_annuity()
sim = Dict(
    "Grid_IN" => grid_in_wh,                       # Wh
    "Grid_Out" => grid_out_wh,                     # Wh feed-in
    "Grid_price" => fill(fixed_price_eur_mwh, N),  # €/MWh
    "Heat" => heat_wh,                             # Wh
    "FeedIn_price" => feedin_price_eur_mwh,        # €/MWh
    "Regelenergie_Wh" => regelenergie_wh,
    "RegelenergiePreis_EUR_MWh" => regelenergie_preis
)


############################################################
# 2. DEFINE TEST COMPONENTS (HEATPUMP, BOILER, BUFFER, BATTERY)
############################################################

components = [
    heatpump_component(25000, 20),      # A0=25k€, lifetime 20y
    boiler_component(5000, 15),         # A0=5k€
    buffertank_component(10000, 25),    # A0=10k€
    battery_component(8000, 12)         # A0=8k€, lifetime 12y
]


############################################################
# 3. RUN ALL THREE SCENARIOS
############################################################

function print_results(label, result)
    println("----- Scenario: $label -----")
    for (k,v) in result
        println(rpad(k, 28), " : ", v)
    end
    println()
end

# Scenario 1: no escalation
res_no = vdi2067_annuity(sim, components, VDI_SCENARIO_NO)
print_results("NO ESCALATION", res_no)

# Scenario 2: moderate escalation
res_mod = vdi2067_annuity(sim, components, VDI_SCENARIO_MOD)
print_results("MODERATE ESCALATION", res_mod)

# Scenario 3: progressive escalation
res_pro = vdi2067_annuity(sim, components, VDI_SCENARIO_PRO)
print_results("PROGRESSIVE ESCALATION", res_pro)

println("\n=== VDI 2067 TEST END ===\n")
