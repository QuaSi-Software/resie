using Statistics
using Plots: Plots

"""
Implementation of a battery component holding electric charge.

For the moment the implementation remains simple with only one state (its charge) and one
parameter (its capacity).
"""
Base.@kwdef mutable struct Battery <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    m_el_in::Symbol
    m_el_out::Symbol
    m_heat_lt_out::Symbol
    model_type::String

    capacity::Float64
    load::Float64
    # simplified
    charge_efficiency::Floathing
    discharge_efficiency::Floathing
    self_discharge_rate::Float64
    max_charge_C_rate::Float64
    max_discharge_C_rate::Float64
    SOC_min::Float64
    SOC_max::Float64
    # detailed
    V_n_bat::Float64
    V_n::Float64
    r_i::Float64
    V_0::Float64
    K::Float64
    A::Float64
    B::Float64
    V_cell_min::Float64
    V_cell_max::Float64
    capacity_cell_Ah::Float64
    m::Float64
    alpha::Float64
    cell_cutoff_current::Float64
    # with aging
    cycles::Float64
    Temp::Float64
    k_qn::Array{Float64,1}
    k_qT::Array{Float64,1}
    k_n::Array{Float64,1}
    k_T::Array{Float64,1}
    I_ref::Float64
    T_ref::Float64

    load_end_of_last_timestep::Float64
    losses::Float64
    max_charge_energy::Float64
    max_discharge_energy::Float64
    V_cell::Float64
    V_cell_last::Float64
    V_cell_charge::Float64
    V_cell_discharge::Float64
    n_cell_p::Float64
    n_cell_s::Float64
    extracted_charge::Float64
    extracted_charge_last::Float64
    SOC::Float64

    heat_lt_is_usable::Bool
    heat_out::Float64

    process_done::Bool
    load_done::Bool

    function Battery(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        heat_lt_is_usable = default(config, "heat_lt_is_usable", false)

        m_el_in = Symbol(default(config, "m_el_in", "m_e_ac_230v"))
        m_el_out = Symbol(default(config, "m_el_out", "m_e_ac_230v"))
        m_heat_lt_out = Symbol(default(config, "m_heat_lt_out", "m_h_w_lt1"))

        register_media([m_el_in, m_el_out, m_heat_lt_out])
        output_interfaces = InterfaceMap(m_el_out => nothing)
        if heat_lt_is_usable
            output_interfaces[m_heat_lt_out] = nothing
        end

        # get model type from input file
        model_type = default(config, "model_type", "simplified")
        model_type_allowed_values = ["simplified", "detailed", "Li-LFP"]
        if !(model_type in model_type_allowed_values)
            @error "Undefined model type \"$(model_type)\" of battery \"$(uac)\". " *
                   "Has to be one of: $(model_type_allowed_values)."
            throw(InputError)
        end

        if model_type == "Li-LFP"
            config["V_n"] =  3.2
            config["r_i"] =  0.00016
            config["V_0"] =  3.36964
            config["K"] =  0.03546
            config["A"] =  0.08165
            config["B"] =  0.1003
            config["V_cell_min"] =  2.0
            config["V_cell_max"] =  3.4
            config["capacity_cell_Ah"] =  1090
            config["m"] = 1.0269
            config["alpha"] = -0.01212
            config["k_qn"] = [-1.27571e-7 1.22095e-11]
            config["k_qT"] = [1.32729e-3 -7.9763e-6]
            config["k_n"] = [9.71249e-6 7.51635e-4 -8.59363e-5 -2.92489e-4]
            config["k_T"] = [1.05135e-3 1.83721e-2 -7.72438e-3 -4.31833e-2]
            config["I_ref"] = 100
            config["T_ref"] = 25
        end

        # check that charge_efficiency and discharge_efficiency are set for simplified model
        charge_efficiency = default(config, "charge_efficiency", nothing)
        discharge_efficiency = default(config, "discharge_efficiency", nothing)
        if model_type == "simplified" && 
           (charge_efficiency == nothing || discharge_efficiency == nothing)
           # end of expression
            @error "If model \"simplified\" is used for battery \"$(uac)\" then " *
                   "\"charge_efficiency\" and \"discharge_efficiency\" must be given."
            throw(InputError)      
        end

        # NMC
        # Q_exp= 2.584
        # Q_full= 3.2
        # Q_nom= 3.126
        # V_exp= 3.53
        # V_full= 4.2
        # V_nom= 3.342
        # V_cut= 2.772
        # I= 0.64
        # V_n= 3.6
        # r_i= 0.001155

        # LFP  
        # Q_exp = 0.04
        # Q_full= 2.25
        # Q_nom = 2.0
        # V_exp = 4.05
        # V_full= 4.1
        # V_nom = 3.4
        # V_cut = 2.706
        # V_n = 3.6
        # I = 0.45
        # r_i = 0.002

        # A = V_full - V_exp
        # B = 3/Q_exp
        # K = ((V_full - V_nom + A*(exp(-B*Q_nom) - 1))*(Q_full - Q_nom))/(Q_nom)
        # V_0 = V_full + K + r_i*I - A

        return new(uac, # uac
                   Controller(default(config, "control_parameters", nothing)),
                   sf_storage, # sys_function
                   InterfaceMap(m_el_in => nothing), # input_interfaces
                   output_interfaces, # output_interfaces
                   m_el_in,
                   m_el_out,
                   m_heat_lt_out,
                   model_type, # model_type
                   config["capacity"], # capacity of the battery in Wh
                   config["load"], # energy load of the battery in Wh
                   default(config, "charge_efficiency", nothing),
                   default(config, "discharge_efficiency", nothing),
                   default(config, "self_discharge_rate", 0.0), # rate of self-discharge including stand-by losses in %/month (month=30days)
                   default(config, "max_charge_C_rate", 1.0), # maximum continuos charge power in W
                   default(config, "max_discharge_C_rate", 1.0), # maximum continuos discharge power in W
                   default(config, "SOC_min", 0),
                   default(config, "SOC_max", 0),
                   config["V_n_bat"], # nominal voltage of the battery pack
                   default(config, "V_n", 0.0), # nominal voltage of the cell in V
                   default(config, "r_i", 0.0), # internal resistance of the cell in Ohm
                   default(config, "V_0", 0.0), # battery constant voltage in V
                   default(config, "K", 0.0), # polarisation voltage in V
                   default(config, "A", 0.0), # exponential zone amplitude in V
                   default(config, "B", 0.0), # exponential zone time constant inverse in 1/(Ah)
                   default(config, "V_cell_min", 0.0), # minimum cell voltage that defines the battery as empty in V
                   default(config, "V_cell_max", 0.0), # maximum cell voltage that defines the battery as full in V
                   default(config, "capacity_cell_Ah", 0.0), # nominal cell capacity in Ah
                   default(config, "m", 1.0),
                   default(config, "alpha", 0.0),
                   default(config, "cell_cutoff_current", 0.003 * default(config, "capacity_cell_Ah", 0.0)),
                   default(config, "cycles", 1.0), # number of cycles
                   default(config, "T_bat", 25.0), # cell temperature in °C
                   default(config, "k_qn", [0.0, 0.0]),
                   default(config, "k_qT", [0.0, 0.0]),
                   default(config, "k_n", [0.0, 0.0, 0.0, 0.0]),
                   default(config, "k_T", [0.0, 0.0, 0.0, 0.0]),
                   default(config, "I_ref", 0.0),
                   default(config, "T_ref", 25.0),
                   0.0, # load_end_of_last_timestep
                   0.0, # losses
                   0.0, # maximum charge energy in current timestep
                   0.0, # maximum discharge energy in current timestep
                   0.0, # battery voltage in V
                   0.0, # battery voltage at the end of the last timestep in V
                   0.0, # battery voltage after charging
                   0.0, # battery voltage after discharging
                   1.0, # number of cells connected in parallel 
                   1.0, # number of cells connected in series
                   0.0, # sum of all added and removed charge to the cell in Ah
                   0.0, # sum of all added and removed charge to the cell in the last timestep in Ah
                   0.0, # SOC in the current time step
                   heat_lt_is_usable,
                   0.0, # heat output if the heat from the losses are usable
                   false, # process_done, bool indicating if the process step has already been performed in the current time step
                   false) # load_done, bool indicating if the load step has already been performed in the current time step
    end
end

function initialise!(unit::Battery, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.input_interfaces[unit.m_el_in],
                          unload_storages(unit.controller, unit.m_el_in))
    set_storage_transfer!(unit.output_interfaces[unit.m_el_out],
                          load_storages(unit.controller, unit.m_el_out))
    if unit.heat_lt_is_usable
        set_storage_transfer!(unit.output_interfaces[unit.m_heat_lt_out],
                              load_storages(unit.controller, unit.m_heat_lt_out))
    end

    unit.load_end_of_last_timestep = copy(unit.load)
    unit.self_discharge_rate *= sim_params["time_step_seconds"] / (30 * 24 * 3600)

    unit.SOC = (unit.load / unit.capacity) * 100
    if unit.SOC < unit.SOC_min || unit.SOC > unit.SOC_max
        @warn "Starting load for battery $(unit.uac) is not in the defined SOC range." *
              "The load is set to the closest valid value."
        if abs(unit.SOC_min - unit.SOC) < abs(unit.SOC_max - unit.SOC)
            unit.SOC = unit.SOC_min
            unit.load = unit.capacity * unit.SOC_min / 100
        else
            unit.SOC = unit.SOC_max
            unit.load = unit.capacity * unit.SOC_max / 100
        end
    end

    if unit.model_type != "simplified"
        # estimate cells in parallel and series in the battery
        unit.n_cell_s = unit.V_n_bat / unit.V_n
        unit.n_cell_p = (unit.capacity / unit.V_n_bat) / unit.capacity_cell_Ah
        if unit.n_cell_p < 1
            @warn "The battery capacity, battery voltage, cell voltage and cell capacity " *
                  "for $(unit.uac) don't fit together. There is less than 1 battery cell " *
                  "in parallel which may lead to wrong results close to maximum and " *
                  "minimum SOC. The suggestion is to reduce the battery voltage V_n_bat."
        end

        n = unit.cycles
        # f_I = (discharge_current / unit.I_ref)^unit.alpha
        f_n = (1 + unit.k_qn[1] * (n - 1) * n / 2 + unit.k_qn[2] * (n - 1) * n * (2 * n - 1) / 2)
        f_T = (1 + unit.k_qT[1] * (unit.Temp - unit.T_ref) + unit.k_qT[2] * (unit.Temp - unit.T_ref)^2)
        Q = unit.capacity_cell_Ah * 1 * f_n * f_T

        # calculate minimum and maximum cell voltages for SOC_max and SOC_min
        extracted_charge_empty = Q * (1 - (unit.SOC_min) / 100)
        unit.V_cell_min = f_V_cell(extracted_charge_empty, unit.I_ref, unit)
        if unit.V_cell_min < 0
            @error "The cells parameters are not correctly defined."
            throw(InputError)
        end

        extracted_charge_full = Q * (1 - unit.SOC_max / 100)
        unit.V_cell_max = f_V_cell(extracted_charge_full, -unit.cell_cutoff_current, unit)
        # unit.V_cell_max = f_V_cell(extracted_charge_full, -unit.max_discharge_C_rate * unit.capacity_cell_Ah, unit)

        # get starting point for simulation
        if unit.SOC == 100
            unit.extracted_charge = 0
            unit.V_cell = unit.V_cell_max
        elseif unit.SOC == 0
            unit.extracted_charge = Q
            unit.V_cell = unit.V_cell_min
        else
            unit.V_cell_last = unit.V_n
            unit.extracted_charge = Q * (1 - unit.SOC / 100)
            unit.V_cell = f_V_cell(unit.extracted_charge, unit.I_ref, unit)
        end
        unit.V_cell_last = copy(unit.V_cell)
        unit.extracted_charge_last = copy(unit.extracted_charge)
    end
end

function reset(unit::Battery)
    invoke(reset, Tuple{Component}, unit)

    unit.max_charge_energy = 0.0
    unit.max_discharge_energy = 0.0
    unit.losses = 0.0
end

function control(unit::Battery,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)
    if discharge_is_allowed(unit.controller, sim_params) && unit.SOC > unit.SOC_min &&
       (unit.V_cell > unit.V_cell_min || unit.model_type == "simplified")
       # end of expression
        discharge_current = unit.max_discharge_C_rate * unit.capacity_cell_Ah
        unit.discharge_efficiency,
        unit.V_cell_discharge,
        unit.max_discharge_energy,
        charge = calc_efficiency_current(discharge_current, unit, sim_params)
    else
        unit.max_discharge_energy = 0.0
    end
    set_max_energy!(unit.output_interfaces[unit.m_el_out], unit.max_discharge_energy)

    charge_current = unit.max_charge_C_rate * unit.capacity_cell_Ah
    if charge_is_allowed(unit.controller, sim_params) && unit.SOC < unit.SOC_max - 0.01
        charge_current = -unit.max_charge_C_rate * unit.capacity_cell_Ah
        unit.charge_efficiency,
        unit.V_cell_charge,
        unit.max_charge_energy,
        charge = calc_efficiency_current(charge_current, unit, sim_params)
        if unit.model_type != "simplified" && 
           -sim_params["wh_to_watts"](charge) < unit.cell_cutoff_current &&
           unit.V_cell_last == unit.V_cell_max
           # end of expression
            unit.max_charge_energy = 0.0
        end
    else
        unit.max_charge_energy = 0.0
    end
    set_max_energy!(unit.input_interfaces[unit.m_el_in], unit.max_charge_energy)

    if unit.heat_lt_is_usable
        set_max_energy!(unit.output_interfaces[unit.m_heat_lt_out], unit.heat_out, nothing, unit.Temp)
    end
end

function calc_efficiency(energy::Number, unit::Battery, sim_params::Dict{String,Any})
    if unit.model_type == "simplified"
        if energy >= 0
            return unit.discharge_efficiency, 1, min(energy, unit.load), 0
        else
            return unit.charge_efficiency, 1, min(abs(energy), unit.capacity - unit.load), 0
        end
    else
        # calculate V_cell
        cell_energy = energy / unit.n_cell_p / unit.n_cell_s
        V_min = unit.V_cell_min
        V_max = unit.V_cell_max
        asymp = 0.0
        # limit Voltage range to reasonable range to prevent find_zero to run into an asymptote
        if energy < 0
            V_min = max(unit.V_cell_last * 0.8, unit.V_cell_min)
        elseif energy > 0
            V_max = min(unit.V_cell_last * 1.2, unit.V_cell_max)
            try
                asymp = find_zero(e_cell -> denom(V_min, e_cell, unit, sim_params),
                                  (unit.SOC - unit.SOC_min) / 100 * unit.capacity_cell_Ah * unit.V_n / 2,
                                  Roots.Order1())
            catch
                asymp = find_zero(e_cell -> denom(V_min, e_cell, unit, sim_params),
                                  (unit.SOC - unit.SOC_min) / 100 * unit.capacity_cell_Ah * unit.V_n / 4,
                                  Roots.Order1())
            end
            if asymp < 0
                asymp = cell_energy * 1.1
            end
        end
        if energy == 0
            V_cell = find_zero(V_cell -> V_cell_function(V_cell, 0.0, sim_params["time_step_seconds"], unit),
                               unit.V_cell_last,
                               Roots.Order1())
            max_energy_cell = 0.0
        # check if the resulting voltage with given energy is between V_min & V_max
        elseif cell_energy < asymp &&
               sign(V_cell_function(V_min, cell_energy, sim_params["time_step_seconds"], unit)) !=
               sign(V_cell_function(V_max, cell_energy, sim_params["time_step_seconds"], unit))
            # end of expression
            V_cell = find_zero(V_cell -> V_cell_function(V_cell, cell_energy, sim_params["time_step_seconds"], unit),
                               (V_min, V_max),
                               Roots.Brent())
            max_energy_cell = cell_energy

        # if very little energy is added the voltage of the battery can drop and give a wrong result
        elseif energy < 0 && abs(cell_energy) < 0.1 * unit.capacity_cell_Ah
            V_cell = unit.V_cell_last
            max_energy_cell = cell_energy
        # if deeply discharged and little energy is added the voltage of the battery can drop below V_min
        elseif unit.V_cell <= unit.V_cell_min && energy < 0 && abs(cell_energy) < 0.1 * unit.capacity_cell_Ah
            V_cell = unit.V_cell_min
            max_energy_cell = cell_energy

        # if battery is discharged but less energy can be removed than given
        elseif energy > 0
            ul = min(asymp * 0.999, cell_energy)
            try
                max_energy_cell = find_zero(e_cell -> V_cell_function(V_min, e_cell,
                                                                      sim_params["time_step_seconds"], unit),
                                            (unit.capacity_cell_Ah * 0.0001, ul),
                                            Roots.Brent())
            catch
                max_energy_cell = 0.0
            end
            V_cell = V_min

        # if battery is charged but less energy can be charged than given
        elseif energy < 0
            if sign(V_cell_function(V_max, cell_energy, sim_params["time_step_seconds"], unit)) !=
               sign(V_cell_function(V_max, -unit.capacity_cell_Ah * 0.0001, sim_params["time_step_seconds"], unit))
               # end fo expression
                max_energy_cell = find_zero(e_cell -> V_cell_function(V_max, e_cell,
                                                                      sim_params["time_step_seconds"], unit),
                                            (cell_energy, -unit.capacity_cell_Ah * 0.0001),
                                            Roots.Brent())
            else
                try
                    max_energy_cell = find_zero(e_cell -> V_cell_function(V_max, e_cell,
                                                                          sim_params["time_step_seconds"], unit),
                                                -unit.capacity_cell_Ah * 0.01,
                                                Roots.Order1())
                    max_energy_cell = min(max_energy_cell, 0)
                catch
                    max_energy_cell = 0.0
                end
            end
            V_cell = V_max
        end

        # calculate losses and efficiency
        V_cell_avg = (unit.V_cell_last + V_cell) / 2
        charge = max_energy_cell / V_cell_avg
        internal_losses = unit.r_i * charge^2
        max_energy_bat = abs(max_energy_cell) * unit.n_cell_p * unit.n_cell_s
        efficiency = ifelse(max_energy_cell == 0, NaN, 1 - (internal_losses / abs(max_energy_cell)))
        return efficiency, V_cell, max_energy_bat, charge
    end
end

# function to find the voltage level, after charging with a constant power with aging 
# processes considered; positive energy = discharging, negative energy = charging
function V_cell_function(V_cell::Number, energy::Number, time_step::Number, unit::Battery)
    charge = 2 * energy / (unit.V_cell_last + V_cell)
    current = charge / (time_step / 3600)
    n = unit.cycles
    T = unit.Temp
    # functions describing the dependence of the parameters on cycles n and Temperature T
    if energy == 0
        f_I = 1
    else
        f_I = (abs(current) / unit.I_ref)^unit.alpha
    end
    f_n = (1 + unit.k_qn[1] * (n - 1) * n / 2 + unit.k_qn[2] * (n - 1) * n * (2 * n - 1) / 2)
    f_T = (1 + unit.k_qT[1] * (T - unit.T_ref) + unit.k_qT[2] * (T - unit.T_ref)^2)
    Q = unit.capacity_cell_Ah * f_I * f_n * f_T

    V_0 = unit.V_0 * (1 + unit.k_n[1] * (n - 1)) * (1 + unit.k_T[1] * (T - unit.T_ref))
    r_i = unit.r_i * (1 + unit.k_n[2] * (n - 1)) * (1 + unit.k_T[2] * (T - unit.T_ref))
    K = unit.K * (1 + unit.k_n[3] * (n - 1)) * (1 + unit.k_T[3] * (T - unit.T_ref))
    A = unit.A * (1 + unit.k_n[4] * (n - 1)) * (1 + unit.k_T[4] * (T - unit.T_ref))

    y = V_0 -
        K * (unit.m * Q) / (unit.m * Q - (unit.extracted_charge + charge)) +
        A * exp(-unit.B * (unit.extracted_charge + charge)) -
        r_i * current -
        V_cell
    return y
end

function denom(V_cell, energy, unit, sim_params)
    n = unit.cycles
    T = unit.Temp
    f_n = (1 + unit.k_qn[1] * (n - 1) * n / 2 + unit.k_qn[2] * (n - 1) * n * (2 * n - 1) / 2)
    f_T = (1 + unit.k_qT[1] * (T - unit.T_ref) + unit.k_qT[2] * (T - unit.T_ref)^2)
    f_I = (abs(2 * energy) / (unit.I_ref * (unit.V_cell_last + V_cell) * (sim_params["time_step_seconds"] / 3600)))^unit.alpha

    unit.m * unit.capacity_cell_Ah * f_I * f_n * f_T -
    (unit.extracted_charge + 2 * energy / (unit.V_cell_last + V_cell))
end

function calc_efficiency_current(current::Number, unit::Battery, sim_params::Dict{String,Any})
    if unit.model_type == "simplified"
        if current >= 0
            energy = sim_params["watt_to_wh"](unit.max_discharge_C_rate * unit.capacity)
            return unit.discharge_efficiency, 1, min(energy, unit.load), 0
        else
            energy = sim_params["watt_to_wh"](unit.max_charge_C_rate * unit.capacity)
            return unit.charge_efficiency, 1, min(abs(energy), unit.capacity - unit.load), 0
        end
    else
        # calculate V_cell
        V_cell = f_V_cell(unit.extracted_charge + sim_params["watt_to_wh"](current), current, unit)
        max_current = current

        # check if the resulting voltage with given energy is between V_min & V_max
        # if very little energy is added the voltage of the battery can drop and give a wrong result
        if current < 0 && abs(current) < 0.1 * unit.capacity_cell_Ah
            V_cell = unit.V_cell_last
        # if deeply discharged and little energy is added the voltage of the battery can drop below V_min
        elseif unit.V_cell <= unit.V_cell_min && current < 0 && abs(current) < 0.1 * unit.capacity_cell_Ah
            V_cell = unit.V_cell_min

        # if battery is discharged but less energy can be removed than given
        elseif V_cell < unit.V_cell_min
            asymp = find_zero(I_cell -> denom_current(I_cell, unit, sim_params),
                              (unit.SOC - unit.SOC_min) / 100 * unit.capacity_cell_Ah / 2,
                              Roots.Order1())
            ul = min(asymp * 0.999, current)
            try
                max_current = find_zero(I_cell -> V_cell_function_current(unit.V_cell_min, I_cell,
                                                                          sim_params["time_step_seconds"], unit),
                                        (unit.capacity_cell_Ah * 0.0001, ul),
                                        Roots.Brent())
            catch
                max_current = 0.0
            end
            V_cell = unit.V_cell_min

        # if battery is charged but less energy can be charged than given
        elseif V_cell > unit.V_cell_max
            if sign(V_cell_function_current(unit.V_cell_max, current, sim_params["time_step_seconds"], unit)) !=
               sign(V_cell_function_current(unit.V_cell_max, -unit.capacity_cell_Ah * 0.0001,
                                            sim_params["time_step_seconds"], unit))
               # end fo expression
                max_current = find_zero(I_cell -> V_cell_function_current(unit.V_cell_max, I_cell,
                                                                          sim_params["time_step_seconds"], unit),
                                        (current, -unit.capacity_cell_Ah * 0.0001),
                                        Roots.Brent())
            else
                try
                    max_current = find_zero(I_cell -> V_cell_function_current(unit.V_cell_max, I_cell,
                                                                              sim_params["time_step_seconds"], unit),
                                            -unit.capacity_cell_Ah * 0.01,
                                            Roots.Order1())
                    max_current = min(max_current, 0)
                catch
                    max_current = 0.0
                end
            end
            V_cell = unit.V_cell_max
        end

        # calculate losses and efficiency
        V_cell_avg = (unit.V_cell_last + V_cell) / 2
        charge = sim_params["watt_to_wh"](max_current)
        internal_losses = unit.r_i * charge^2
        max_energy_cell = charge * V_cell_avg
        max_energy_bat = abs(max_energy_cell) * unit.n_cell_p * unit.n_cell_s
        efficiency = ifelse(max_energy_cell == 0, NaN, 1 - (internal_losses / abs(max_energy_cell)))
        return efficiency, V_cell, max_energy_bat, charge
    end
end

# function to find the voltage level, after charging with a constant power with aging 
# processes considered; positive energy = discharging, negative energy = charging
function V_cell_function_current(V_cell::Number, current::Number, time_step::Number,
                                 unit::Battery)
    charge = current * (time_step / 3600)
    n = unit.cycles
    T = unit.Temp
    # functions describing the dependence of the parameters on cycles n and Temperature T
    if current == 0
        f_I = 1
    else
        f_I = (abs(current) / unit.I_ref)^unit.alpha
    end
    f_n = (1 + unit.k_qn[1] * (n - 1) * n / 2 + unit.k_qn[2] * (n - 1) * n * (2 * n - 1) / 2)
    f_T = (1 + unit.k_qT[1] * (T - unit.T_ref) + unit.k_qT[2] * (T - unit.T_ref)^2)
    Q = unit.capacity_cell_Ah * f_I * f_n * f_T

    V_0 = unit.V_0 * (1 + unit.k_n[1] * (n - 1)) * (1 + unit.k_T[1] * (T - unit.T_ref))
    r_i = unit.r_i * (1 + unit.k_n[2] * (n - 1)) * (1 + unit.k_T[2] * (T - unit.T_ref))
    K = unit.K * (1 + unit.k_n[3] * (n - 1)) * (1 + unit.k_T[3] * (T - unit.T_ref))
    A = unit.A * (1 + unit.k_n[4] * (n - 1)) * (1 + unit.k_T[4] * (T - unit.T_ref))

    y = V_0 -
        K * (unit.m * Q) / (unit.m * Q - (unit.extracted_charge + charge)) +
        A * exp(-unit.B * (unit.extracted_charge + charge)) -
        r_i * current -
        V_cell
    return y
end

function denom_current(current, unit, sim_params)
    n = unit.cycles
    T = unit.Temp
    f_n = (1 + unit.k_qn[1] * (n - 1) * n / 2 + unit.k_qn[2] * (n - 1) * n * (2 * n - 1) / 2)
    f_T = (1 + unit.k_qT[1] * (T - unit.T_ref) + unit.k_qT[2] * (T - unit.T_ref)^2)
    f_I = f_I = (abs(current) / unit.I_ref)^unit.alpha

    unit.m * unit.capacity_cell_Ah * f_I * f_n * f_T -
    (unit.extracted_charge + current * (sim_params["time_step_seconds"] / 3600))
end

function process(unit::Battery, sim_params::Dict{String,Any})
    if unit.max_discharge_energy >= sim_params["epsilon"]
        exchanges_out = balance_on(unit.output_interfaces[unit.m_el_out], unit.output_interfaces[unit.m_el_out].target)
        energy_demand = balance(exchanges_out) + energy_potential(exchanges_out)

        if energy_demand >= 0.0 # process is only concerned with moving energy to the target
            set_max_energy!(unit.output_interfaces[unit.m_el_out], 0.0)
            if unit.model_type != "simplified"
                unit.discharge_efficiency = NaN
            end
        else
            if abs(abs(energy_demand) - unit.max_discharge_energy) > sim_params["epsilon"]
                unit.discharge_efficiency,
                unit.V_cell,
                discharge_energy_bat,
                charge = calc_efficiency(abs(energy_demand), unit, sim_params)
            else
                discharge_energy_bat = unit.max_discharge_energy
                unit.V_cell = unit.V_cell_discharge
                charge = discharge_energy_bat / unit.n_cell_p / unit.n_cell_s * 2 / (unit.V_cell_last + unit.V_cell)
            end
            unit.losses += discharge_energy_bat * (1 - unit.discharge_efficiency)
            unit.extracted_charge += charge

            unit.load -= discharge_energy_bat
            add!(unit.output_interfaces[unit.m_el_out], discharge_energy_bat)
        end
    else
        set_max_energy!(unit.output_interfaces[unit.m_el_out], 0.0)
        if unit.model_type != "simplified"
            unit.discharge_efficiency = NaN
        end
    end
    handle_component_update!(unit, "process", sim_params)
end

function load(unit::Battery, sim_params::Dict{String,Any})
    if unit.max_charge_energy >= sim_params["epsilon"]
        exchanges_in = balance_on(unit.input_interfaces[unit.m_el_in], unit.input_interfaces[unit.m_el_in].source)
        energy_available = balance(exchanges_in) + energy_potential(exchanges_in)

        if energy_available <= 0.0 # load is only concerned with receiving energy from the source
            set_max_energy!(unit.input_interfaces[unit.m_el_in], 0.0)
            if unit.model_type != "simplified"
                unit.charge_efficiency = NaN
            end
        else
            if abs(energy_available - unit.max_charge_energy) > sim_params["epsilon"] || 
               unit.extracted_charge != unit.extracted_charge_last
               # end of expression
                unit.charge_efficiency,
                unit.V_cell,
                charge_energy_bat,
                charge = calc_efficiency(-energy_available, unit, sim_params)
            else
                charge_energy_bat = unit.max_charge_energy
                unit.V_cell = unit.V_cell_charge
                charge = -charge_energy_bat / unit.n_cell_p / unit.n_cell_s * 2 / (unit.V_cell_last + unit.V_cell)
            end
            unit.losses += charge_energy_bat * (1.0 / unit.charge_efficiency - 1)
            unit.extracted_charge += charge

            unit.load += charge_energy_bat
            sub!(unit.input_interfaces[unit.m_el_in], charge_energy_bat)
        end
    else
        set_max_energy!(unit.input_interfaces[unit.m_el_in], 0.0)
        if unit.model_type != "simplified"
            unit.charge_efficiency = NaN
        end
    end
    handle_component_update!(unit, "load", sim_params)
end

function handle_component_update!(unit::Battery, step::String, sim_params::Dict{String,Any})
    if step == "process"
        unit.process_done = true
    elseif step == "load"
        unit.load_done = true
    end
    if unit.process_done && unit.load_done
        # update component
        self_loss = 0.0
        if unit.load_end_of_last_timestep == unit.load ||
           unit.extracted_charge_last == unit.extracted_charge && unit.SOC > unit.SOC_min
            self_loss = unit.capacity * unit.self_discharge_rate # TODO fix loss calculation with respect to SOC and charge
            unit.losses += self_loss
        end

        if unit.model_type == "simplified"
            unit.load = max(unit.load - unit.losses, 0)
            unit.SOC = unit.load / unit.capacity * 100
        else
            # calculate current and check for different conditions
            charge_diff = unit.extracted_charge - unit.extracted_charge_last
            current = sim_params["wh_to_watts"](charge_diff)
            unit.cycles += abs(charge_diff) / unit.capacity_cell_Ah / 2
            # if no current is flowing just consider self_discharge
            if current == 0
                if (unit.load - unit.losses) < 0
                    unit.load = 0
                else
                    unit.load -= unit.losses
                    unit.extracted_charge += self_loss / unit.n_cell_p / unit.n_cell_s / unit.V_n
                end
                unit.SOC = unit.load / unit.capacity * 100

            # check for constant voltage charging cutoff point 
            # CV charging current cutoff is set at 0.3% of the nominal cell capacity
            elseif unit.V_cell == unit.V_cell_max && abs(current) < unit.cell_cutoff_current
                unit.SOC = unit.SOC_max
                unit.load = unit.capacity * unit.SOC / 100
                if unit.extracted_charge < 0
                    unit.extracted_charge = 0
                end

            # normal calculation if current is not 0
            else
                n = unit.cycles
                f_I = (abs(current) / unit.I_ref)^unit.alpha
                f_n = (1 + unit.k_qn[1] * (n - 1) * n / 2 + unit.k_qn[2] * (n - 1) * n * (2 * n - 1) / 2)
                f_T = (1 + unit.k_qT[1] * (unit.Temp - unit.T_ref) + unit.k_qT[2] * (unit.Temp - unit.T_ref)^2)
                Q = unit.capacity_cell_Ah * f_I * f_n * f_T

                if unit.extracted_charge > Q
                    unit.extracted_charge = Q
                elseif unit.extracted_charge < 0
                    unit.extracted_charge = 0
                end

                unit.SOC = ((unit.m * Q) - unit.extracted_charge) / (unit.m * Q) * 100
                unit.load = unit.capacity * unit.SOC / 100
            end
            unit.V_cell_last = copy(unit.V_cell)
        end

        if unit.heat_lt_is_usable
            add!(unit.output_interfaces[unit.m_heat_lt_out], unit.heat_out, nothing, unit.Temp)
            unit.heat_out = unit.losses
            unit.losses = 0.0
        end

        unit.load_end_of_last_timestep = copy(unit.load)
        unit.extracted_charge_last = copy(unit.extracted_charge)
        # reset
        unit.process_done = false
        unit.load_done = false
    end
end

function plot_optional_figures_begin(unit::Battery, output_path::String, output_formats::Vector{String},
                                     sim_params::Dict{String,Any})::Bool
    # discharge curves of the cell chemistry at different discharge rates
    Q_max = find_zero(Q -> f_V_cell(Q, unit.max_discharge_C_rate * unit.capacity_cell_Ah, unit, 1, 50) - unit.V_cell_min,
                      unit.capacity_cell_Ah) * 1.001
    x = 0:(Q_max / 500):Q_max
    y1 = fill(NaN, length(x))
    y2 = fill(NaN, length(x))
    y3 = fill(NaN, length(x))
    for (idx, Q) in enumerate(x)
        y1[idx] = f_V_cell(Q, unit.max_discharge_C_rate * unit.capacity_cell_Ah, unit, unit.cycles, unit.T_ref)
        if y1[idx] <= unit.V_cell_min
            break
        end
    end
    for (idx, Q) in enumerate(x)
        y2[idx] = f_V_cell(Q, unit.max_discharge_C_rate * 2 * unit.capacity_cell_Ah, unit, unit.cycles, unit.T_ref)
        if y2[idx] <= unit.V_cell_min
            break
        end
    end
    for (idx, Q) in enumerate(x)
        y3[idx] = f_V_cell(Q, unit.max_discharge_C_rate / 2 * unit.capacity_cell_Ah, unit, unit.cycles, unit.T_ref)
        if y3[idx] <= unit.V_cell_min
            break
        end
    end
    p = Plots.plot(x, [y1, y2, y3], 
                   labels=["$(unit.max_discharge_C_rate)C" "$(unit.max_discharge_C_rate*2)C" "$(unit.max_discharge_C_rate/2)C"], 
                   lw=3)
    Plots.plot!(p, ; title="Battery discharge curve for $(unit.uac) at different C-rates",
                xlabel="Removed Cell Charge / Ah",
                ylabel="Cell Voltage / V",
                ylims=(unit.V_cell_min, unit.V_cell_max),
                xlims=(0, Q_max),
                minorticks=10,
                size=(1800, 1200),
                titlefontsize=30,
                guidefontsize=24,
                tickfontsize=24,
                legendfontsize=24,
                grid=true,
                minorgrid=true,
                gridlinewidth=1,
                margin=15Plots.mm)
    fig_name = "Battery_$(unit.uac)_discharge_curves_current"
    for output_format in output_formats
        Plots.savefig(p, output_path * "/" * fig_name * "." * output_format)
    end

    # discharge curves of the cell chemistry at different Temperatures
    y1 = fill(NaN, length(x))
    y2 = fill(NaN, length(x))
    y3 = fill(NaN, length(x))
    for (idx, Q) in enumerate(x)
        y1[idx] = f_V_cell(Q, unit.max_discharge_C_rate * unit.capacity_cell_Ah, unit, unit.cycles, unit.T_ref)
        if y1[idx] <= unit.V_cell_min
            break
        end
    end
    for (idx, Q) in enumerate(x)
        y2[idx] = f_V_cell(Q, unit.max_discharge_C_rate * unit.capacity_cell_Ah, unit, unit.cycles, 0)
        if y2[idx] <= unit.V_cell_min
            break
        end
    end
    for (idx, Q) in enumerate(x)
        y3[idx] = f_V_cell(Q, unit.max_discharge_C_rate * unit.capacity_cell_Ah, unit, unit.cycles, 50)
        if y3[idx] <= unit.V_cell_min
            break
        end
    end
    p = Plots.plot(x, [y1, y2, y3], 
                   labels=["$(unit.T_ref) °C" "0 °C" "50 °C"], 
                   lw=3)
    Plots.plot!(p, ; title="Battery discharge curve for $(unit.uac) at different C-rates",
                xlabel="Removed Cell Charge / Ah",
                ylabel="Cell Voltage / V",
                ylims=(unit.V_cell_min, unit.V_cell_max),
                xlims=(0, Q_max),
                minorticks=10,
                size=(1800, 1200),
                titlefontsize=30,
                guidefontsize=24,
                tickfontsize=24,
                legendfontsize=24,
                grid=true,
                minorgrid=true,
                gridlinewidth=1,
                margin=15Plots.mm)
    fig_name = "Battery_$(unit.uac)_discharge_curves_temperature"
    for output_format in output_formats
        Plots.savefig(p, output_path * "/" * fig_name * "." * output_format)
    end

    # discharge curves of the cell chemistry after different number of cycles
    y1 = fill(NaN, length(x))
    y2 = fill(NaN, length(x))
    y3 = fill(NaN, length(x))
    for (idx, Q) in enumerate(x)
        y1[idx] = f_V_cell(Q, unit.max_discharge_C_rate * unit.capacity_cell_Ah, unit, unit.cycles, unit.T_ref)
        if y1[idx] <= unit.V_cell_min
            break
        end
    end
    for (idx, Q) in enumerate(x)
        y2[idx] = f_V_cell(Q, unit.max_discharge_C_rate * unit.capacity_cell_Ah, unit, 3000, unit.T_ref)
        if y2[idx] <= unit.V_cell_min
            break
        end
    end
    for (idx, Q) in enumerate(x)
        y3[idx] = f_V_cell(Q, unit.max_discharge_C_rate * unit.capacity_cell_Ah, unit, 2000, unit.T_ref)
        if y3[idx] <= unit.V_cell_min
            break
        end
    end
    p = Plots.plot(x, [y1, y2, y3], 
                   labels=["$(unit.cycles) cycles" "3000 cycles" "2000 cycles"], 
                   lw=3)
    Plots.plot!(p, ; title="Battery discharge curve for $(unit.uac) after different number of cycles",
                xlabel="Removed Cell Charge / Ah",
                ylabel="Cell Voltage / V",
                ylims=(unit.V_cell_min, unit.V_cell_max),
                xlims=(0, Q_max),
                minorticks=10,
                size=(1800, 1200),
                titlefontsize=30,
                guidefontsize=24,
                tickfontsize=24,
                legendfontsize=24,
                grid=true,
                minorgrid=true,
                gridlinewidth=1,
                margin=15Plots.mm)
    fig_name = "Battery_$(unit.uac)_discharge_curves_cycles"
    for output_format in output_formats
        Plots.savefig(p, output_path * "/" * fig_name * "." * output_format)
    end

    return true
end

"""
I_1, I_2
V_full (P1)
V_2, Q_2 in exponential zone (P2)
V_3, Q_3 in exponential zone Q_3 = 2*Q_2 (P3)
V_4, Q_4 near end of nominal zone (P4)
V_cut, Q_full_1 end of curve for I_1 (P5)
Q_full_1, Q_full_2 capacities for I_1 and I_2 at end of curve
n_1, n_2 for n_2 > n_1 > 1
Q_full_1, Q_full_n_2 capacities at n_1 and n_2 cycles at reference Temperature at I_n and n_2 > n_1 > 1
V_full, V_full_n_2 Voltages at n_1 and n_2 cycles at reference Temperature at I_n and n_2 > n_1 > 1
V_nom_n_2, Q_nom_n_2 in nominal zone of curve for n_2 (P8)
T_1, T_2, T_ref for T_1 != T_2 != T_ref
Q_full_1, Q_full_T_2 capacities at T_1 and T_2 temperature at 1 cycle at I_T
V_full, V_full_T_2 voltages at T_1 and T_2 temperature at 1 cycle at I_T
V_nom_T_2, Q_nom_T_2 in nominal zone of curve for T_2 (P11)
r_i, 
V_6, Q_6 in nominal zone of curve I_2 only necessary if r_i not given

example Paper SP-LFP1000AHA:
calc_cell_values(300, 3.4, 3.318, 7.56, 3.301, 15.12, 3.172, 820, 2.0, 1070.4,
calc_cell_values(300, 3.4, 3.301, 15, 3.284, 30, 3.172, 820, 2.0, 1070.4,
calc_cell_values(100, 3.4, 3.346, 10, 3.332, 20, 3.216, 820, 2.0, 1087.5,
                 1000, 1061.5,
                 nothing, 3.141, 550,
                 500, 1000, 1026,
                 3000, 3.2, 3.15, 400, 806,
                 500, 25, 1002,         
                 55, 3.31, 3.277, 696, 1035,
                 25)

paper reversed:
V_full = 3.4
I_1 = 100
Q_full_1 = 1090
AB_5 = 1.31835
V_cut = 2.0
V_3 = 3.3293
Q_3 = 20
V_2 = 3.3483
Q_2 = 10
V_4 = 3.2295
Q_4 = 800
I_2 = 1000
Q_full_2 = 1060
calc_cell_values(100, 3.4, 3.3483, 10, 3.3293, 20, 3.2295, 800, 2.0, 1090,
                 1000, 1060,
                 nothing, 3.141, 550,
                 500, 1000, 1027,
                 3000, 3.2, 3.15, 400, 809,
                 500, -20, 987,
                 55, 3.4, 3.2782, 757.57, 1102, 25)
"""
function calc_cell_values(I_1, V_full, V_2, Q_2, V_3, Q_3, V_4, Q_4, V_cut, Q_full_1,
                          I_2, Q_full_2,
                          r_i::Union{Number,Nothing}=nothing, V_6=0.0, Q_6=0.0,
                          I_n=0.0, n_1=0.0, Q_full_n_1=0.0,
                          n_2::Union{Number,Nothing}=nothing, V_full_n_2=0.0, V_nom_n_2=0.0, Q_nom_n_2=0.0,
                          Q_full_n_2=0.0,
                          I_T=0.0, T_1=0.0, Q_full_T_1=0.0,
                          T_2::Union{Number,Nothing}=nothing, V_full_T_2=0.0, V_nom_T_2=0.0, Q_nom_T_2=0.0,
                          Q_full_T_2=0.0,
                          T_ref=25)

    # basic values new
    B = -1 / Q_2 * log((V_full - V_3) / (V_full - V_2) - 1)
    A = (V_full - V_3) / (1 - exp(-B * Q_3))

    AB_4 = V_full - V_4 - A * (1 - exp(-B * Q_4)) # exp(-B*Q_4) -> 0
    AB_5 = V_full - V_cut - A * (1 - exp(-B * Q_full_1)) # exp(-B*Q_full) -> 0

    m = (1 - AB_5 / AB_4) / (1 - (AB_5 * Q_4) / (AB_4 * Q_full_1)) * (Q_4 / Q_full_1)
    K = AB_5 * ((m * Q_full_1 / Q_full_1) - 1)
    if isnothing(r_i)
        r_i = 1 / (I_1 - I_2) * (V_6 + K * (m * Q_full_2 / (m * Q_full_2 - Q_6)) - A * exp(-B * Q_6) -
                                 V_3 - K * (m * Q_full_1 / (m * Q_full_1 - Q_3)) + A * exp(-B * Q_3))
    end
    V_0 = V_full + K + r_i * I_1 - A
    alpha = log(Q_full_2 / Q_full_1) / log(I_2 / I_1)
    I_ref = I_1
    Q_ref = Q_full_1

    # basic values old
    # A = V_full - V_exp
    # B = 3/Q_exp
    # K = ((V_full - V_nom + A*(exp(-B*Q_nom) - 1))*(Q_full_1 - Q_nom))/(Q_nom)
    # V_0 = V_full + K + r_i*I_1 - A

    # for aging if the needed values are provided
    k_qn = [0.0, 0.0]
    k_n = [0.0, 0.0, 0.0, 0.0]
    if !isnothing(n_2)
        N = [(n_1 - 1) * n_1/2 (n_1 - 1) * n_1 * (2 * n_1 - 1)/2
             (n_2 - 1) * n_2/2 (n_2 - 1) * n_2 * (2 * n_2 - 1)/2]
        q = [Q_full_n_1 / (Q_ref * (I_n / I_ref)^alpha) - 1
             Q_full_n_2 / (Q_ref * (I_n / I_ref)^alpha) - 1]
        # k_qn[1] = (q[2] - q[1]*N[2,2]/N[1,2]) / (N[2,1] - N[1,1]*N[2,2]/N[1,2])
        # k_qn[2] = (q[1] - k_qn[1]*N[1,1]) / N[1,2]
        k_qn = N \ q

        Q_n = Q_ref * (I_n / I_ref)^alpha *
              (1 + k_qn[1] * (n_2 - 1) * n_2 / 2 + k_qn[2] * (n_2 - 1) * n_2 * (2 * n_2 - 1) / 2)
        b_n = [V_full - (V_0 - r_i * I_1 - K + A)
               V_full_n_2 - (V_0 - r_i * I_n - K + A)
               V_nom_n_2 - (V_0 - r_i * I_n - K * (m * Q_n / (m * Q_n - Q_nom_n_2)) + A * exp(-B * Q_nom_n_2))
               V_cut - (V_0 - r_i * I_n - K * (m * Q_n / (m * Q_n - Q_full_n_2)) + A * exp(-B * Q_full_n_2))]
        G_n = [V_0 -r_i*I_1 -K A
               V_0*(n_2 - 1) -r_i*I_n*(n_2-1) -K*(n_2 - 1) A*(n_2 - 1)
               V_0*(n_2 - 1) -r_i*I_n*(n_2-1) -K*(n_2-1)*(m * Q_n/(m * Q_n - Q_nom_n_2)) A*(n_2-1)*exp(-B * Q_nom_n_2)
               V_0*(n_2 - 1) -r_i*I_n*(n_2-1) -K*(n_2-1)*(m * Q_n/(m * Q_n - Q_full_n_2)) A*(n_2-1)*exp(-B * Q_full_n_2)]
        k_n = G_n \ b_n
    end

    k_qT = [0.0, 0.0]
    k_T = [0.0, 0.0, 0.0, 0.0]
    if !isnothing(T_2)
        T = [T_1-T_ref (T_1 - T_ref)^2
             T_2-T_ref (T_2 - T_ref)^2]
        q = [Q_full_T_1 / (Q_ref * (I_T / I_ref)^alpha) - 1
             Q_full_T_2 / (Q_ref * (I_T / I_ref)^alpha) - 1]
        # k_qT[1] = (q[2] - q[1]*T[2,2]/T[1,2]) / (T[2,1] - T[1,1]*T[2,2]/T[1,2])
        # k_qT[2] = (q[1] - k_qn[1]*T[1,1]) / T[1,2]
        k_qT = T \ q

        Q_T = Q_ref * (I_T / I_ref)^alpha *
              (1 + k_qT[1] * (T_2 - T_ref) + k_qT[2] * (T_2 - T_ref)^2)
        b_T = [V_full - (V_0 - r_i * I_1 - K + A)
               V_full_T_2 - (V_0 - r_i * I_T - K + A)
               V_nom_T_2 - (V_0 - r_i * I_T - K * (m * Q_T / (m * Q_T - Q_nom_T_2)) + A * exp(-B * Q_nom_T_2))
               V_cut - (V_0 - r_i * I_T - K * (m * Q_T / (m * Q_T - Q_full_T_2)) + A * exp(-B * Q_full_T_2))]
        G_T = [V_0 -r_i*I_1 -K A
               V_0*(T_2 - T_ref) -r_i*I_T*(T_2-T_ref) -K*(T_2 - T_ref) A*(T_2 - T_ref)
               V_0*(T_2 - T_ref) -r_i*I_T*(T_2-T_ref) -K*(T_2-T_ref)*(m * Q_T/(m * Q_T - Q_nom_T_2)) A*(T_2-T_ref)*exp(-B * Q_nom_T_2)
               V_0*(T_2 - T_ref) -r_i*I_T*(T_2-T_ref) -K*(T_2-T_ref)*(m * Q_T/(m * Q_T - Q_full_T_2)) A*(T_2-T_ref)*exp(-B * Q_full_T_2)]
        k_T = G_T \ b_T
    end

    return V_0, K, A, B, r_i, m, alpha, k_qn, k_qT, k_n, k_T, I_ref, T_ref
end

# calculate the voltage level, after charging with a constant current defined from charging 
# power and voltage level before charging, positive energy = discharging, negative energy = charging
# changing voltage levels during charging are ignored
function calc_V_cell_cc_aging(charge::Number, I::Number, n, T, unit::Battery)
    # functions describing the dependence of the parameters on cycles n and Temperature T
    Q = unit.capacity_cell_Ah * (abs(I) / unit.I_ref)^unit.alpha *
        (1 + unit.k_qn[1] * (n - 1) * n / 2 + unit.k_qn[2] * (n - 1) * n * (2 * n - 1) / 2) *
        (1 + unit.k_qT[1] * (T - unit.T_ref) + unit.k_qT[2] * (T - unit.T_ref)^2)
    V_0 = unit.V_0 * (1 + unit.k_n[1] * (n - 1)) * (1 + unit.k_T[1] * (T - unit.T_ref))
    r_i = unit.r_i * (1 + unit.k_n[2] * (n - 1)) * (1 + unit.k_T[2] * (T - unit.T_ref))
    K = unit.K * (1 + unit.k_n[3] * (n - 1)) * (1 + unit.k_T[3] * (T - unit.T_ref))
    A = unit.A * (1 + unit.k_n[4] * (n - 1)) * (1 + unit.k_T[4] * (T - unit.T_ref))

    V_0 -
    K * (unit.m * Q) / (unit.m * Q - (unit.extracted_charge + charge)) +
    A * exp(-unit.B * (unit.extracted_charge + charge)) -
    r_i * I
end

# calculate the voltage level, after charging with a constant current defined from charging 
# power and voltage level before charging, positive energy = discharging, negative energy = charging
# changing voltage levels during charging are ignored
function f_V_cell(extracted_charge::Number, current::Number, unit::Battery, n::Number=1,
                  T::Number=25)
    if current == 0
        f_I = 1
    else
        f_I = (abs(current) / unit.I_ref)^unit.alpha
    end
    f_n = (1 + unit.k_qn[1] * (n - 1) * n / 2 + unit.k_qn[2] * (n - 1) * n * (2 * n - 1) / 2)
    f_T = (1 + unit.k_qT[1] * (T - unit.T_ref) + unit.k_qT[2] * (T - unit.T_ref)^2)
    Q = unit.capacity_cell_Ah * f_I * f_n * f_T

    V_0 = unit.V_0 * (1 + unit.k_n[1] * (n - 1)) * (1 + unit.k_T[1] * (T - unit.T_ref))
    r_i = unit.r_i * (1 + unit.k_n[2] * (n - 1)) * (1 + unit.k_T[2] * (T - unit.T_ref))
    K = unit.K * (1 + unit.k_n[3] * (n - 1)) * (1 + unit.k_T[3] * (T - unit.T_ref))
    A = unit.A * (1 + unit.k_n[4] * (n - 1)) * (1 + unit.k_T[4] * (T - unit.T_ref))

    V_0 -
    K * (unit.m * Q) / (unit.m * Q - extracted_charge) +
    A * exp(-unit.B * extracted_charge) -
    r_i * current
end

function output_values(unit::Battery)::Vector{String}
    channels = [string(unit.m_el_in) * ":IN",
                string(unit.m_el_out) * ":OUT",
                "Load",
                "Load%",
                "Capacity",
                "LossesGains",
                "charge_efficiency",
                "discharge_efficiency",
                "CellVoltage",
                "SOC",
                "ExtractedCharge",
                "Cycles",
                "Temperature"]

    if unit.heat_lt_is_usable
        append!(channels, [string(unit.m_heat_lt_out) * ":OUT"])
        return channels
    else
        return channels
    end
end

function output_value(unit::Battery, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "Load"
        return unit.load
    elseif key.value_key == "Load%"
        return 100 * unit.load / unit.capacity
    elseif key.value_key == "Capacity"
        return unit.capacity
    elseif key.value_key == "LossesGains"
        return -unit.losses
    elseif key.value_key == "charge_efficiency"
        return unit.charge_efficiency
    elseif key.value_key == "discharge_efficiency"
        return unit.discharge_efficiency
    elseif key.value_key == "CellVoltage"
        return unit.V_cell
    elseif key.value_key == "SOC"
        return unit.SOC
    elseif key.value_key == "ExtractedCharge"
        return unit.extracted_charge
    elseif key.value_key == "Cycles"
        return unit.cycles
    elseif key.value_key == "Temperature"
        return unit.Temp
    end
    throw(KeyError(key.value_key))
end

export calc_V_cell_cc_aging # TODO remove
export f_V_cell # TODO remove
export Battery
