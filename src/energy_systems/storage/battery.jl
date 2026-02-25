using Plots: Plots

"""
Implementation of a battery component holding electric charge.
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
    cell_cutoff_current::Float64
    cycles::Float64
    Temp::Float64
    V_n::Float64
    r_i::Float64
    V_0::Float64
    K::Float64
    A::Float64
    B::Float64
    capacity_cell_Ah::Float64
    m::Float64
    alpha::Float64
    # with aging
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
    V_cell_min::Float64
    V_cell_max::Float64
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
            config["V_n"] = 3.2
            config["r_i"] = 0.00016
            config["V_0"] = 3.36964
            config["K"] = 0.03546
            config["A"] = 0.08165
            config["B"] = 0.1003
            config["capacity_cell_Ah"] = 1090
            config["m"] = 1.0269
            config["alpha"] = -0.01212
            config["k_qn"] = [-1.27571e-7, 1.22095e-11]
            config["k_qT"] = [1.32729e-3, -7.9763e-6]
            config["k_n"] = [9.71249e-6, 7.51635e-4, -8.59363e-5, -2.92489e-4]
            config["k_T"] = [1.05135e-3, 1.83721e-2, -7.72438e-3, -4.31833e-2]
            config["I_ref"] = 100
            config["T_ref"] = 25
        end

        # check that charge_efficiency and discharge_efficiency are set for simplified model
        charge_efficiency = default(config, "charge_efficiency", nothing)
        discharge_efficiency = default(config, "discharge_efficiency", nothing)
        if model_type == "simplified" &&
           (charge_efficiency === nothing || discharge_efficiency === nothing)
            # end of expression
            @error "If model \"simplified\" is used for battery \"$(uac)\" then " *
                   "\"charge_efficiency\" and \"discharge_efficiency\" must be given."
            throw(InputError())
        end

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
                   default(config, "initial_load", 0.0) * config["capacity"], # energy load of the battery in Wh
                   default(config, "charge_efficiency", nothing),
                   default(config, "discharge_efficiency", nothing),
                   default(config, "self_discharge_rate", 0.0), # rate of self-discharge including stand-by losses in %/month (month=30days)
                   default(config, "max_charge_C_rate", 1.0), # maximum continuos charge power in W
                   default(config, "max_discharge_C_rate", 1.0), # maximum continuos discharge power in W
                   default(config, "SOC_min", 0),
                   default(config, "SOC_max", 100),
                   config["V_n_bat"], # nominal voltage of the battery pack
                   default(config, "cell_cutoff_current", 0.003 * default(config, "capacity_cell_Ah", 0.0)),
                   default(config, "cycles", 1.0), # number of cycles
                   default(config, "Temp", 25.0), # cell temperature in °C
                   default(config, "V_n", 3.2), # nominal voltage of the cell in V
                   default(config, "r_i", 0.00016), # internal resistance of the cell in Ohm
                   default(config, "V_0", 3.36964), # battery constant voltage in V
                   default(config, "K", 0.03546), # polarisation voltage in V
                   default(config, "A", 0.08165), # exponential zone amplitude in V
                   default(config, "B", 0.1003), # exponential zone time constant inverse in 1/(Ah)
                   default(config, "capacity_cell_Ah", 1090), # nominal cell capacity in Ah
                   default(config, "m", 1.0269),
                   default(config, "alpha", -0.01212),
                   default(config, "k_qn", [-1.27571e-7, 1.22095e-11]),
                   default(config, "k_qT", [1.32729e-3, -7.9763e-6]),
                   default(config, "k_n", [9.71249e-6, 7.51635e-4, -8.59363e-5, -2.92489e-4]),
                   default(config, "k_T", [1.05135e-3, 1.83721e-2, -7.72438e-3, -4.31833e-2]),
                   default(config, "I_ref", 100.0),
                   default(config, "T_ref", 25.0),
                   0.0, # load_end_of_last_timestep
                   0.0, # losses
                   0.0, # maximum charge energy in current timestep
                   0.0, # maximum discharge energy in current timestep
                   0.0, # battery voltage in V
                   0.0, # battery voltage at the end of the last timestep in V
                   0.0, # minimum cell voltage that defines the battery as empty in V
                   0.0, # maximum cell voltage that defines the battery as full in V
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
        f_n = (1 + unit.k_qn[1] * (n - 1) * n / 2 + unit.k_qn[2] * (n - 1) * n * (2 * n - 1) / 2)
        f_T = (1 + unit.k_qT[1] * (unit.Temp - unit.T_ref) + unit.k_qT[2] * (unit.Temp - unit.T_ref)^2)
        Q = unit.capacity_cell_Ah * 1 * f_n * f_T

        # calculate minimum and maximum cell voltages for SOC_max and SOC_min
        extracted_charge_empty = Q * (1 - (unit.SOC_min) / 100)
        unit.V_cell_min = f_V_cell(extracted_charge_empty, unit.I_ref, unit, n, unit.Temp)
        if unit.V_cell_min < 0
            @error "The battery cell parameters of $(unit.uac) are not correctly defined."
            throw(InputError)
        end

        extracted_charge_full = Q * (1 - unit.SOC_max / 100)
        unit.V_cell_max = f_V_cell(extracted_charge_full, -unit.cell_cutoff_current, unit, n, unit.Temp)

        # get starting point for simulation
        if unit.SOC == 100
            unit.extracted_charge = 0
            unit.V_cell = unit.V_cell_max
        elseif unit.SOC == 0
            unit.extracted_charge = Q
            unit.V_cell = unit.V_cell_min
        else
            unit.V_cell_last = unit.V_n
            # unit.m keeps the SOC consistent to later calculation
            unit.extracted_charge = unit.m * Q * (1 - unit.SOC / 100)
            unit.V_cell = f_V_cell(unit.extracted_charge, unit.I_ref, unit, n, unit.Temp)
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
        set_max_energy!(unit.output_interfaces[unit.m_heat_lt_out], unit.heat_out, nothing,
                        unit.Temp)
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
        # limit Voltage to reasonable range to prevent find_zero to run into an asymptote
        if energy < 0
            V_min = max(unit.V_cell_last * 0.8, unit.V_cell_min)
        elseif energy > 0
            #TODO increased limit on V_max; with very small discharge the Voltage can jump 
            # above V_max and fall through the raster
            V_max = min(unit.V_cell_last * 1.3, unit.V_cell_max)
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
            V_cell = find_zero(V_cell -> V_cell_function(V_cell, 0.0,
                                                         sim_params["time_step_seconds"],
                                                         unit),
                               unit.V_cell_last,
                               Roots.Order1())
            max_energy_cell = 0.0
        # check if the resulting voltage with given energy is between V_min & V_max
        elseif cell_energy < asymp &&
               sign(V_cell_function(V_min, cell_energy, sim_params["time_step_seconds"], unit)) !=
               sign(V_cell_function(V_max, cell_energy, sim_params["time_step_seconds"], unit))
            # end of expression
            V_cell = find_zero(V_cell -> V_cell_function(V_cell, cell_energy,
                                                         sim_params["time_step_seconds"],
                                                         unit),
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
                                                                      sim_params["time_step_seconds"],
                                                                      unit),
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
        max_current = sim_params["wh_to_watts"](max_energy_cell) / V_cell_avg
        max_energy_bat = abs(max_energy_cell) * unit.n_cell_p * unit.n_cell_s
        efficiency = ifelse(max_energy_cell == 0, NaN, 1 - abs(unit.r_i * max_current / V_cell_avg))
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
        max_energy_cell = sim_params["watt_to_wh"](V_cell_avg * max_current)
        max_energy_bat = abs(max_energy_cell) * unit.n_cell_p * unit.n_cell_s
        efficiency = ifelse(max_energy_cell == 0, NaN, 1 - abs(unit.r_i * max_current / V_cell_avg))
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
                # ignore effect of current on capacity for SOC calculation to stabilize the 
                # SOC calculation when current changes a lot between time steps
                f_I = 1
                f_n = (1 + unit.k_qn[1] * (n - 1) * n / 2 + unit.k_qn[2] * (n - 1) * n * (2 * n - 1) / 2)
                f_T = (1 + unit.k_qT[1] * (unit.Temp - unit.T_ref) + unit.k_qT[2] * (unit.Temp - unit.T_ref)^2)
                Q = unit.capacity_cell_Ah * f_I * f_n * f_T

                if unit.extracted_charge > unit.m * Q
                    unit.extracted_charge = unit.m * Q
                elseif unit.extracted_charge < 0
                    unit.extracted_charge = 0
                end
                unit.SOC = ((unit.m * Q) - unit.extracted_charge) / (unit.m * Q) * 100
                unit.load = unit.capacity * unit.SOC / 100
            end
            unit.cycles += abs(unit.extracted_charge - unit.extracted_charge_last) / unit.capacity_cell_Ah / 2
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
    Q_max = unit.capacity_cell_Ah * 1.3
    try
        Q_max = find_zero(Q -> f_V_cell(Q, 
                                        unit.max_discharge_C_rate * unit.capacity_cell_Ah, 
                                        unit, 1, 50) -
                               unit.V_cell_min,
                          unit.capacity_cell_Ah) * 1.001
    catch
    end
    x = 0:(Q_max / 500):Q_max
    y1 = fill(NaN, length(x))
    y2 = fill(NaN, length(x))
    y3 = fill(NaN, length(x))
    for (idx, Q) in enumerate(x)
        y1[idx] = f_V_cell(Q, unit.max_discharge_C_rate * unit.capacity_cell_Ah, unit, 
                           unit.cycles, unit.T_ref)
        if y1[idx] <= unit.V_cell_min
            break
        end
    end
    for (idx, Q) in enumerate(x)
        y2[idx] = f_V_cell(Q, unit.max_discharge_C_rate * 2 * unit.capacity_cell_Ah, unit, 
                           unit.cycles, unit.T_ref)
        if y2[idx] <= unit.V_cell_min
            break
        end
    end
    for (idx, Q) in enumerate(x)
        y3[idx] = f_V_cell(Q, unit.max_discharge_C_rate / 2 * unit.capacity_cell_Ah, unit, 
                           unit.cycles, unit.T_ref)
        if y3[idx] <= unit.V_cell_min
            break
        end
    end
    p = Plots.plot(x, [y1, y2, y3];
                   labels=hcat("$(unit.max_discharge_C_rate)C",
                               "$(unit.max_discharge_C_rate*2)C",
                               "$(unit.max_discharge_C_rate/2)C"),
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
        y1[idx] = f_V_cell(Q, unit.max_discharge_C_rate * unit.capacity_cell_Ah, unit, 
                           unit.cycles, unit.T_ref)
        if y1[idx] <= unit.V_cell_min
            break
        end
    end
    for (idx, Q) in enumerate(x)
        y2[idx] = f_V_cell(Q, unit.max_discharge_C_rate * unit.capacity_cell_Ah, unit, 
                           unit.cycles, 0)
        if y2[idx] <= unit.V_cell_min
            break
        end
    end
    for (idx, Q) in enumerate(x)
        y3[idx] = f_V_cell(Q, unit.max_discharge_C_rate * unit.capacity_cell_Ah, unit, 
                           unit.cycles, 50)
        if y3[idx] <= unit.V_cell_min
            break
        end
    end
    p = Plots.plot(x, [y1, y2, y3];
                   labels=["$(unit.T_ref) °C" "0 °C" "50 °C"],
                   lw=3)
    Plots.plot!(p; 
                title="Battery discharge curve for $(unit.uac) at different Temperatures",
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
        y1[idx] = f_V_cell(Q, unit.max_discharge_C_rate * unit.capacity_cell_Ah, unit, 
                           unit.cycles, unit.T_ref)
        if y1[idx] <= unit.V_cell_min
            break
        end
    end
    for (idx, Q) in enumerate(x)
        y2[idx] = f_V_cell(Q, unit.max_discharge_C_rate * unit.capacity_cell_Ah, unit, 
                           1000, unit.T_ref)
        if y2[idx] <= unit.V_cell_min
            break
        end
    end
    for (idx, Q) in enumerate(x)
        y3[idx] = f_V_cell(Q, unit.max_discharge_C_rate * unit.capacity_cell_Ah, unit, 
                           3000, unit.T_ref)
        if y3[idx] <= unit.V_cell_min
            break
        end
    end
    p = Plots.plot(x, [y1, y2, y3];
                   labels=["$(unit.cycles) cycles" "1000 cycles" "3000 cycles"],
                   lw=3)
    Plots.plot!(p;
                title="Battery discharge curve for $(unit.uac) after different number of cycles",
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

# calculate the voltage level, after charging with a constant current defined from charging 
# power and voltage level before charging, positive energy = discharging, negative e
# nergy = charging; changing voltage levels during charging are ignored
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

    if unit.m * Q < extracted_charge
        V_cell = unit.V_cell_min - 0.1
    elseif extracted_charge < 0
        V_cell = unit.V_cell_max + 0.1
    else
        V_0 = unit.V_0 * (1 + unit.k_n[1] * (n - 1)) * (1 + unit.k_T[1] * (T - unit.T_ref))
        r_i = unit.r_i * (1 + unit.k_n[2] * (n - 1)) * (1 + unit.k_T[2] * (T - unit.T_ref))
        K = unit.K * (1 + unit.k_n[3] * (n - 1)) * (1 + unit.k_T[3] * (T - unit.T_ref))
        A = unit.A * (1 + unit.k_n[4] * (n - 1)) * (1 + unit.k_T[4] * (T - unit.T_ref))

        V_cell = V_0 -
                 K * (unit.m * Q) / (unit.m * Q - extracted_charge) +
                 A * exp(-unit.B * extracted_charge) -
                 r_i * current
    end
    return V_cell
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

export Battery
