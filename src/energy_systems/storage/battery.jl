using Statistics

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

    medium::Symbol
    model_type::String

    capacity::Float64
    load::Float64
    
    charge_efficiency::Floathing
    discharge_efficiency::Floathing
    self_discharge_rate::Float64
    max_charge_power::Float64
    max_discharge_power::Float64
    SOC_min::Float64
    SOC_max::Float64
    V_n_bat::Float64
    # detailed
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
    # with aging
    n::Float64
    T::Float64
    k_qn::Array{Float64, 1}
    k_qT::Array{Float64, 1}
    k_n::Array{Float64, 1}
    k_T::Array{Float64, 1}
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
    r_sd::Float64
    n_cell_p::Float64
    n_cell_s::Float64
    charge_sum::Float64
    charge_sum_last::Float64
    SOC::Float64

    process_done::Bool
    load_done::Bool

    time_array_1::Any
    time_array_2::Any
    time_array_3::Any
    max_capacity_last::Float64

    function Battery(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        medium = Symbol(default(config, "medium", "m_e_ac_230v"))
        register_media([medium])
        
        # get model type from input file
        model_type = default(config, "model_type", "simplified")
        model_type_allowed_values = ["simplified", "no_aging", "with_aging"]
        if !(model_type in model_type_allowed_values)
            @error "Undefined model type \"$(model_type)\" of battery \"$(uac)\". Has to be one of: $(model_type_allowed_values)."
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
                   InterfaceMap(medium => nothing), # input_interfaces
                   InterfaceMap(medium => nothing), # output_interfaces
                   medium,
                   model_type, # model_type
                   config["capacity"], # capacity of the battery in Wh
                   config["load"], # energy load of the battery in Wh
                   default(config, "charge_efficiency", nothing),
                   default(config, "discharge_efficiency", nothing),
                   default(config, "self_discharge_rate", 0.0), # rate of self-discharge inlcuding stand-by losses in %/month (month=30days)
                   default(config, "max_charge_power", config["capacity"]), # maximum continuos charge power in W
                   default(config, "max_discharge_power", config["capacity"]), # maximum continuos discharge power in W
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
                   default(config, "m", 0.0),
                   default(config, "alpha", 0.0),
                   default(config, "cycles", 1.0), # number of cycles
                   default(config, "T_bat", 25.0), # cell temperature in °C
                   default(config, "k_qn", [0.0,0.0]),
                   default(config, "k_qT", [0.0,0.0]),
                   default(config, "k_n", [0.0,0.0,0.0,0.0]),
                   default(config, "k_T", [0.0,0.0,0.0,0.0]),
                   default(config, "I_ref", 0.0),
                   default(config, "T_ref", 25.0),
                   0.0, # load_end_of_last_timestep
                   0.0, # losses
                   0.0, # maximum charge energy in current timestep
                   0.0, # maximum discharge energy in current timestep
                   0.0, # battery voltage in V
                   0.0, # battery voltage at the end of the last timestep in V
                   0.0, # battery voltage after charging
                   0.0, # battery volatge after discharging
                   0.0, # self-discharge resistance used to calculate the self discharge in Ohm
                   0.0, # number of cells connected in parallel 
                   0.0, # number of cells connected in series
                   0.0, # sum of all added and removed charge to the cell in Ah
                   0.0, # sum of all added and removed charge to the cell in the last timestep in Ah
                   0.0, # SOC in the current time step
                   false, # process_done, bool indicating if the process step has already been performed in the current time step
                   false, # load_done, bool indicating if the load step has already been performed in the current time step
                   [],
                   [],
                   [],
                   default(config, "capacity_cell_Ah", 0.0)) 
    end
end

function initialise!(unit::Battery, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.input_interfaces[unit.medium],
                          unload_storages(unit.controller, unit.medium))
    set_storage_transfer!(unit.output_interfaces[unit.medium],
                          load_storages(unit.controller, unit.medium))

    unit.load_end_of_last_timestep = copy(unit.load)
    unit.self_discharge_rate *= sim_params["time_step_seconds"] / (30 * 24 * 3600)

    if unit.model_type != "simplified"
        # estimate cells in parallel and series in the battery
        unit.n_cell_s = unit.V_n_bat / unit.V_n
        unit.n_cell_p = (unit.capacity / unit.V_n) / unit.capacity_cell_Ah

        f_I = 1
        f_n = (1 + unit.k_qn[1]*(unit.n-1)*unit.n / 2 + unit.k_qn[2]*(unit.n-1)*unit.n*(2*unit.n-1) / 2) 
        f_T = (1 + unit.k_qT[1]*(unit.T-unit.T_ref) + unit.k_qT[2]*(unit.T-unit.T_ref)^2)
        Q = unit.capacity_cell_Ah * f_I * f_n * f_T
        # calculate minimum and maximum cell voltages for SOC_max and SOC_min
        if unit.SOC_min > 0
            unit.V_cell_last = unit.V_n
            unit.charge_sum = unit.m*Q * (1 - unit.SOC_min/100)
            unit.V_cell_min = find_zero(V_cell -> V_cell_function(V_cell, 0, unit, sim_params),
                                        unit.V_n,
                                        Order1())
        end

        unit.V_cell_last = unit.V_n
        unit.charge_sum = unit.m*Q * (1 - unit.SOC_max/100)
        unit.V_cell_max = find_zero(V_cell -> V_cell_function(V_cell, 0, unit, sim_params),
                                    unit.V_n,
                                    Order1())

        unit.SOC = (unit.load / unit.capacity) * 100
        if unit.SOC < unit.SOC_min || unit.SOC > unit.SOC_max
            @info "Starting load is not in the defined SOC range. The load is set to the closest valid value"
            if abs(unit.SOC_min - unit.SOC) < abs(unit.SOC_max - unit.SOC)
                unit.SOC = unit.SOC_min
                unit.load = unit.capacity * unit.SOC_min/100
            else
                unit.SOC = unit.SOC_max
                unit.load = unit.capacity * unit.SOC_max/100
            end
        end
        if unit.SOC == 100 || false
            unit.charge_sum = 0
            unit.V_cell = unit.V_cell_max
        elseif unit.SOC == 0 || false
            unit.charge_sum = unit.m*Q
            unit.V_cell = unit.V_cell_min
        else
            unit.V_cell_last = unit.V_n
            unit.charge_sum = unit.m*Q * (1 - unit.SOC/100)
            unit.V_cell = find_zero(V_cell -> V_cell_function(V_cell, 0, unit, sim_params),
                                    unit.V_n,
                                    Order1())
        end
        unit.V_cell_last = copy(unit.V_cell)
        unit.charge_sum_last = copy(unit.charge_sum)
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

    if discharge_is_allowed(unit.controller, sim_params) && unit.SOC > unit.SOC_min && unit.V_cell > unit.V_cell_min
        unit.discharge_efficiency, 
        unit.V_cell_discharge, 
        unit.max_discharge_energy = calc_efficiency(sim_params["watt_to_wh"](unit.max_discharge_power), unit, sim_params)
    else
        unit.max_discharge_energy = 0.0
        unit.discharge_efficiency = NaN
    end
    set_max_energy!(unit.output_interfaces[unit.medium], unit.max_discharge_energy)

    if charge_is_allowed(unit.controller, sim_params) && unit.SOC < unit.SOC_max
        unit.charge_efficiency, 
        unit.V_cell_charge, 
        charge_energy_bat = calc_efficiency(sim_params["watt_to_wh"](-unit.max_charge_power), unit, sim_params)
        unit.max_charge_energy = charge_energy_bat
    else
        unit.max_charge_energy = 0.0
        unit.charge_efficiency = NaN
    end
    set_max_energy!(unit.input_interfaces[unit.medium], unit.max_charge_energy)
end

function calc_efficiency(energy::Number, unit::Battery, sim_params::Dict{String,Any})
    if unit.model_type == "simplified"
        if energy >= 0
            return unit.discharge_efficiency, 0, energy
        else
            return unit.charge_efficiency, 0, energy
        end
    else
        cell_energy = energy / unit.n_cell_p
        # calculate V_cell
        V_min = unit.V_cell_min
        V_max = unit.V_cell_max
        # if energy < 0
        #     V_min = unit.V_cell_last - 0.5
        # elseif energy > 0
        #     V_max = unit.V_cell_last + 0.5
        # end
        if sign(V_cell_function(V_min, cell_energy, unit, sim_params)) != sign(V_cell_function(V_max, cell_energy, unit, sim_params))
            #TODO fails if deeply discharged and little energy is added
            V_cell = find_zero(V_cell -> V_cell_function(V_cell, cell_energy, unit, sim_params),
                               (V_min, V_max),
                               Roots.Brent())
            max_energy_cell = cell_energy

        elseif energy > 0
            as = find_zero(e_cell -> denom(unit.V_cell_min, e_cell, unit, sim_params), 
                           unit.capacity_cell_Ah * unit.V_n / 2, 
                           Roots.Order1())
            ul = min(as*0.999, cell_energy)
            try max_energy_cell = find_zero(e_cell -> V_cell_function(unit.V_cell_min, e_cell, 
                                                                  unit, sim_params), 
                                        (0.01, ul), 
                                        Roots.Brent())
            catch
                @infiltrate
            end
            V_cell = unit.V_cell_min

        elseif energy < 0
            # max_energy_cell = find_zero((e_cell -> V_cell_function(unit.V_cell_max, 
            #                                                             e_cell, 
            #                                                             unit, sim_params),
            #                              e_cell -> energy_cell_derivative(unit.V_cell_max,
            #                                                               e_cell, 
            #                                                               unit, 
            #                                                               sim_params)),
            #                             -unit.capacity * (1 - unit.SOC/100) / unit.n_cell_p,
            #                             Roots.Newton()) 
            try max_energy_cell = find_zero(e_cell -> V_cell_function(unit.V_cell_max, e_cell, 
                                                                      unit, sim_params), 
                                        (cell_energy, -0.001), 
                                        Roots.Brent())
            catch
                try max_energy_cell = find_zero(e_cell -> V_cell_function(unit.V_cell_max, e_cell, 
                                                                      unit, sim_params), 
                                        cell_energy, 
                                        Roots.Order1())
                catch
                    @infiltrate
                end
            end
            V_cell = unit.V_cell_max
        end

        V_cell_avg = (unit.V_cell_last + V_cell) / 2
        charge = max_energy_cell / V_cell_avg
        internal_losses = unit.r_i * charge^2        
        max_energy_bat = abs(max_energy_cell) * unit.n_cell_p
        efficiency = ifelse(max_energy_cell == 0, NaN, 1- (internal_losses / abs(max_energy_cell)))
        if V_cell > unit.V_cell_max + 0.01  || V_cell < unit.V_cell_min
            @infiltrate
        end
        return efficiency, V_cell, max_energy_bat
    end
end

# function to find the voltage level, after charging with a constant power with aging 
# processes considered; positive energy = discharging, negative energy = charging
function V_cell_function(V_cell::Number, energy::Number, unit::Battery, 
                         sim_params::Dict{String,Any}, n::Number=1, T::Number=25)
    charge = 2*energy / (unit.V_cell_last + V_cell)
    current = charge / (sim_params["time_step_seconds"] / 3600)
    # functions describing the dependence of the parameters on cycles n and Temperature T
    if energy == 0
        f_I = 1
    else        
        f_I = (abs(current) / unit.I_ref)^unit.alpha
    end
    f_n = (1 + unit.k_qn[1]*(n-1)*n / 2 + unit.k_qn[2]*(n-1)*n*(2*n-1) / 2) 
    f_T = (1 + unit.k_qT[1]*(T-unit.T_ref) + unit.k_qT[2]*(T-unit.T_ref)^2)
    Q = unit.capacity_cell_Ah * f_I * f_n * f_T                     # wird kleiner mit größerer energy  

    if (unit.m * Q - (unit.charge_sum + charge)) < 0
        y = 1
    else
        V_0 = unit.V_0 * (1 + unit.k_n[1]*(n-1)) * (1 + unit.k_T[1]*(T-unit.T_ref))
        r_i = unit.r_i * (1 + unit.k_n[2]*(n-1)) * (1 + unit.k_T[2]*(T-unit.T_ref))
        K = unit.K * (1 + unit.k_n[3]*(n-1)) * (1 + unit.k_T[3]*(T-unit.T_ref))
        A = unit.A * (1 + unit.k_n[4]*(n-1)) * (1 + unit.k_T[4]*(T-unit.T_ref))

        y = V_0 - 
            K * (unit.m * Q) / (unit.m * Q - (unit.charge_sum + charge)) +  # verhält sich seltsam (vorzeichen wechsel wenn charge zu groß)
            A * exp(-unit.B * (unit.charge_sum + charge)) -                 # irrelevant in den meisten Fällen
            r_i * current -                                                                             # wird größer mit größerer energy
            V_cell
    end
    return y
end

function denom(V_cell, energy, unit, sim_params, n::Number=1, T::Number=25)
    f_n = (1 + unit.k_qn[1]*(n-1)*n / 2 + unit.k_qn[2]*(n-1)*n*(2*n-1) / 2) 
    f_T = (1 + unit.k_qT[1]*(T-unit.T_ref) + unit.k_qT[2]*(T-unit.T_ref)^2)
    f_I = (abs(2*energy) / (unit.I_ref * (unit.V_cell_last + V_cell) * (sim_params["time_step_seconds"] / 3600)))^unit.alpha

    unit.m * unit.capacity_cell_Ah * f_I * f_n * f_T  - 
    (unit.charge_sum + 2*energy / (unit.V_cell_last + V_cell))
end

function denom_derivative(V_cell, energy, unit, sim_params, n::Number=1, T::Number=25)
    f_n = (1 + unit.k_qn[1]*(n-1)*n / 2 + unit.k_qn[2]*(n-1)*n*(2*n-1) / 2) 
    f_T = (1 + unit.k_qT[1]*(T-unit.T_ref) + unit.k_qT[2]*(T-unit.T_ref)^2)
    df_I = unit.alpha * (abs(2*energy) / (unit.I_ref * (unit.V_cell_last + V_cell) * (sim_params["time_step_seconds"] / 3600)))^(unit.alpha - 1) * 
           (2.0 / (unit.I_ref * (sim_params["time_step_seconds"] / 3600.0))) * 
           sign(energy) * (unit.V_cell_last + V_cell)^(-unit.alpha)

    unit.m * unit.capacity_cell_Ah * df_I * f_n * f_T  - 2 / (unit.V_cell_last + V_cell)
end

function V_cell_derivative(V_cell::Number, energy::Number, unit::Battery, 
                           sim_params::Dict{String,Any}, n::Number=1, T::Number=25)
    charge = 2*energy / (unit.V_cell_last+V_cell)
    current = charge / (sim_params["time_step_seconds"] / 3600) 
    if energy == 0
        return -1
    elseif V_cell < 0
        return 0
    else      
        f_I = (abs(current) / unit.I_ref)^unit.alpha
        df_I = -unit.alpha * (abs(current) / unit.I_ref)^unit.alpha * (unit.V_cell_last + V_cell)^(-unit.alpha - 1)
    end
    f_n = (1 + unit.k_qn[1]*(n-1)*n / 2 + unit.k_qn[2]*(n-1)*n*(2*n-1) / 2) 
    f_T = (1 + unit.k_qT[1]*(T-unit.T_ref) + unit.k_qT[2]*(T-unit.T_ref)^2)
    Q = unit.capacity_cell_Ah * f_I * f_n * f_T
    dQ = unit.capacity_cell_Ah * df_I * f_n * f_T

    r_i = unit.r_i * (1 + unit.k_n[2]*(n-1)) * (1 + unit.k_T[2]*(T-unit.T_ref))
    K = unit.K * (1 + unit.k_n[3]*(n-1)) * (1 + unit.k_T[3]*(T-unit.T_ref))
    A = unit.A * (1 + unit.k_n[4]*(n-1)) * (1 + unit.k_T[4]*(T-unit.T_ref))

    C = unit.charge_sum + charge

    K * unit.m * (dQ*C + Q*(2*energy) / (unit.V_cell_last+V_cell)^2) / (unit.m*Q - C)^2 +
    A * (2*energy*unit.B) / (unit.V_cell_last+V_cell)^2 * exp(-unit.B * C) +
    r_i * (2*energy) / ((unit.V_cell_last+V_cell)^2 * (sim_params["time_step_seconds"] / 3600)) -
    1
end

function energy_cell_derivative(V_cell::Number, energy::Number, unit::Battery, 
                           sim_params::Dict{String,Any}, n::Number=1, T::Number=25)
    charge = 2*energy / (unit.V_cell_last+V_cell)
    current = charge / (sim_params["time_step_seconds"] / 3600) 
    if energy == 0
        f_I = 1
        df_I = 1
    else      
        f_I = (abs(current) / unit.I_ref)^unit.alpha
        df_I = unit.alpha * (abs(current) / unit.I_ref)^(unit.alpha - 1) * 
               (2.0 / (unit.I_ref * (sim_params["time_step_seconds"] / 3600.0))) * 
               sign(energy) * (unit.V_cell_last + V_cell)^(-unit.alpha)
    end
    f_n = (1 + unit.k_qn[1]*(n-1)*n / 2 + unit.k_qn[2]*(n-1)*n*(2*n-1) / 2) 
    f_T = (1 + unit.k_qT[1]*(T-unit.T_ref) + unit.k_qT[2]*(T-unit.T_ref)^2)
    Q = unit.capacity_cell_Ah * f_I * f_n * f_T
    dQ = unit.capacity_cell_Ah * df_I * f_n * f_T

    r_i = unit.r_i * (1 + unit.k_n[2]*(n-1)) * (1 + unit.k_T[2]*(T-unit.T_ref))
    K = unit.K * (1 + unit.k_n[3]*(n-1)) * (1 + unit.k_T[3]*(T-unit.T_ref))
    A = unit.A * (1 + unit.k_n[4]*(n-1)) * (1 + unit.k_T[4]*(T-unit.T_ref))

    C = unit.charge_sum + charge

    K * unit.m * (dQ*C - Q * 2/(unit.V_cell_last+V_cell)) / (unit.m*Q - C)^2 -
    A * (2*unit.B) / (unit.V_cell_last+V_cell) * exp(-unit.B * C) -
    r_i * 2/((unit.V_cell_last+V_cell) * (sim_params["time_step_seconds"] / 3600))
end

function process(unit::Battery, sim_params::Dict{String,Any})
    if unit.max_discharge_energy >= sim_params["epsilon"]
        outface = unit.output_interfaces[unit.medium]
        exchanges = balance_on(outface, outface.target)
        energy_demand = balance(exchanges) + energy_potential(exchanges)

        if energy_demand >= 0.0 # process is only concerned with moving energy to the target
            set_max_energy!(unit.output_interfaces[unit.medium], 0.0)
            unit.discharge_efficiency = NaN
        else
            #TODO discharge_energy_bat is wrong after recalcuting discharge_efficiency leads to load < 0
            # discharge_energy_bat = min((unit.SOC - unit.SOC_min)/100 * unit.capacity, abs(energy_demand) / unit.discharge_efficiency)
            # if discharge_energy_bat != unit.max_discharge_energy / unit.discharge_efficiency
            #     unit.discharge_efficiency, unit.V_cell, discharge_energy_bat_temp = calc_efficiency(discharge_energy_bat, unit, sim_params)

            if abs(energy_demand) < unit.max_discharge_energy
                unit.discharge_efficiency, unit.V_cell, discharge_energy_bat = calc_efficiency(abs(energy_demand), unit, sim_params)
            else
                discharge_energy_bat = unit.max_discharge_energy
                unit.V_cell = unit.V_cell_discharge
            end       
            unit.losses += discharge_energy_bat * (1 - unit.discharge_efficiency)
            unit.charge_sum += discharge_energy_bat / unit.n_cell_p * 2 / (unit.V_cell_last + unit.V_cell)

            unit.load -= discharge_energy_bat
            add!(outface, discharge_energy_bat)
        end
    else
        set_max_energy!(unit.output_interfaces[unit.medium], 0.0)
    end
    handle_component_update!(unit, "process", sim_params)

    # if sim_params["time"] > 670*900
    #     unit.time_array_1 = unit.time_array_1[2:end]
    #     unit.time_array_2 = unit.time_array_2[2:end]
    #     unit.time_array_3 = unit.time_array_3[2:end]
    #     @info unit.time_array_1
    #     @info unit.time_array_2
    #     @info unit.time_array_3
    #     @info "Mean:$(round(mean(unit.time_array_1), digits=3)) us; std: $(round(std(unit.time_array_1), digits=3)) us; max:$(round(maximum(unit.time_array_1), digits=3)) us"
    #     @info "Mean:$(round(mean(unit.time_array_2), digits=3)) us; std: $(round(std(unit.time_array_2), digits=3)) us; max:$(round(maximum(unit.time_array_2), digits=3)) us"
    #     @info "Mean:$(round(mean(unit.time_array_3), digits=3)) us; std: $(round(std(unit.time_array_3), digits=3)) us; max:$(round(maximum(unit.time_array_3), digits=3)) us"
    # end
end

function load(unit::Battery, sim_params::Dict{String,Any})
    if unit.max_charge_energy >= sim_params["epsilon"]
        inface = unit.input_interfaces[unit.medium]
        exchanges = balance_on(inface, inface.source)
        energy_available = balance(exchanges) + energy_potential(exchanges)

        if energy_available <= 0.0 # load is only concerned with receiving energy from the source
            set_max_energy!(unit.input_interfaces[unit.medium], 0.0)
            unit.charge_efficiency = NaN
        else
            if energy_available < unit.max_charge_energy
                unit.charge_efficiency, unit.V_cell, charge_energy_bat = calc_efficiency(-energy_available, unit, sim_params)
            else
                charge_energy_bat = unit.max_charge_energy
                unit.V_cell = unit.V_cell_charge
            end 

            unit.losses += charge_energy_bat  * (1.0 / unit.charge_efficiency - 1)
            unit.charge_sum -= charge_energy_bat / unit.n_cell_p * 2 / (unit.V_cell_last + unit.V_cell)

            unit.load += charge_energy_bat
            sub!(inface, charge_energy_bat)
        end
    else
        set_max_energy!(unit.input_interfaces[unit.medium], 0.0)
    end
    handle_component_update!(unit, "load", sim_params)
end

function handle_component_update!(unit::Battery, step::String, sim_params)
    if step == "process"
        unit.process_done = true
    elseif step == "load"
        unit.load_done = true
    end
    if unit.process_done && unit.load_done
        # update component
        self_loss = 0.0
        if unit.load_end_of_last_timestep == unit.load || unit.charge_sum_last == unit.charge_sum && unit.SOC > unit.SOC_min # TODO SOC doesnt change with self discharge for detailted model
            self_loss = unit.capacity * unit.self_discharge_rate # TODO fix loss calculation with respect to SOC and charge
        end

        if unit.model_type == "simplified"
            unit.losses += self_loss
            unit.load -= unit.losses
            unit.SOC = unit.load / unit.capacity * 100
        else
            # calculated current cell capacity
            current = (unit.charge_sum - unit.charge_sum_last) * 2 / (unit.V_cell_last + unit.V_cell)
            if unit.V_cell == unit.V_cell_max
                @infiltrate
            end
            # CV charging current cutoff is set at 3% of the nominal cell capacity
            # if unit.V_cell == unit.V_cell_max && abs(current) < (0.03 * unit.capacity_cell_Ah)# && unit.SOC < 100
            # if current < 0 && abs(current) < (0.03 * unit.max_charge_power/unit.V_n) 
            # if current < 0 && abs(current) < (0.03 * unit.capacity_cell_Ah) 
            # if sim_params["time"] >= 41400
            #     @infiltrate
            # end
            if current == 0
                #TODO gives weird behaviour because SOC and charge_sum are recalculated differently especially if SOC was set to SOC_max in time step before
                # f_I = (abs(unit.current_last) / unit.I_ref)^unit.alpha
                # f_n = (1 + unit.k_qn[1]*(unit.n-1)*unit.n / 2 + unit.k_qn[2]*(unit.n-1)*unit.n*(2*unit.n-1) / 2) 
                # f_T = (1 + unit.k_qT[1]*(unit.T-unit.T_ref) + unit.k_qT[2]*(unit.T-unit.T_ref)^2)
                # Q = unit.capacity_cell_Ah * f_I * f_n * f_T

                unit.SOC = max(0, (unit.load - self_loss) / unit.capacity * 100)
                unit.charge_sum = unit.max_capacity_last * (1 - unit.SOC/100) # TODO überprüfen ob das sinnvoll ist
                # TODO for testing
                # _, unit.V_cell, self_loss = calc_efficiency(self_loss, unit, sim_params)
                # unit.charge_sum += self_loss  / unit.n_cell_p * 2 / (unit.V_cell_last + unit.V_cell) 
        
            elseif unit.V_cell == unit.V_cell_max && abs(current) < (0.03 * unit.capacity_cell_Ah) && unit.SOC >= unit.SOC_max*0.99
                unit.SOC = unit.SOC_max

            elseif current != 0
                f_I = (abs(current) / unit.I_ref)^unit.alpha
                f_n = (1 + unit.k_qn[1]*(unit.n-1)*unit.n / 2 + unit.k_qn[2]*(unit.n-1)*unit.n*(2*unit.n-1) / 2) 
                f_T = (1 + unit.k_qT[1]*(unit.T-unit.T_ref) + unit.k_qT[2]*(unit.T-unit.T_ref)^2)
                Q = unit.capacity_cell_Ah * f_I * f_n * f_T
                
                if unit.charge_sum > unit.m * Q
                    unit.charge_sum = unit.m * Q
                elseif unit.charge_sum < 0
                    unit.charge_sum = 0
                end
                unit.SOC = ((unit.m * Q) - unit.charge_sum) / (unit.m * Q) * 100
                unit.max_capacity_last = unit.m*Q
                unit.load = unit.capacity * unit.SOC / 100 
            end
            
            unit.V_cell_last = copy(unit.V_cell)
        end

        unit.load_end_of_last_timestep = copy(unit.load)
        unit.charge_sum_last = copy(unit.charge_sum)
        # reset
        unit.process_done = false
        unit.load_done = false
    end
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
                          r_i::Union{Number, Nothing}=nothing, V_6=0.0, Q_6=0.0,
                          I_n=0.0, n_1=0.0, Q_full_n_1=0.0,
                          n_2::Union{Number, Nothing}=nothing, V_full_n_2=0.0, V_nom_n_2=0.0, Q_nom_n_2=0.0, Q_full_n_2=0.0,
                          I_T=0.0, T_1=0.0, Q_full_T_1=0.0,
                          T_2::Union{Number, Nothing}=nothing, V_full_T_2=0.0, V_nom_T_2=0.0, Q_nom_T_2=0.0, Q_full_T_2=0.0,
                          T_ref=25)

    # basic values new
    B = -1/Q_2 * log((V_full - V_3)/(V_full - V_2) - 1) 
    A = (V_full - V_3)/(1 - exp(-B*Q_3))

    AB_4 = V_full - V_4 - A*(1 - exp(-B*Q_4)) # exp(-B*Q_4) -> 0
    AB_5 = V_full - V_cut - A*(1 - exp(-B*Q_full_1)) # exp(-B*Q_full) -> 0

    m = (1-AB_5/AB_4) / (1-(AB_5*Q_4)/(AB_4*Q_full_1)) * (Q_4/Q_full_1)
    K = AB_5 * ((m * Q_full_1/Q_full_1) - 1)
    if isnothing(r_i)
        r_i = 1/(I_1-I_2) * (V_6 + K*(m*Q_full_2/(m*Q_full_2-Q_6)) - A*exp(-B*Q_6) -
                             V_3 - K*(m*Q_full_1/(m*Q_full_1-Q_3)) + A*exp(-B*Q_3))
    end 
    V_0 = V_full + K + r_i*I_1 - A
    alpha = log(Q_full_2 / Q_full_1) / log(I_2 / I_1)
    I_ref = I_1
    Q_ref = Q_full_1

    # basic values old
    # A = V_full - V_exp
    # B = 3/Q_exp
    # K = ((V_full - V_nom + A*(exp(-B*Q_nom) - 1))*(Q_full_1 - Q_nom))/(Q_nom)
    # V_0 = V_full + K + r_i*I_1 - A

    # for aging if the needed values are provided
    k_qn = [0.0 ,0.0]
    k_n = [0.0, 0.0, 0.0, 0.0]
    if !isnothing(n_2)
        N = [(n_1-1)*n_1/2   (n_1-1)*n_1*(2*n_1-1)/2
             (n_2-1)*n_2/2   (n_2-1)*n_2*(2*n_2-1)/2]
        q = [Q_full_n_1/(Q_ref*(I_n/I_ref)^alpha) - 1
             Q_full_n_2/(Q_ref*(I_n/I_ref)^alpha) - 1]
        # k_qn[1] = (q[2] - q[1]*N[2,2]/N[1,2]) / (N[2,1] - N[1,1]*N[2,2]/N[1,2])
        # k_qn[2] = (q[1] - k_qn[1]*N[1,1]) / N[1,2]
        k_qn = N \ q

        Q_n = Q_ref * (I_n / I_ref)^alpha * 
              (1 + k_qn[1]*(n_2-1)*n_2 / 2 + k_qn[2]*(n_2-1)*n_2*(2*n_2-1) / 2) 
        b_n = [V_full - (V_0 - r_i*I_1 - K + A)
               V_full_n_2 - (V_0 - r_i*I_n - K + A)
               V_nom_n_2 - (V_0 - r_i*I_n - K*(m*Q_n/(m*Q_n-Q_nom_n_2)) + A*exp(-B*Q_nom_n_2))
               V_cut - (V_0 - r_i*I_n - K*(m*Q_n/(m*Q_n-Q_full_n_2)) + A*exp(-B*Q_full_n_2))]
        G_n = [V_0          -r_i*I_1            -K                                      A
               V_0*(n_2-1)  -r_i*I_n*(n_2-1)    -K*(n_2-1)                              A*(n_2-1)
               V_0*(n_2-1)  -r_i*I_n*(n_2-1)    -K*(n_2-1)*(m*Q_n/(m*Q_n-Q_nom_n_2))    A*(n_2-1)*exp(-B*Q_nom_n_2)
               V_0*(n_2-1)  -r_i*I_n*(n_2-1)    -K*(n_2-1)*(m*Q_n/(m*Q_n-Q_full_n_2))   A*(n_2-1)*exp(-B*Q_full_n_2)] 
        k_n = G_n \ b_n
    end
    
    k_qT = [0.0 ,0.0]
    k_T = [0.0, 0.0, 0.0, 0.0]
     if !isnothing(T_2)
        T = [T_1-T_ref   (T_1-T_ref)^2
             T_2-T_ref   (T_2-T_ref)^2]
        q = [Q_full_T_1/(Q_ref*(I_T/I_ref)^alpha) - 1
             Q_full_T_2/(Q_ref*(I_T/I_ref)^alpha) - 1]
        # k_qT[1] = (q[2] - q[1]*T[2,2]/T[1,2]) / (T[2,1] - T[1,1]*T[2,2]/T[1,2])
        # k_qT[2] = (q[1] - k_qn[1]*T[1,1]) / T[1,2]
        k_qT = T \ q
     
        Q_T = Q_ref * (I_T / I_ref)^alpha * 
              (1 + k_qT[1]*(T_2-T_ref) + k_qT[2]*(T_2-T_ref)^2) 
        b_T = [V_full - (V_0 - r_i*I_1 - K + A)
               V_full_T_2 - (V_0 - r_i*I_T - K + A)
               V_nom_T_2 - (V_0 - r_i*I_T - K*(m*Q_T/(m*Q_T-Q_nom_T_2)) + A*exp(-B*Q_nom_T_2))
               V_cut - (V_0 - r_i*I_T - K*(m*Q_T/(m*Q_T-Q_full_T_2)) + A*exp(-B*Q_full_T_2))]
        G_T = [V_0              -r_i*I_1                -K                                          A
               V_0*(T_2-T_ref)  -r_i*I_T*(T_2-T_ref)    -K*(T_2-T_ref)                              A*(T_2-T_ref)
               V_0*(T_2-T_ref)  -r_i*I_T*(T_2-T_ref)    -K*(T_2-T_ref)*(m*Q_T/(m*Q_T-Q_nom_T_2))    A*(T_2-T_ref)*exp(-B*Q_nom_T_2)
               V_0*(T_2-T_ref)  -r_i*I_T*(T_2-T_ref)    -K*(T_2-T_ref)*(m*Q_T/(m*Q_T-Q_full_T_2))   A*(T_2-T_ref)*exp(-B*Q_full_T_2)] 
        k_T = G_T \ b_T
    end

    return V_0, K, A, B, r_i, m, alpha, k_qn, k_qT, k_n, k_T, I_ref, T_ref
end


# calculate the voltage level, after charging with a constant current defined from charging 
# power and voltage level before charging, positive energy = discharging, negative energy = charging
# changing voltage levels during charging are ignored
function calc_V_cell_cc_aging(charge::Number, I::Number, n ,T, unit::Battery)
    # functions describing the dependence of the parameters on cycles n and Temperature T
    Q = unit.capacity_cell_Ah * (abs(I) / unit.I_ref)^unit.alpha * 
        (1 + unit.k_qn[1]*(n-1)*n / 2 + unit.k_qn[2]*(n-1)*n*(2*n-1) / 2) *
        (1 + unit.k_qT[1]*(T-unit.T_ref) + unit.k_qT[2]*(T-unit.T_ref)^2) 
    V_0 = unit.V_0 * (1 + unit.k_n[1]*(n-1)) * (1 + unit.k_T[1]*(T-unit.T_ref))
    r_i = unit.r_i * (1 + unit.k_n[2]*(n-1)) * (1 + unit.k_T[2]*(T-unit.T_ref))
    K   = unit.K   * (1 + unit.k_n[3]*(n-1)) * (1 + unit.k_T[3]*(T-unit.T_ref))
    A   = unit.A   * (1 + unit.k_n[4]*(n-1)) * (1 + unit.k_T[4]*(T-unit.T_ref))

    V_0 - 
    K * (unit.m * Q) / (unit.m * Q - (unit.charge_sum + charge)) +
    A * exp(-unit.B * (unit.charge_sum + charge)) -
    r_i * I
end

# calculate the voltage level, after charging with a constant current defined from charging 
# power and voltage level before charging, positive energy = discharging, negative energy = charging
# changing voltage levels during charging are ignored
function calc_V_cell_cc_sum_aging(sum_charge::Number, I::Number, n, T, unit::Battery)
    # functions describing the dependence of the parameters on cycles n and Temperature T
    Q = unit.capacity_cell_Ah * (abs(I) / unit.I_ref)^unit.alpha * 
        (1 + unit.k_qn[1]*(n-1)*n / 2 + unit.k_qn[2]*(n-1)*n*(2*n-1) / 2) *
        (1 + unit.k_qT[1]*(T-unit.T_ref) + unit.k_qT[2]*(T-unit.T_ref)^2) 

    V_0 = unit.V_0 * (1 + unit.k_n[1]*(n-1)) * (1 + unit.k_T[1]*(T-unit.T_ref))
    r_i = unit.r_i * (1 + unit.k_n[2]*(n-1)) * (1 + unit.k_T[2]*(T-unit.T_ref))
    K   = unit.K   * (1 + unit.k_n[3]*(n-1)) * (1 + unit.k_T[3]*(T-unit.T_ref))
    A   = unit.A   * (1 + unit.k_n[4]*(n-1)) * (1 + unit.k_T[4]*(T-unit.T_ref))

    V_0 - 
    K * (unit.m * Q) / (unit.m * Q - sum_charge) +
    A * exp(-unit.B * sum_charge) -
    r_i * I
end

function output_values(unit::Battery)::Vector{String}
    return [string(unit.medium) * " IN",
            string(unit.medium) * " OUT",
            "Load",
            "Load%",
            "Capacity",
            "LossesGains",
            "charge_efficiency",
            "discharge_efficiency",
            "CellVoltage",
            "SOC",
            "Charge"]
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
    elseif key.value_key == "Charge"
        return unit.charge_sum
    end
    throw(KeyError(key.value_key))
end

export calc_V_cell_cc # TODO remove
export calc_V_cell_cc_sum # TODO remove
export calc_V_cell_cc_aging # TODO remove
export calc_V_cell_cc_sum_aging # TODO remove
export Battery

