"""
Implementation of a boiler producing heat (as hot water) from the chemical energy in a fuel.

While not technically correct, this implementation can also be used to model an electric
boiler, as the input medium can be freely chosen and therefore can be chosen to be
electricity instead of a chemical fuel.
"""
mutable struct FuelBoiler <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    m_fuel_in::Symbol
    m_heat_out::Symbol

    power_th::Float64
    is_plr_dependant::Bool
    max_consumable_fuel::Float64
    plr_to_expended_energy::Vector{Tuple{Float64,Float64}}
    min_power_fraction::Float64 # Minimum amount of power so that the component can be
    min_run_time::UInt          # operated
    output_temperature::Temperature

    function FuelBoiler(uac::String, config::Dict{String,Any})
        m_fuel_in = Symbol(config["m_fuel_in"])
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_ht1"))
        register_media([m_fuel_in, m_heat_out])

        max_consumable_fuel = watt_to_wh(float(config["power_th"])) /
                             default(config, "max_thermal_efficiency", 1.0)
        # lookup table for conversion of part load ratio to expended energy
        plr_to_expended_energy = []
        # fill up plr_to_expended_energy lookup table
        start_value = 0.0
        end_value = 1.0
        step_size = 0.1 # set the discretization
        for plr in collect(start_value:step_size:end_value)
            # Create a tuple: (part load ratio value, expended energy); hard coded function
            # for expended energy = useful_energy / thermal_efficiency
            plr_expended_energy_pair = (plr,
                (plr / (-0.9117 * plr^2 + 1.8795 * plr + 0.0322)) * 1000)
            # Append the tuple to the lookup table
            push!(plr_to_expended_energy, plr_expended_energy_pair)
        end

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
            ),
            sf_transformer, # sys_function
            InterfaceMap( # input_interfaces
                m_fuel_in => nothing
            ),
            InterfaceMap( # output_interfaces
                m_heat_out => nothing
            ),
            m_fuel_in,
            m_heat_out,
            config["power_th"], # power_th
            default(config, "is_plr_dependant", false), # toggles PLR-dependant efficiency
            max_consumable_fuel,
            plr_to_expended_energy,
            default(config, "min_power_fraction", 0.1),
            default(config, "min_run_time", 0),
            default(config, "output_temperature", nothing)
        )
    end
end

function control(
    unit::FuelBoiler,
    components::Grouping,
    parameters::Dict{String,Any}
)
    move_state(unit, components, parameters)
    unit.output_interfaces[unit.m_heat_out].temperature = highest_temperature(
        unit.output_temperature,
        unit.output_interfaces[unit.m_heat_out].temperature
    )
end

"""
Set maximum energies that can be taken in and put out by the unit
"""
function set_max_energies!(unit::FuelBoiler, fuel_in::Float64, heat_out::Float64)
    set_max_energy!(unit.input_interfaces[unit.m_fuel_in], fuel_in)
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], heat_out)
end

function check_fuel_in(
    unit::FuelBoiler,
    parameters::Dict{String,Any}
)
    if unit.controller.parameter["consider_m_fuel_in"] == true
        if (
            unit.input_interfaces[unit.m_fuel_in].source.sys_function == sf_transformer
            &&
            unit.input_interfaces[unit.m_fuel_in].max_energy === nothing
        )
            return (Inf, Inf)
        else
            exchange = balance_on(
                unit.input_interfaces[unit.m_fuel_in],
                unit.input_interfaces[unit.m_fuel_in].source
            )
            potential_energy_fuel = exchange.balance + exchange.energy_potential
            potential_storage_fuel = exchange.storage_potential
            if (
                unit.controller.parameter["unload_storages"] ?
                potential_energy_fuel + potential_storage_fuel :
                potential_energy_fuel
            ) <= parameters["epsilon"]
                return (nothing, nothing)
            end
            return (potential_energy_fuel, potential_storage_fuel, exchange.temperature)
        end
    else
        return (Inf, Inf)
    end
end

function check_heat_out(
    unit::FuelBoiler,
    parameters::Dict{String,Any}
)
    if unit.controller.parameter["consider_m_heat_out"] == true
        exchange = balance_on(
            unit.output_interfaces[unit.m_heat_out],
            unit.output_interfaces[unit.m_heat_out].target
        )
        potential_energy_heat_out = exchange.balance + exchange.energy_potential
        potential_storage_heat_out = exchange.storage_potential
        if (
            unit.controller.parameter["load_storages"] ?
            potential_energy_heat_out + potential_storage_heat_out :
            potential_energy_heat_out
        ) >= -parameters["epsilon"]
            return (nothing, nothing)
        end
        return (potential_energy_heat_out, potential_storage_heat_out)
    else
        return (-Inf, -Inf)
    end
end

"""
This function, with set magic-numbers, serves as a temporary implementation of an
efficiency/PLR curve until a more generalised framework for such calculations is
implemented.
"""
function calculate_thermal_efficiency(
    plr::Float64, # part load ratio
)
    return -0.9117 * plr^2 + 1.8795 * plr + 0.0322
end

function calculate_plr_through_inversed_expended_energy(
    intake_fuel::Float64,
    unit::FuelBoiler
)
    # Variables to define interval where provided value falls in
    lower_x = nothing
    lower_y = nothing
    upper_x = nothing
    upper_y = nothing

    # Iterate through the lookup table to identify values of interval
    for (x, y) in unit.plr_to_expended_energy
        if y <= intake_fuel
            lower_x = x
            lower_y = y
        else
            upper_x = x
            upper_y = y
            break  # Once we exceed the target value, exit the loop
        end
    end

    # If intake_fuel equals max_consumable_fuel, then fuel boiler is working at full mode
    if isnothing(upper_x) && unit.plr_to_expended_energy[end][2] == intake_fuel
        return 1.0
    end

    try  # Check if intake_fuel is within the range of the lookup table
        if lower_x !== nothing && upper_x !== nothing
            # perform linear interpolation
            return lower_x + (intake_fuel - lower_y) *
                             (upper_x - lower_x) / (upper_y - lower_y)
        end
    catch
        println("The intake_fuel value is not within the range of the lookup table.")
    end
end

function calculate_energies(
    unit::FuelBoiler,
    parameters::Dict{String,Any},
    potentials::Vector{Float64}
)
    potential_energy_fuel_in = potentials[1]
    potential_storage_fuel_in = potentials[2]
    potential_energy_heat_out = potentials[3]
    potential_storage_heat_out = potentials[4]

    if unit.is_plr_dependant == false
        max_produce_heat = watt_to_wh(unit.power_th)
        max_consume_fuel = max_produce_heat
    elseif unit.is_plr_dependant == true # part load ratio
        if unit.controller.strategy == "demand_driven"
            # max_produce_heat should equal rated power_th, else when using demand heat, which
            # can be inf, a non-physical state might occur
            max_produce_heat = watt_to_wh(unit.power_th)
            demand_heat = -(unit.controller.parameter["load_storages"] ?
                            potential_energy_heat_out + potential_storage_heat_out :
                            potential_energy_heat_out)
            if demand_heat >= max_produce_heat
                # set demand_heat to rated power_th -> prevent max_consume_fuel to equal inf
                demand_heat = max_produce_heat
            elseif demand_heat < 0
                # ensure that part-load-ratio is greater equal 0
                demand_heat = 0
            end
            part_load_ratio = demand_heat / watt_to_wh(unit.power_th)
            thermal_efficiency = calculate_thermal_efficiency(part_load_ratio)
            max_consume_fuel = max_produce_heat / thermal_efficiency
        elseif unit.controller.strategy == "supply_driven"
            intake_fuel = unit.input_interfaces[unit.m_fuel_in].max_energy
            max_consume_fuel = unit.max_consumable_fuel
            part_load_ratio = calculate_plr_through_inversed_expended_energy(
                intake_fuel, unit
            )
            thermal_efficiency = calculate_thermal_efficiency(part_load_ratio)
            max_produce_heat = thermal_efficiency * max_consume_fuel
        end
    end

    # get usage fraction of external profile (normalized from 0 to 1)
    # when no profile is provided, then 100 % of unit is used
    usage_fraction_operation_profile =
        unit.controller.parameter["operation_profile_path"] === nothing ?
        1.0 :
        value_at_time(unit.controller.parameter["operation_profile"], parameters["time"])
    if usage_fraction_operation_profile <= 0.0
        return # no operation allowed from external profile
    end

    # all three standard operating strategies behave the same, but it is better to be
    # explicit about the behaviour rather than grouping all together
    if (
        unit.controller.strategy == "storage_driven" &&
        unit.controller.state_machine.state == 2
    )
        usage_fraction_heat_out = -((unit.controller.parameter["load_storages"] ?
                                     potential_energy_heat_out + potential_storage_heat_out :
                                     potential_energy_heat_out) / max_produce_heat)
        usage_fraction_fuel_in = +((unit.controller.parameter["unload_storages"] ?
                                   potential_energy_fuel_in + potential_storage_fuel_in :
                                   potential_energy_fuel_in) / max_consume_fuel)

    elseif unit.controller.strategy == "storage_driven"
        return (false, nothing, nothing, nothing)

    elseif unit.controller.strategy == "supply_driven"
        usage_fraction_heat_out = -((unit.controller.parameter["load_storages"] ?
                                     potential_energy_heat_out + potential_storage_heat_out :
                                     potential_energy_heat_out) / max_produce_heat)
        usage_fraction_fuel_in = +((unit.controller.parameter["unload_storages"] ?
                                   potential_energy_fuel_in + potential_storage_fuel_in :
                                   potential_energy_fuel_in) / max_consume_fuel)

    elseif unit.controller.strategy == "demand_driven"
        usage_fraction_heat_out = -((unit.controller.parameter["load_storages"] ?
                                     potential_energy_heat_out + potential_storage_heat_out :
                                     potential_energy_heat_out) / max_produce_heat)
        usage_fraction_fuel_in = +((unit.controller.parameter["unload_storages"] ?
                                   potential_energy_fuel_in + potential_storage_fuel_in :
                                   potential_energy_fuel_in) / max_consume_fuel)
    end

    # limit actual usage by limits of inputs, outputs and profile
    usage_fraction = min(
        1.0,
        usage_fraction_heat_out,
        usage_fraction_fuel_in,
        usage_fraction_operation_profile
    )

    # Component is not used if usage_fraction is less than min_power_fraction that is
    # required to use the component
    if usage_fraction < unit.min_power_fraction
        return (false, nothing, nothing)
    end

    return (
        true,
        max_consume_fuel * usage_fraction,
        max_produce_heat * usage_fraction
    )
end

function potential(
    unit::FuelBoiler,
    parameters::Dict{String,Any}
)
    potential_energy_fuel_in, potential_storage_fuel_in = check_fuel_in(unit, parameters)
    if potential_energy_fuel_in === nothing && potential_storage_fuel_in === nothing
        set_max_energies!(unit, 0.0, 0.0)
        return
    end

    potential_energy_heat_out, potential_storage_heat_out = check_heat_out(unit, parameters)
    if potential_energy_heat_out === nothing && potential_storage_heat_out === nothing
        set_max_energies!(unit, 0.0, 0.0)
        return
    end

    energies = calculate_energies(
        unit, parameters,
        [
            potential_energy_fuel_in, potential_storage_fuel_in,
            potential_energy_heat_out, potential_storage_heat_out
        ]
    )

    if !energies[1]
        set_max_energies!(unit, 0.0, 0.0)
    else
        set_max_energies!(unit, energies[2], energies[3])
    end
end

function process(unit::FuelBoiler, parameters::Dict{String,Any})
    potential_energy_fuel_in, potential_storage_fuel_in = check_fuel_in(unit, parameters)
    if potential_energy_fuel_in === nothing && potential_storage_fuel_in === nothing
        set_max_energies!(unit, 0.0, 0.0)
        return
    end

    potential_energy_heat_out, potential_storage_heat_out = check_heat_out(unit, parameters)
    if potential_energy_heat_out === nothing && potential_storage_heat_out === nothing
        set_max_energies!(unit, 0.0, 0.0)
        return
    end

    energies = calculate_energies(
        unit, parameters,
        [
            potential_energy_fuel_in, potential_storage_fuel_in,
            potential_energy_heat_out, potential_storage_heat_out
        ]
    )
    if energies[1]
        sub!(unit.input_interfaces[unit.m_fuel_in], energies[2])
        add!(unit.output_interfaces[unit.m_heat_out], energies[3])
    else
        set_max_energies!(unit, 0.0, 0.0)
    end
end

export FuelBoiler