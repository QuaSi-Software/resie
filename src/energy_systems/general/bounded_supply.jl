"""
Implementation of a component modeling an abstract bounded supply of some medium.

This is particularly useful for testing, but can also be used to model any bounded
component or other equipment unit that processes energy in a given medium. The component
might still have a maximum power draw in a single time step, but can provide any fraction
of this to connected components.
"""
mutable struct BoundedSupply <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction
    medium::Symbol

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    max_power_profile::Union{Profile,Nothing}
    temperature_profile::Union{Profile,Nothing}
    scaling_factor::Float64

    max_energy::Float64
    temperature::Temperature
    constant_power::Union{Nothing,Float64}
    constant_temperature::Temperature

    function BoundedSupply(uac::String, config::Dict{String,Any}, parameters::Dict{String,Any})
        max_power_profile = "max_power_profile_file_path" in keys(config) ?
                            Profile(config["max_power_profile_file_path"], parameters) :
                            nothing

        # check input
        if (    haskey(config, "temperature_profile_file_path") + 
                haskey(config, "temperature_from_global_file") + 
                haskey(config, "constant_temperature")
            ) > 1
            println("Warning: Two or more temperature profile sources for $(uac) have been specified in the input file!")
        end

        # read temperature file
        if haskey(config,"temperature_profile_file_path")
            temperature_profile = Profile(config["temperature_profile_file_path"], parameters) 
            # println("Info: For bounded supply '$uac', the temperature profile is taken from the user-defined .prf file.")
        elseif haskey(config, "constant_temperature")
            temperature_profile = nothing
            # println("Info: For bounded supply '$uac', a constant temperature of $(config["constant_temperature"]) °C is set.")
        elseif haskey(config, "temperature_from_global_file") && haskey(parameters, "weatherdata")
            if any(occursin(config["temperature_from_global_file"], string(field_name)) for field_name in fieldnames(typeof(parameters["weatherdata"])))
                temperature_profile = getfield(parameters["weatherdata"], Symbol(config["temperature_from_global_file"]))
                # println("Info: For bounded supply '$uac', the temperature profile is taken from the project-wide weather file: $(config["temperature_from_global_file"])")
            else
                print("Error: For bounded supply '$uac', the'temperature_from_global_file' has to be one of: $(join(string.(fieldnames(typeof(parameters["weatherdata"]))), ", ")).")
                exit()
            end
        else
            temperature_profile = nothing
            # println("Info: For bounded supply '$uac', no temperature is set.")
        end

        medium = Symbol(config["medium"])
        register_media([medium])

        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"], parameters
            ),
            sf_bounded_source, # sys_function
            medium, # medium
            InterfaceMap( # input_interfaces
                medium => nothing
            ),
            InterfaceMap( # output_interfaces
                medium => nothing
            ),
            max_power_profile, # max_power_profile
            temperature_profile, #temperature_profile
            config["scale"], # scaling_factor
            0.0, # max_energy
            nothing, # temperature
            default(config, "constant_power", nothing), # constant_power
            default(config, "constant_temperature", nothing), # constant_temperature
        )
    end
end

function control(
    unit::BoundedSupply,
    components::Grouping,
    parameters::Dict{String,Any}
)
    move_state(unit, components, parameters)

    if unit.constant_power !== nothing
        unit.max_energy = watt_to_wh(unit.constant_power)
    elseif unit.max_power_profile !== nothing
        unit.max_energy = unit.scaling_factor * Profiles.work_at_time(
            unit.max_power_profile, parameters["time"]
        )
    else
        unit.max_energy = 0.0
    end
    set_max_energy!(unit.output_interfaces[unit.medium], unit.max_energy)

    if unit.constant_temperature !== nothing
        unit.temperature = unit.constant_temperature
    elseif unit.temperature_profile !== nothing
        unit.temperature = Profiles.value_at_time(
            unit.temperature_profile, parameters["time"]
        )
    end
    unit.output_interfaces[unit.medium].temperature = highest_temperature(
        unit.temperature,
        unit.output_interfaces[unit.medium].temperature
    )
end

function process(unit::BoundedSupply, parameters::Dict{String,Any})
    outface = unit.output_interfaces[unit.medium]
    # 1. @TODO: if disp. sources should be allowed to load storage components, then the potential
    # must be handled here instead of being ignored
    # 2. we also ignore the temperature of the interface as the source defines that itself
    exchange = balance_on(outface, outface.target)
    if exchange.balance < 0.0
        add!(
            outface,
            min(abs(exchange.balance), unit.max_energy),
            unit.temperature
        )
    end
end

function output_values(unit::BoundedSupply)::Vector{String}
    if unit.temperature_profile === nothing && unit.constant_temperature === nothing
        return [string(unit.medium)*" OUT",
                "Max_Energy"]
    else
        return [string(unit.medium)*" OUT",
                "Max_Energy",
                "Temperature"]
    end
end

function output_value(unit::BoundedSupply, key::OutputKey)::Float64
    if key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "Max_Energy"
        return unit.max_energy
    elseif key.value_key == "Temperature"
        return unit.temperature
    end
    throw(KeyError(key.value_key))
end

export BoundedSupply