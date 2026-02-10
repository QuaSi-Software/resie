using ..Resie: get_run

"""
Thermal Booster

"""
mutable struct ThermalBooster <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    m_el_in::Symbol
    m_heat_out::Symbol
    m_heat_in::Symbol

    output_temperature::Temperature
    input_temperature::Temperature

    power_losses_factor::Float64
    heat_losses_factor::Float64

    power_el::Float64
    cp_medium_out::Float64
    terminal_dT:Float64
    demand_input_temperature:Float64

    function ThermalBooster(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        m_el_in = Symbol(default(config, "m_el_in", "m_e_ac_230v"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_ht1"))
        m_heat_in = Symbol(default(config, "m_heat_in", "m_h_w_lt1"))
        register_media([m_el_in, m_heat_out, m_heat_in])

        return new(uac,
                   Controller(default(config, "control_parameters", nothing)),
                   sf_transformer,
                   InterfaceMap(m_heat_in => nothing,
                                m_el_in => nothing),
                   InterfaceMap(m_heat_out => nothing),
                   m_el_in,
                   m_heat_out,
                   m_heat_in,
                   default(config, "output_temperature", nothing),    # [°C]
                   default(config, "input_temperature", nothing),     # [°C]
                   default(config, "power_losses_factor", 0.97),      # [-]
                   default(config, "heat_losses_factor", 0.95),       # [-]
                   config["power_th"],                                # [W]
                   default(config, "cp_medium_out", 4180),            # [J/(kgK)]
                   default(config, "terminal_dT", 2),                 # [K]
                   default(config, "demand_input_temperature", 12.0)) # [°C] TODO Could also be a profile
    end
end

mutable struct TBEnergies
    potential_el_in::Floathing
    potential_el_in_layered::Vector{<:Floathing}
    potentials_heat_in::Vector{<:Floathing}
    potentials_heat_out::Vector{<:Floathing}
    available_el_in::Float64
    in_uacs_el::Vector{<:Stringing}
    available_heat_in::Vector{<:Floathing}
    available_heat_out::Vector{<:Floathing}
    in_temps_min::Vector{<:Temperature}
    in_temps_max::Vector{<:Temperature}
    in_uacs_heat::Vector{<:Stringing}
    out_temps_min::Vector{<:Temperature}
    out_temps_max::Vector{<:Temperature}
    out_uacs::Vector{<:Stringing}
    heat_in_has_inf_energy::Bool
    heat_out_has_inf_energy::Bool
    slices_temp_el_in::Vector{Floathing}
    slices_temp_heat_in::Vector{Floathing}
    slices_temp_heat_in_temperature::Vector{Temperature}
    slices_temp_heat_in_uac::Vector{Stringing}
    slices_temp_heat_out::Vector{Floathing}
    slices_temp_heat_out_temperature::Vector{Temperature}
    slices_temp_heat_out_uac::Vector{Stringing}
    slices_temp_times::Vector{Floathing}
    slices_el_in::Vector{Floathing}
    slices_heat_in::Vector{Floathing}
    slices_heat_in_temperature::Vector{Temperature}
    slices_heat_in_uac::Vector{Stringing}
    slices_heat_out::Vector{Floathing}
    slices_heat_out_temperature::Vector{Temperature}
    slices_heat_out_uac::Vector{Stringing}

    function TBEnergies()
        return new(0.0,
                   Vector{Floathing}(),
                   Vector{Floathing}(),
                   Vector{Floathing}(),
                   0.0,
                   Vector{Stringing}(),
                   Vector{Floathing}(),
                   Vector{Floathing}(),
                   Vector{Temperature}(),
                   Vector{Temperature}(),
                   Vector{Stringing}(),
                   Vector{Temperature}(),
                   Vector{Temperature}(),
                   Vector{Stringing}(),
                   false,
                   false,
                   Vector{Floathing}(),
                   Vector{Floathing}(),
                   Vector{Temperature}(),
                   Vector{Stringing}(),
                   Vector{Floathing}(),
                   Vector{Temperature}(),
                   Vector{Stringing}(),
                   Vector{Floathing}(),
                   Vector{Floathing}(),
                   Vector{Floathing}(),
                   Vector{Temperature}(),
                   Vector{Stringing}(),
                   Vector{Floathing}(),
                   Vector{Temperature}(),
                   Vector{Stringing}())
    end
end

function initialise!(unit::ThermalBooster, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.input_interfaces[unit.m_heat_in],
                          unload_storages(unit.controller, unit.m_heat_in))
    set_storage_transfer!(unit.input_interfaces[unit.m_el_in],
                          unload_storages(unit.controller, unit.m_el_in))
    set_storage_transfer!(unit.output_interfaces[unit.m_heat_out],
                          load_storages(unit.controller, unit.m_heat_out))
end

function control(unit::ThermalBooster, components::Grouping, sim_params::Dict{String,Any})
    update(unit.controller)

    if unit.output_temperature !== nothing
        set_max_energy!(unit.output_interfaces[unit.m_heat_out],
                        nothing,
                        nothing,
                        unit.output_temperature)
    end
    if unit.input_temperature !== nothing
        set_max_energy!(unit.input_interfaces[unit.m_heat_in],
                        nothing,
                        unit.input_temperature,
                        nothing)
    end
end

function set_max_energies!(unit::ThermalBooster,
                           el_in::Union{Floathing,Vector{<:Floathing}},
                           heat_in::Union{Floathing,Vector{<:Floathing}},
                           heat_out::Union{Floathing,Vector{<:Floathing}},
                           slices_heat_in_temperature::Union{Temperature,Vector{<:Temperature}}=nothing,
                           slices_heat_out_temperature::Union{Temperature,Vector{<:Temperature}}=nothing,
                           purpose_uac_heat_in::Union{Stringing,Vector{Stringing}}=nothing,
                           purpose_uac_heat_out::Union{Stringing,Vector{Stringing}}=nothing,
                           has_calculated_all_maxima_heat_in::Bool=false,
                           has_calculated_all_maxima_heat_out::Bool=false)
    set_max_energy!(unit.input_interfaces[unit.m_el_in], el_in)
    set_max_energy!(unit.input_interfaces[unit.m_heat_in], heat_in, slices_heat_in_temperature, nothing,
                    purpose_uac_heat_in, has_calculated_all_maxima_heat_in)
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], heat_out, nothing, slices_heat_out_temperature,
                    purpose_uac_heat_out, has_calculated_all_maxima_heat_out)
end

function calculate_energies(unit::ThermalBooster, sim_params::Dict{String,Any})::TBEnergies
    energies = TBEnergies()

    # get electricity potential and reduce it by constant power draw (or however much
    # is available)
    energies.potential_el_in_layered, energies.in_uacs_el = check_el_in_layered(unit, sim_params)
    energies.potential_el_in = sum(energies.potential_el_in_layered)

    # shortcut if we're limited by zero input electricity
    if energies.potential_el_in <= 0.0
        do_calculation = false
    end

    if do_calculation
        # get vectored values for the input and output heat potentials
        energies.potentials_heat_in,
        energies.in_temps_min,
        energies.in_temps_max,
        energies.in_uacs_heat = check_heat_in_layered(unit, sim_params)

        energies.potentials_heat_out,
        energies.out_temps_min,
        energies.out_temps_max,
        energies.out_uacs = check_heat_out_layered(unit, sim_params)

        # in the following we want to work with positive values as it is easier
        energies.potentials_heat_in = abs.(energies.potentials_heat_in)
        energies.potentials_heat_out = abs.(energies.potentials_heat_out)

        # reduce available input energies by the power/heat losses that would occur if the
        # sources would be fully utilised. since the actual usage is equal or less than that,
        # it works out even if the PLR for that source is not 1.0
        energies.potentials_heat_in .*= unit.heat_losses_factor
        energies.potential_el_in *= unit.power_losses_factor

        # reorder inputs and outputs according to control modules
        index = reorder_inputs(unit.controller, energies.in_temps_min, energies.in_temps_max)
        energies.in_temps_min = energies.in_temps_min[index]
        energies.in_temps_max = energies.in_temps_max[index]
        energies.in_uacs_heat = energies.in_uacs_heat[index]
        energies.potentials_heat_in = energies.potentials_heat_in[index]

        index = reorder_outputs(unit.controller, energies.out_temps_min, energies.out_temps_max)
        energies.out_temps_min = energies.out_temps_min[index]
        energies.out_temps_max = energies.out_temps_max[index]
        energies.out_uacs = energies.out_uacs[index]
        energies.potentials_heat_out = energies.potentials_heat_out[index]

        # there are three different cases of how to handle the layered approach of operating the
        # thermal booster, depending on wether or not any input heat or output heat transformer has a
        # value of infinite as the potential. if this is the case for both the input and output,
        # the calculation cannot be resolved.
        energies.heat_in_has_inf_energy = any(isinf, filter_by_transformer(energies, sim_params; heat_in=true))
        energies.heat_out_has_inf_energy = any(isinf, filter_by_transformer(energies, sim_params; heat_in=false))

        if energies.heat_in_has_inf_energy && energies.heat_out_has_inf_energy
            @warn "The thermal booster $(unit.uac) has unknown energies in both its inputs and " *
                  "outputs. This cannot be resolved. Please check the order of operation and " *
                  "make sure that either the inputs or the outputs have been fully calculated " *
                  "before the thermal booster $(unit.uac) has its potential step."
            return energies
        end

        if energies.heat_in_has_inf_energy
            for heat_in_idx in eachindex(energies.potentials_heat_in)
                energies = calculate_booster(unit, sim_params, energies; fixed_heat_in=heat_in_idx)
            end
        elseif energies.heat_out_has_inf_energy
            for heat_out_idx in eachindex(energies.potentials_heat_out)
                energies = calculate_booster(unit, sim_params, energies; fixed_heat_out=heat_out_idx)
            end
        else
            energies = calculate_booster(unit, sim_params, energies)
        end
    end

    # now set losses of the thermal booster and add the losses to the actually consumed
    # power / heat for the slices
    el_in = sum(energies.slices_el_in; init=0.0)
    heat_in = sum(energies.slices_heat_in; init=0.0)
    unit.losses_power = -1.0 * el_in
    unit.losses_heat = -1.0 * heat_in
    energies.slices_el_in ./= unit.power_losses_factor
    energies.slices_heat_in ./= unit.heat_losses_factor

    el_in = sum(energies.slices_el_in; init=0.0)
    heat_in = sum(energies.slices_heat_in; init=0.0)
    unit.losses_power += el_in
    unit.losses_heat += heat_in

    return energies
end

function calculate_booster(unit::ThermalBooster,
                           sim_params::Dict{String,Any},
                           energies::TBEnergies;
                           fixed_heat_in=Union{Nothing,Integer}=nothing,
                           fixed_heat_out=Union{Nothing,Integer}=nothing)::TBEnergies

    ## INPUTS
    # energies.potential_el_in_layered
    # energies.in_uacs_el 
    # energies.potential_el_in 

    # energies.potentials_heat_in, # positive, reduced by losses
    # energies.in_temps_min,
    # energies.in_temps_max,
    # energies.in_uacs_heat

    # energies.potentials_heat_out,  # positive, reduced by losses
    # energies.out_temps_min,
    # energies.out_temps_max,
    # energies.out_uacs 

    # energies.heat_in_has_inf_energy
    # energies.heat_out_has_inf_energy

    ##OUTPUTS
    # energies.slices_el_in,
    # energies.slices_heat_in,
    # energies.slices_heat_out,
    # energies.slices_heat_in_temperature,
    # energies.slices_heat_out_temperature,
    # energies.slices_heat_in_uac,
    # energies.slices_heat_out_uac,
    # energies.heat_in_has_inf_energy,
    # energies.heat_out_has_inf_energy

    ##PARAMETER
    # power_el::Float64
    # cp_medium_out::Float64
    # terminal_dT:Float64
    # demand_input_temperature:Float64

end

function potential(unit::ThermalBooster, sim_params::Dict{String,Any})
    energies = calculate_energies(unit, sim_params)

    if sum(energies.slices_heat_out; init=0.0) < sim_params["epsilon"]
        set_max_energies!(unit, sum(energies.slices_el_in; init=0.0), 0.0, 0.0)
        return
    end

    set_max_energies!(unit,
                      energies.slices_el_in,
                      energies.slices_heat_in,
                      energies.slices_heat_out,
                      energies.slices_heat_in_temperature,
                      energies.slices_heat_out_temperature,
                      energies.slices_heat_in_uac,
                      energies.slices_heat_out_uac,
                      energies.heat_in_has_inf_energy,
                      energies.heat_out_has_inf_energy)
end

function process(unit::ThermalBooster, sim_params::Dict{String,Any})
    energies = calculate_energies(unit, sim_params)

    el_in = sum(energies.slices_el_in; init=0.0)
    heat_out = sum(energies.slices_heat_out; init=0.0)

    if heat_out < sim_params["epsilon"]
        set_max_energies!(unit, el_in, 0.0, 0.0, 0.0)
        sub!(unit.input_interfaces[unit.m_el_in], el_in)
    else
        sub!(unit.input_interfaces[unit.m_el_in], el_in)
        sub!(unit.input_interfaces[unit.m_heat_in],
             energies.slices_heat_in,
             energies.slices_heat_in_temperature,
             [nothing for _ in energies.slices_heat_in_temperature],
             energies.slices_heat_in_uac)
        add!(unit.output_interfaces[unit.m_heat_out],
             energies.slices_heat_out,
             [nothing for _ in energies.slices_heat_out_temperature],
             energies.slices_heat_out_temperature,
             energies.slices_heat_out_uac)
    end
    # calculate total losses
    unit.losses = unit.losses_power + unit.losses_hea
end

# has its own reset function as here more parameters are present that need to be reset in
# every timestep
function reset(unit::ThermalBooster)
    invoke(reset, Tuple{Component}, unit)

    # reset other parameter
    unit.losses = 0.0
end

function output_values(unit::ThermalBooster)::Vector{String}
    output_vals = [string(unit.m_el_in) * ":IN",
                   string(unit.m_heat_in) * ":IN",
                   string(unit.m_heat_out) * ":OUT"]
    append!(output_vals, ["Losses", "InputTemperature", "OutputTemperature"])
    return output_vals
end

function output_value(unit::ThermalBooster, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "Losses"
        return -unit.losses
    end
    throw(KeyError(key.value_key))
end

export ThermalBooster
