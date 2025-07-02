"""
Control module to negotiate temperatures between a source and a target when both
have temperature-dependent energy input/output.
This control module is located always at the source component.

Currently implemented for
    Sources: GeothermalProbes
    Targets: SeasonalThermalStorage

Source components need to have the following functions:
    - get_output_temperature_bounds()
    - check_temperature_and_get_max_energy()
    - calculate_output_energy_from_output_temperature() (for optimize_for_max_energy only)

Target components need to have the following functions:
    - get_input_temperature_bounds()
    - calculate_input_energy_from_input_temperature() (for optimize_for_max_energy only)

"""

import Optim

mutable struct CM_Negotiate_Temperature <: ControlModule
    name::String
    parameters::Dict{String,Any}

    function CM_Negotiate_Temperature(parameters::Dict{String,Any},
                                      components::Grouping,
                                      sim_params::Dict{String,Any},
                                      source_uac::String)
        default_parameters = Dict{String,Any}(
            "name" => "negotiate_temperature",
            "source_uac" => source_uac,
            "target_uac" => nothing,
            "temperature_mode" => "mean",
            "limit_max_output_energy_to_avoid_pulsing" => true,
            "constant_output_temperature" => nothing,
            "optim_temperature_rel_tol" => 1e-5,   # Looser relative tolerance: 1e-3
            "optim_temperature_abs_tol" => 0.001,  # Looser absolute tolerance: 0.1
        )
        params = Base.merge(default_parameters, parameters)

        # check if the temperature mode is allowed
        allowed_parameter_temperature_mode = ["constant_temperature",
                                              "optimize_for_max_energy",
                                              "mean",
                                              "upper",
                                              "lower"]
        if !(params["temperature_mode"] in allowed_parameter_temperature_mode)
            @error "The temperature mode of control module negotiate_temperature of component $(params["source_uac"]) " *
                   "has to be one of $(allowed_parameter_temperature_mode)."
            throw(InputError)
        end

        # for temperature mode constant_temperature, check if a constant_output_temperature is provided
        if params["temperature_mode"] == "constant_temperature" && params["constant_output_temperature"] === nothing
            @error "The control module negotiate_temperature of component $(params["source_uac"]) " *
                   "is set to constant_temperature, but no `constant_output_temperature` is given!"
            throw(InputError)
        end

        # check if a valid target uac is specified
        if !(params["target_uac"] !== nothing
             && params["target_uac"] in keys(components)
             && components[params["target_uac"]] isa TemperatureNegotiateTarget)
            @error "The target of the control module negotiate_temperature of component $(params["source_uac"]) is " *
                   "not a valid component for this control module!"
            throw(InputError)
        end

        # check if a valid source uac is specified
        if !(params["source_uac"] !== nothing
             && params["source_uac"] in keys(components)
             && components[params["source_uac"]] isa TemperatureNegotiateSource)
            @error "The source of the control module negotiate_temperature of component $(params["source_uac"]) " *
                   "is not a valid source for this control module"
            throw(InputError)
        end

        return new("negotiate_temperature", params)
    end
end

function has_method_for(mod::CM_Negotiate_Temperature, func::ControlModuleFunction)::Bool
    return func == cmf_negotiate_temperature
end

function update(mod::CM_Negotiate_Temperature)
    # nothing to do
end

"""
    determine_temperature_and_energy(mod::CM_Negotiate_Temperature,
                                     components::Grouping,
                                     source_uac::String,
                                     target_uac::String,
                                     sim_params::Dict{String,Any})

Depending on the control module, this function determines the temperauture and the maximum energy that can be extracted
from a given source. Several options are available:
- Optimisation: Here, the temperature is optimized (within a given bound) to get the maximum possible energy between 
                the source and the target component.
- Constant temperature: The desired output temperature is set to a constant value by the input parameter.
- Dynamic temperature: The desired output temperature is set to a dynamic value, e.g. the temperature of a given layer
                       in a connected STES.
- Mean temperature: The desired output temperature is set to the mean value of the upper and lower temperature bound,
                    calculated from the source and the target component. This is computationally much cheaper than 
                    optimization but can lead to quite good results.

# Arguments
- `mod::CM_Negotiate_Temperature`: The current control module
- `components::Grouping`: All components of the energy system.
- `source_uac::String`: The source unit uac.
- `target_uac::String`: The target unit uac.
- `sim_params::Dict{String,Any}`: Simulation parameters.

# Returns
- The resulting maximum possible energy and the corresponding temperature.
"""
function determine_temperature_and_energy(mod::CM_Negotiate_Temperature,
                                          components::Grouping,
                                          source_uac::String,
                                          target_uac::Stringing,
                                          sim_params::Dict{String,Any})::Tuple{Temperature,Float64}
    # get temperature bounds from source and target
    # would be the same as from exchange, but control might not performed yet
    target_min_in_temperature, target_max_in_temperature = get_input_temperature_bounds(components[target_uac])
    source_min_out_temperature, source_max_out_temperature = get_output_temperature_bounds(components[source_uac], sim_params)

    # get temperature bounds for connection from source_uac to target
    lower_temperature_bound = highest(source_min_out_temperature, target_min_in_temperature)
    upper_temperature_bound = lowest(source_max_out_temperature, target_max_in_temperature)

    if upper_temperature_bound < lower_temperature_bound # both should be a number here
        # temperatures do not match, no energy can be delivered
        return nothing, 0.0
    end

    if mod.parameters["temperature_mode"] == "optimize_for_max_energy"
        return find_best_temperature_and_get_energy(mod,
                                                    calculate_output_energy_from_output_temperature,
                                                    calculate_input_energy_from_input_temperature,
                                                    lower_temperature_bound,
                                                    upper_temperature_bound,
                                                    components[source_uac],
                                                    components[target_uac],
                                                    sim_params)
    else
        # call the control module to calculate the temperature
        if mod.parameters["temperature_mode"] == "constant_temperature"
            temperature_output = Float64(mod.parameters["constant_output_temperature"])
            # elseif mod.parameters["temperature_mode"] == "STES_layer"
            #     temperature_output = components[dynamic_reference[1]].temperature_segments[dynamic_reference[2]]
        elseif mod.parameters["temperature_mode"] == "mean"
            temperature_output = (lower_temperature_bound + upper_temperature_bound) / 2
        elseif mod.parameters["temperature_mode"] == "upper"
            temperature_output = upper_temperature_bound
        elseif mod.parameters["temperature_mode"] == "lower"
            temperature_output = lower_temperature_bound
        end

        # check if the temperature is within the bounds and calculate the corresponding energy
        return check_temperature_and_get_max_energy(components[source_uac],
                                                    sim_params,
                                                    temperature_output,
                                                    mod.parameters["limit_max_output_energy_to_avoid_pulsing"])
    end
end

function find_best_temperature_and_get_energy(mod::CM_Negotiate_Temperature,
                                              func_output::Function,
                                              func_input::Function,
                                              temp_min::Temperature,
                                              temp_max::Temperature,
                                              unit_output::Component,
                                              unit_input::Component,
                                              sim_params::Dict{String,Any})::Tuple{Temperature,Float64}
    function f_min(temperature, unit_output, unit_input, sim_params)
        temp::Temperature = temperature
        return_1 = func_output(unit_output, temp, sim_params)
        return_2 = func_input(unit_input, temp, sim_params)
        return min(return_1, return_2)
    end

    # Find indices where energy is significantly non-zero to get proper bounds
    step_size = 2
    n_steps = max(6, Int(ceil((temp_max - temp_min) / step_size)))
    temps = range(temp_min, temp_max; length=n_steps)
    energies = [f_min(T, unit_output, unit_input, sim_params) for T in temps]
    active_indices = findall(energy -> energy > sim_params["epsilon"], energies)

    if isempty(active_indices)
        return nothing, 0.0
    else
        temp_min_current = max(temps[first(active_indices)] - step_size, temp_min)
        temp_max_current = min(temps[last(active_indices)] + step_size, temp_max)

        result = Optim.optimize(temperature -> -f_min(temperature, unit_output, unit_input, sim_params),
                          temp_min_current,
                          temp_max_current,
                          Optim.Brent();
                          rel_tol=mod.parameters["optim_temperature_rel_tol"],
                          abs_tol=mod.parameters["optim_temperature_abs_tol"],
                          show_trace=false)
        temperature_opt = Optim.minimizer(result)
        max_energy = f_min(temperature_opt, unit_output, unit_input, sim_params)
        return temperature_opt, max_energy
    end
end
