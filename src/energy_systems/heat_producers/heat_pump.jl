using Optim: optimize, minimizer, Options, NelderMead
using ..Resie: get_run

#! format: off
const HEAT_PUMP_PARAMETERS = Dict(
    "m_el_in" => (
        default="m_e_ac_230v",
        description="Electrical input medium",
        display_name="Medium el_in",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "m_heat_in" => (
        default="m_h_w_lt1",
        description="Heat input medium",
        display_name="Medium heat_in",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "m_heat_out" => (
        default="m_h_w_ht1",
        description="Heat output medium",
        display_name="Medium heat_out",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "m_heat_out_secondary" => (
        default="m_h_w_ht1",
        description="Heat output medium on secondary interface",
        display_name="Medium heat_out_secondary",
        required=false,
        type=String,
        json_type="string",
        conditionals=[("has_secondary_interface", "is_true")],
        unit="-"
    ),
    "has_secondary_interface" => (
        default=false,
        description="Toggles the secondary heat output interface on/off",
        display_name="Has secondary interface?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "primary_el_sources" => (
        default=Vector{String}(),
        description="Electricity sources for the primary heat output",
        display_name="Primary electricity sources",
        required=false,
        type=Vector{String},
        json_type="array",
        conditionals=[("has_secondary_interface", "is_true")],
        unit="-"
    ),
    "secondary_el_sources" => (
        default=Vector{String}(),
        description="Electricity sources for the secondary heat output",
        display_name="Secondary electricity sources",
        required=false,
        type=Vector{String},
        json_type="array",
        conditionals=[("has_secondary_interface", "is_true")],
        unit="-"
    ),
    "cop_function" => (
        default="carnot:0.4",
        description="COP function definition (e.g., 'carnot:0.4' or 'const:3.0')",
        display_name="COP function",
        required=false,
        type=String,
        json_type="string",
        function_type="cop",
        unit="-"
    ),
    "plf_function" => (
        default="const:1.0",
        description="Part-load factor function definition",
        display_name="PLF function",
        required=false,
        type=String,
        json_type="string",
        function_type="1dim",
        unit="-"
    ),
    "max_power_function" => (
        default="const:1.0",
        description="Maximum power function as fraction of design power",
        display_name="Max. power function",
        required=false,
        type=String,
        json_type="string",
        function_type="2dim",
        unit="-"
    ),
    "min_power_function" => (
        default="const:0.2",
        description="Minimum power function as fraction of design power",
        display_name="Min. power function",
        required=false,
        type=String,
        json_type="string",
        function_type="2dim",
        unit="-"
    ),
    "icing_coefficients" => (
        default="3,-0.42,15,2,30",
        description="Icing loss coefficients (comma-separated)",
        display_name="Icing coefficients",
        required=false,
        type=String,
        json_type="string",
        unit="-"
    ),
    "model_type" => (
        default="simplified",
        description="Operation model: 'simplified', 'on-off', or 'inverter'",
        display_name="Model type",
        required=false,
        type=String,
        json_type="string",
        options=["simplified", "on-off", "inverter"],
        unit="-"
    ),
    "power_th" => (
        default=nothing,
        description="Design thermal power",
        display_name="Thermal power",
        required=true,
        type=Float64,
        json_type="number",
        unit="W"
    ),
    "bypass_cop" => (
        default=15.0,
        description="COP in bypass operation (when input >= output temperature)",
        display_name="Bypass COP",
        required=false,
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "min_usage_fraction" => (
        default=0.0,
        description="Minimum part-load ratio to operate",
        display_name="Min. usage fraction",
        required=false,
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "consider_icing" => (
        default=false,
        description="Toogle to account for icing losses on heat input side",
        display_name="Consider icing?",
        required=false,
        type=Bool,
        json_type="boolean",
        unit="-"
    ),
    "output_temperature" => (
        default=nothing,
        description="Fixed output temperature, or nothing for auto-detection",
        display_name="Output temperature",
        required=false,
        type=Floathing,
        json_type="number",
        unit="°C"
    ),
    "input_temperature" => (
        default=nothing,
        description="Fixed input temperature, or nothing for auto-detection",
        display_name="Input temperature",
        required=false,
        type=Floathing,
        json_type="number",
        unit="°C"
    ),
    "nr_optimisation_passes" => (
        default=20,
        description="Maximum iterations for PLR optimisation",
        display_name="Nr. optimisation passes",
        required=false,
        type=UInt,
        json_type="number",
        conditionals=[("model_type", "is_one_of", ("inverter", "on-off"))],
        unit="-"
    ),
    "eval_factor_heat" => (
        default=5.0,
        description="Weight for heat demand matching in optimization",
        display_name="Eval. factor heat",
        required=false,
        type=Float64,
        json_type="number",
        conditionals=[("model_type", "is_one_of", ("inverter", "on-off"))],
        unit="-"
    ),
    "eval_factor_time" => (
        default=1.0,
        description="Weight for time usage in optimisation",
        display_name="Eval. factor time",
        required=false,
        type=Float64,
        json_type="number",
        conditionals=[("model_type", "is_one_of", ("inverter", "on-off"))],
        unit="-"
    ),
    "eval_factor_elec" => (
        default=1.0,
        description="Weight for electricity minimization in optimisation",
        display_name="Eval. factor electricity",
        required=false,
        type=Float64,
        json_type="number",
        conditionals=[("model_type", "is_one_of", ("inverter", "on-off"))],
        unit="-"
    ),
    "fudge_factor" => (
        default=1.001,
        description="Overestimation factor for available heat output in optimisation",
        display_name="Fudge factor",
        required=false,
        type=Float64,
        json_type="number",
        conditionals=[("model_type", "is_one_of", ("inverter", "on-off"))],
        unit="-"
    ),
    "x_abstol" => (
        default=0.01,
        description="Absolute tolerance for PLR optimization",
        display_name="Abs. tolerance PLR",
        required=false,
        type=Float64,
        json_type="number",
        conditionals=[("model_type", "is_one_of", ("inverter", "on-off"))],
        unit="-"
    ),
    "f_abstol" => (
        default=0.001,
        description="Objective function tolerance for PLR optimization",
        display_name="Abs. tolerance obj. function",
        required=false,
        type=Float64,
        json_type="number",
        conditionals=[("model_type", "is_one_of", ("inverter", "on-off"))],
        unit="-"
    ),
    "power_losses_factor" => (
        default=0.97,
        description="Fraction of input power used (losses: 1 - factor)",
        display_name="Power losses factor",
        required=false,
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "heat_losses_factor" => (
        default=0.95,
        description="Fraction of input heat used (losses: 1 - factor)",
        display_name="Heat losses factor",
        required=false,
        type=Float64,
        json_type="number",
        unit="-"
    ),
    "constant_loss_power" => (
        default=0.0,
        description="Constant power loss, e.g. standby power",
        display_name="Constant power loss",
        required=false,
        type=Float64,
        json_type="number",
        unit="W"
    ),
)
#! format: on

"""
Implementation of a heat pump component, elevating heat to a higher temperature using
electricity.

This models includes several effects that complicate the calculation of efficiency (as COP),
some of which involve the temperatures at which the heat input is provided and heat output
is requested. While this is similar to the calculation of efficiencies in other transformers
such as a CHPP, the implementation does not numerically inverse the given input functions
and determine energies from those. Instead the efficiencies are calculated "forward", which
incurs the problem of having to determine kappa (the PLR). This might be improved in the
future for heat pumps serving exactly one demand using exactly one source.

The heat pump can also be configured to use on of three models: Inverter-driven, on-off and
simplified heat pumps. The first two require an optimisation algorithm to determine the
operation point, while the simplified skips this process and always uses kappa = 1.0.
"""
mutable struct HeatPump <: Component
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap

    m_el_in::Symbol
    m_heat_out::Symbol
    m_heat_out_secondary::Symbol
    m_heat_in::Symbol

    has_secondary_interface::Bool
    primary_el_sources::Vector{String}
    secondary_el_sources::Vector{String}

    model_type::String

    design_power_th::Float64
    max_power_function::Function
    min_power_function::Function
    plf_function::Function
    constant_cop::Floathing
    dynamic_cop::Function
    bypass_cop::Float64
    min_usage_fraction::Float64

    consider_icing::Bool
    icing_coefficients::Vector{Float64}
    output_temperature::Temperature
    input_temperature::Temperature

    nr_optimisation_passes::UInt
    eval_factor_heat::Float64
    eval_factor_time::Float64
    eval_factor_elec::Float64
    fudge_factor::Float64
    x_abstol::Float64
    f_abstol::Float64

    power_losses_factor::Float64
    heat_losses_factor::Float64
    constant_loss_energy::Float64
    current_constant_loss::Float64
    losses_power::Float64
    losses_heat::Float64
    losses::Float64

    cop::Float64
    effective_cop::Float64
    mix_temp_input::Float64
    mix_temp_output::Float64
    avg_plr::Float64
    time_active::Float64

    function HeatPump(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        constructor_errored = false

        # extract all parameters using the parameter dictionary as the source of truth
        extracted_params = Dict{String,Any}()
        for (param_name, param_def) in HEAT_PUMP_PARAMETERS
            try
                extracted_params[param_name] = extract_parameter(HeatPump, config, param_name, param_def, sim_params)
            catch e
                @error "$(sprint(showerror, e))"
                constructor_errored = true
            end
        end

        try
            # extract control_parameters, which is essentially a subconfig
            extracted_params["control_parameters"] = extract_control_parameters(Component, config)

            # validate configuration, e.g. for interdependencies and allowed values
            validate_config(HeatPump, config, extracted_params, uac, sim_params)
        catch e
            @error "$(sprint(showerror, e))"
            constructor_errored = true
        end

        # we delayed throwing the errors upwards so that many errors are caught at once
        if constructor_errored
            throw(InputError("Can't construct component $uac because of errors parsing " *
                             "the parameter config. Check the error log for more " *
                             "information on each error."))
        end

        # initialize and construct the object
        init_values = init_from_params(HeatPump, uac, extracted_params, config)
        return new(init_values...)
    end
end

function component_parameters(x::Type{HeatPump})::Dict{String,NamedTuple}
    return deepcopy(HEAT_PUMP_PARAMETERS) # Return a copy to prevent external modification
end

function extract_parameter(x::Type{HeatPump}, config::Dict{String,Any}, param_name::String,
                           param_def::NamedTuple, sim_params::Dict{String,Any})
    if param_name == "icing_coefficients"
        value = default(config, param_name, param_def.default)
        return parse.(Float64, split(value, ","))
    end

    if param_name == "constant_loss_power"
        value = default(config, param_name, param_def.default)
        return sim_params["watt_to_wh"](value)
    end

    return extract_parameter(Component, config, param_name, param_def, sim_params)
end

function validate_config(x::Type{HeatPump}, config::Dict{String,Any}, extracted::Dict{String,Any},
                         uac::String, sim_params::Dict{String,Any})
    validate_config(Component, extracted, uac, sim_params, component_parameters(HeatPump))

    model_type = extracted["model_type"]
    cop_func_def = default(config, "cop_function", HEAT_PUMP_PARAMETERS["cop_function"].default)
    plf_func_def = default(config, "plf_function", HEAT_PUMP_PARAMETERS["plf_function"].default)

    if model_type in ("inverter", "on-off") && occursin("const", cop_func_def)
        @error "Heat pump $(uac) is configured to use optimisation for inverter-driven" *
               "or on-off operation, but has a constant COP. Toggle optimisation off " *
               "by switching to simplified model type as the algorithm is unstable " *
               "this case."
        throw(InputError())
    end

    if model_type == "simplified" && !occursin("const", plf_func_def)
        @error "Heat pump $(uac) has model type simplified and a non-constant PLF " *
               "function. The simplified model cannot handle this correctly. Please " *
               "use a different model type or switch to a constant PLF function."
        throw(InputError())
    end
end

function init_from_params(x::Type{HeatPump}, uac::String, params::Dict{String,Any}, raw_params::Dict{String,Any})::Tuple
    # turn media names into Symbol and register them
    m_el_in = Symbol(params["m_el_in"])
    m_heat_out = Symbol(params["m_heat_out"])
    m_heat_out_secondary = Symbol(adjust_name_if_secondary(m_heat_out, true))
    m_heat_in = Symbol(params["m_heat_in"])
    register_media([m_el_in, m_heat_out, m_heat_in, m_heat_out_secondary])

    # set constant COP if the function is defined that way
    constant_cop = nothing
    func_def = default(raw_params, "cop_function", HEAT_PUMP_PARAMETERS["cop_function"].default)
    if occursin("const", func_def)
        constant_cop = params["cop_function"](0.0, 0.0)
    end

    output_interfaces = InterfaceMap(m_heat_out => nothing)
    if params["has_secondary_interface"]
        output_interfaces[m_heat_out_secondary] = nothing
    end

    # return tuple in the order expected by new()
    return (uac,
            Controller(params["control_parameters"]),
            sf_transformer,
            InterfaceMap(m_heat_in => nothing, m_el_in => nothing),
            output_interfaces,
            m_el_in,
            m_heat_out,
            m_heat_out_secondary,
            m_heat_in,
            params["has_secondary_interface"],
            params["primary_el_sources"],
            params["secondary_el_sources"],
            params["model_type"],
            params["power_th"],
            params["max_power_function"],
            params["min_power_function"],
            params["plf_function"],
            constant_cop,
            params["cop_function"],
            params["bypass_cop"],
            params["min_usage_fraction"],
            params["consider_icing"],
            params["icing_coefficients"],
            params["output_temperature"],
            params["input_temperature"],
            params["nr_optimisation_passes"],
            params["eval_factor_heat"],
            params["eval_factor_time"],
            params["eval_factor_elec"],
            params["fudge_factor"],
            params["x_abstol"],
            params["f_abstol"],
            params["power_losses_factor"],
            params["heat_losses_factor"],
            params["constant_loss_power"],
            0.0,  # current_constant_loss
            0.0,  # losses_power
            0.0,  # losses_heat
            0.0,  # losses
            0.0,  # cop
            0.0,  # effective_cop
            0.0,  # mix_temp_input
            0.0,  # mix_temp_output
            0.0,  # avg_plr
            0.0)  # time_active
end

"""
A struct containing lots of temporary values for internal calculations of the heat pump.
This is used to avoid excessively long argument and return lists and to improve readability.
Unless you're changing the code for heat pumps, this can be ignored entirely.
"""
mutable struct HPEnergies
    potential_el_in::Floathing
    potential_el_in_layered::Vector{<:Floathing}
    potentials_heat_in::Vector{<:Floathing}
    potentials_heat_out::Vector{<:Floathing}
    available_el_in::Float64
    in_uacs_el::Vector{<:Stringing}
    available_heat_in::Vector{<:Floathing}
    available_heat_out::Vector{<:Floathing}
    max_usage_fraction::Float64
    in_indices::Vector{Integer}
    in_temps_min::Vector{<:Temperature}
    in_temps_max::Vector{<:Temperature}
    in_uacs_heat::Vector{<:Stringing}
    out_indices::Vector{Integer}
    out_temps_min::Vector{<:Temperature}
    out_temps_max::Vector{<:Temperature}
    out_uacs::Vector{<:Stringing}
    heat_in_has_inf_energy::Bool
    heat_out_has_inf_energy::Bool
    double_to_single_idx::Dict{Tuple{Integer,Integer},Tuple{Integer,Integer}}
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
    slices_times::Vector{Floathing}
    used_plrs::Vector{Floathing}

    function HPEnergies()
        return new(0.0,
                   Vector{Floathing}(),
                   Vector{Floathing}(),
                   Vector{Floathing}(),
                   0.0,
                   Vector{Stringing}(),
                   Vector{Floathing}(),
                   Vector{Floathing}(),
                   0.0,
                   Vector{Integer}(),
                   Vector{Temperature}(),
                   Vector{Temperature}(),
                   Vector{Stringing}(),
                   Vector{Integer}(),
                   Vector{Temperature}(),
                   Vector{Temperature}(),
                   Vector{Stringing}(),
                   false,
                   false,
                   Dict{Tuple{Integer,Integer},Tuple{Integer,Integer}}(),
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
                   Vector{Stringing}(),
                   Vector{Floathing}(),
                   Vector{Floathing}())
    end
end

"""
Reset the temporary slices of an energies container.
"""
function reset_temp_slices!(energies::HPEnergies)::HPEnergies
    energies.slices_temp_el_in = Vector{Floathing}()
    energies.slices_temp_heat_in = Vector{Floathing}()
    energies.slices_temp_heat_in_temperature = Vector{Temperature}()
    energies.slices_temp_heat_in_uac = Vector{Stringing}()
    energies.slices_temp_heat_out = Vector{Floathing}()
    energies.slices_temp_heat_out_temperature = Vector{Temperature}()
    energies.slices_temp_heat_out_uac = Vector{Stringing}()
    energies.slices_temp_times = Vector{Floathing}()
    energies.double_to_single_idx = Dict{Tuple{Integer,Integer},Integer}()
    energies.used_plrs = Vector{Floathing}()
    return energies
end

"""
Add a temporary slice to the given energies container.
"""
function add_temp_slice!(energies::HPEnergies, el_in::Floathing, heat_in::Floathing, heat_out::Floathing,
                         heat_in_temp::Temperature, heat_out_temp::Temperature, heat_in_uac::Stringing,
                         heat_out_uac::Stringing, time::Floathing)
    push!(energies.slices_temp_el_in, el_in)
    push!(energies.slices_temp_heat_in, heat_in)
    push!(energies.slices_temp_heat_in_temperature, heat_in_temp)
    push!(energies.slices_temp_heat_in_uac, heat_in_uac)
    push!(energies.slices_temp_heat_out, heat_out)
    push!(energies.slices_temp_heat_out_temperature, heat_out_temp)
    push!(energies.slices_temp_heat_out_uac, heat_out_uac)
    push!(energies.slices_temp_times, time)
end

"""
Resets the available energies of an energies container based on the potentials.

If a specific input or output layer index is given, that layer will be set as infinite and
the calculation will be done with only that layer on the input or output side.

# Arguments
- `energies::HPEnergies`: The energies container
- `fixed_heat_in::Union{Nothing,Integer}`: If a layer index is given, sets that input layer
    to infinite and selects it exclusively for sources.
- `fixed_heat_out::Union{Nothing,Integer}`: If a layer index is given, sets that output
    layer to infinite and selects it exclusively for sinks.
# Returns
- `HPEnergies`: The energies container
"""
function reset_available!(energies::HPEnergies;
                          fixed_heat_in::Union{Nothing,Integer}=nothing,
                          fixed_heat_out::Union{Nothing,Integer}=nothing)::HPEnergies
    energies.available_el_in = copy(energies.potential_el_in)
    energies.available_heat_in = copy(energies.potentials_heat_in)
    energies.available_heat_out = copy(energies.potentials_heat_out)
    if fixed_heat_in !== nothing
        energies.available_heat_in[fixed_heat_in] = Inf
    elseif fixed_heat_out !== nothing
        energies.available_heat_out[fixed_heat_out] = Inf
    end
    energies.in_indices = fixed_heat_in !== nothing ? [fixed_heat_in] :
                          [Int(i) for i in eachindex(energies.available_heat_in)]
    energies.out_indices = fixed_heat_out !== nothing ? [fixed_heat_out] :
                           [Int(i) for i in eachindex(energies.available_heat_out)]
    return energies
end

"""
Adds the heat_in temp slices to the existing slices output fields.
"""
function add_heat_in_temp_to_slices(energies::HPEnergies)::HPEnergies
    append!(energies.slices_heat_in, energies.slices_temp_heat_in)
    append!(energies.slices_heat_in_temperature, energies.slices_temp_heat_in_temperature)
    append!(energies.slices_heat_in_uac, energies.slices_temp_heat_in_uac)
    return energies
end

"""
Adds the heat_out temp slices to the existing slices output fields.
"""
function add_heat_out_temp_to_slices(energies::HPEnergies)::HPEnergies
    append!(energies.slices_heat_out, energies.slices_temp_heat_out)
    append!(energies.slices_heat_out_temperature, energies.slices_temp_heat_out_temperature)
    append!(energies.slices_heat_out_uac, energies.slices_temp_heat_out_uac)
    return energies
end

"""
Copies the electricity and heat_in temp slices to the slices output fields and overwrites
existing values based on whether the sums over the temp slices are smaller than the existing
values.
"""
function copy_heat_in_temp_to_slices(energies::HPEnergies)::HPEnergies
    if length(energies.slices_heat_in) == 0
        energies.slices_el_in = copy(energies.slices_temp_el_in)
        energies.slices_heat_in = copy(energies.slices_temp_heat_in)
        energies.slices_heat_in_temperature = copy(energies.slices_temp_heat_in_temperature)
        energies.slices_heat_in_uac = copy(energies.slices_temp_heat_in_uac)
    else
        if _isless(_sum(energies.slices_el_in), _sum(energies.slices_temp_el_in))
            energies.slices_el_in = energies.slices_temp_el_in
        end
        if _isless(_sum(energies.slices_heat_in), _sum(energies.slices_temp_heat_in))
            energies.slices_heat_in = energies.slices_temp_heat_in
            energies.slices_heat_in_temperature = energies.slices_temp_heat_in_temperature
            energies.slices_heat_in_uac = energies.slices_temp_heat_in_uac
        end
    end
    return energies
end

"""
Copies the electricity and heat_out temp slices to the slices output fields and overwrites
existing values based on whether the sums over the temp slices are smaller than the existing
values.
"""
function copy_heat_out_temp_to_slices(energies::HPEnergies)::HPEnergies
    if length(energies.slices_heat_out) == 0
        energies.slices_el_in = copy(energies.slices_temp_el_in)
        energies.slices_heat_out = copy(energies.slices_temp_heat_out)
        energies.slices_heat_out_temperature = copy(energies.slices_temp_heat_out_temperature)
        energies.slices_heat_out_uac = copy(energies.slices_temp_heat_out_uac)
    else
        if _isless(_sum(energies.slices_el_in), _sum(energies.slices_temp_el_in))
            energies.slices_el_in = energies.slices_temp_el_in
        end
        if _isless(_sum(energies.slices_heat_out), _sum(energies.slices_temp_heat_out))
            energies.slices_heat_out = energies.slices_temp_heat_out
            energies.slices_heat_out_temperature = energies.slices_temp_heat_out_temperature
            energies.slices_heat_out_uac = energies.slices_temp_heat_out_uac
        end
    end
    return energies
end

"""
Copies the temporary slices of the given energy container to the actual slices.

This overwrites any existing actual slices.
"""
function copy_temp_to_slices!(energies::HPEnergies)::HPEnergies
    energies.slices_el_in = copy(energies.slices_temp_el_in)
    energies.slices_heat_in = copy(energies.slices_temp_heat_in)
    energies.slices_heat_in_temperature = copy(energies.slices_temp_heat_in_temperature)
    energies.slices_heat_in_uac = copy(energies.slices_temp_heat_in_uac)
    energies.slices_heat_out = copy(energies.slices_temp_heat_out)
    energies.slices_heat_out_temperature = copy(energies.slices_temp_heat_out_temperature)
    energies.slices_heat_out_uac = copy(energies.slices_temp_heat_out_uac)
    energies.slices_times = copy(energies.slices_temp_times)
    return energies
end

function initialise!(unit::HeatPump, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.input_interfaces[unit.m_heat_in],
                          unload_storages(unit.controller, unit.m_heat_in))
    set_storage_transfer!(unit.input_interfaces[unit.m_el_in],
                          unload_storages(unit.controller, unit.m_el_in))
    set_storage_transfer!(unit.output_interfaces[unit.m_heat_out],
                          load_storages(unit.controller, unit.m_heat_out))
    if unit.has_secondary_interface
        if unit.primary_el_sources == [] || unit.secondary_el_sources == []
            @error "In heat pump $(unit.uac), a secondary interface is requested. Provide both `primary_el_sources` and " *
                   "`secondary_el_sources` to specify which electricity source should be used for the secondary interface."
        end
        common = intersect(unit.primary_el_sources, unit.secondary_el_sources)
        if !isempty(common)
            @error "In heat pump $(unit.uac), a secondary interface is requested. The following uac(s) are both in " *
                   "`primary_el_sources` and `secondary_el_sources`, but they need to be unique: $(common)"
        end
        if unit.output_interfaces[unit.m_heat_out_secondary].target.sys_function !== sf_bus ||
           unit.output_interfaces[unit.m_heat_out].target.sys_function !== sf_bus ||
           unit.output_interfaces[unit.m_heat_out_secondary].target !== unit.output_interfaces[unit.m_heat_out].target
            @error "In heat pump $(unit.uac), a secondary interface is requested. Then, for both output interfaces, " *
                   "the same bus has to be connected, which is not the case!"
        end
        set_storage_transfer!(unit.output_interfaces[unit.m_heat_out_secondary],
                              load_storages(unit.controller, unit.m_heat_out_secondary))
    end
end

function control(unit::HeatPump, components::Grouping, sim_params::Dict{String,Any})
    update(unit.controller)
    # for fixed input/output temperatures, overwrite the interface with those. otherwise
    # highest will choose the interface's temperature (including nothing)
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

function set_max_energies!(unit::HeatPump,
                           el_in::Union{Floathing,Vector{<:Floathing}},
                           heat_in::Union{Floathing,Vector{<:Floathing}},
                           heat_out::Union{Floathing,Vector{<:Floathing}},
                           heat_out_secondary::Union{Floathing,Vector{<:Floathing}},
                           slices_heat_in_temperature::Union{Temperature,Vector{<:Temperature}}=nothing,
                           slices_heat_out_temperature::Union{Temperature,Vector{<:Temperature}}=nothing,
                           slices_heat_out_secondary_temperature::Union{Temperature,Vector{<:Temperature}}=nothing,
                           purpose_uac_heat_in::Union{Stringing,Vector{Stringing}}=nothing,
                           purpose_uac_heat_out::Union{Stringing,Vector{Stringing}}=nothing,
                           purpose_uac_heat_out_secondary::Union{<:Stringing,Vector{<:Stringing}}=nothing,
                           has_calculated_all_maxima_heat_in::Bool=false,
                           has_calculated_all_maxima_heat_out::Bool=false)
    set_max_energy!(unit.input_interfaces[unit.m_el_in], el_in)
    set_max_energy!(unit.input_interfaces[unit.m_heat_in], heat_in, slices_heat_in_temperature, nothing,
                    purpose_uac_heat_in, has_calculated_all_maxima_heat_in)
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], heat_out, nothing, slices_heat_out_temperature,
                    purpose_uac_heat_out, has_calculated_all_maxima_heat_out)
    if unit.has_secondary_interface
        set_max_energy!(unit.output_interfaces[unit.m_heat_out_secondary], heat_out_secondary, nothing,
                        slices_heat_out_secondary_temperature, purpose_uac_heat_out_secondary,
                        has_calculated_all_maxima_heat_out)
    end
end

"""
Determines the temperature of the input/output layer and if it should be skipped.

Layers are skipped if they have a temperature, but it does fall into the temperature band
defined by the heat pump input/output temperatures (if any).

# Arguments
- `unit::HeatPump`: The heat pump.
- `current_idx::Integer`: The index of the layer.
- `temps_min::Vector{<:Temperature}`: The minimum temperatures of all layers.
- `temps_max::Vector{<:Temperature}`: The maximum temperatures of all layers.
- `input::Bool`: If the layer is an input or not. Defaults to true.
# Returns
- `Bool`: If the layer should be skipped.
- `Temperature`: The layer temperature.
"""
function get_layer_temperature(unit::HeatPump,
                               current_idx::Integer,
                               temps_min::Vector{<:Temperature},
                               temps_max::Vector{<:Temperature};
                               input::Bool=true)::Tuple{Bool,Temperature}
    static_temperature = input ? unit.input_temperature : unit.output_temperature
    if static_temperature === nothing
        layer_temp = highest(temps_min[current_idx], temps_max[current_idx])
        if layer_temp === nothing && unit.constant_cop === nothing
            term = input ? "input" : "output"
            @error "Error: The $(term) temperature for $(unit.uac) could not be detected. " *
                   "Please specify one with the parameter '$(term)_temperature' or check " *
                   "the connected components."
            throw(InputError())
        end
        return false, layer_temp
    end

    # skip slice if given static temperature is not within the temperature band of the layer
    if (temps_min[current_idx] !== nothing && temps_min[current_idx] > static_temperature
        ||
        temps_max[current_idx] !== nothing && temps_max[current_idx] < static_temperature)
        # end of condition
        return true, nothing
    else
        return false, static_temperature
    end
end

"""
Calculates the new COP due to icing losses.

# Arguments
`unit::HeatPump`: The heat pump.
`cop::Floathing`: The COP without icing losses.
`in_temp::Temperature`: The input temperature of the current layer.
# Returns
- `Floathing`: The COP corrected due to icing losses. Returns nothing if the input COP is
    nothing.
"""
function icing_correction(unit::HeatPump, cop::Floathing, in_temp::Temperature)::Floathing
    if cop === nothing
        return nothing
    end

    lin_factor = unit.icing_coefficients[1] + unit.icing_coefficients[2] * in_temp
    exp_factor = unit.icing_coefficients[3] * exp(-1.0
                                                  * (in_temp - unit.icing_coefficients[4])^2
                                                  /
                                                  unit.icing_coefficients[5])
    return cop * (1 - 0.01 * (max(0.0, lin_factor) + exp_factor))
end

"""
    filter_by_transformer(energies, sim_param; heat_in=true)

Filters the heat input or heat output potentials by system function so only transformers are
listed.

# Args
- `energie::HPEnergies`: The energies container
- `sim_params::Dict{String,Any}`: Simulation parameters
- `heat_in::Bool`: (Optional) If true, filters the inputs. Otherwise the outputs. Defaults
    to true.
# Returns
- `Vector{Floathing}`: The filtered input or output potentials.
"""
function filter_by_transformer(energies::HPEnergies,
                               sim_params::Dict{String,Any};
                               heat_in::Bool=true)::Vector{Floathing}
    components = get_run(sim_params["run_ID"]).components
    filtered = []
    for (idx, uac) in enumerate(heat_in ? energies.in_uacs_heat : energies.out_uacs)
        if uac !== nothing && components[uac].sys_function === sf_transformer
            if heat_in
                push!(filtered, energies.potentials_heat_in[idx])
            else
                push!(filtered, energies.potentials_heat_out[idx])
            end
        end
    end
    return filtered
end

"""
Evaluates the results of the slicing algorithm during the optimisation of PLRs.

This serves as the objective function to be minimised and includes the following factors:
* Barrier term to exclude implausible PLRs
* Meeting demands exactly (higher and lower sums are evaluated less favourably)
* For inverter heat pumps: Minismising electricity input
* For on-off heat pumps: Using up the entire time step

# Arguments
- `energies::HPEnergies`: The energies container
- `unit::HeatPump`: The heat pump
- `plrs::Vector{Float64}`: The chosen PLRs
- `sim_params::Dict{String,Any}`: Simulation parameters
# Returns
- `Float64`: The evaluation of the results. A lower number is better.
"""
function evaluate(energies::HPEnergies, unit::HeatPump, plrs::Vector{Float64}, sim_params::Dict{String,Any})::Float64
    # barrier term if implausible PLRs are given
    for plr in plrs
        if plr > 1.0 || plr < 0.0
            return 1000000
        end
    end

    # heat_out term with a (usually) higher weight factor so meeting demands exactly is
    # preferred by the optimisation
    heat_out_sum = sum(energies.slices_temp_heat_out; init=0.0)
    demand_sum = sum(energies.potentials_heat_out; init=0.0)
    heat_part = abs(demand_sum - heat_out_sum) / demand_sum

    if unit.model_type == "on-off"
        # for on-off heat pumps we "optimise" the time to use up the whole timestep. this is not
        # actually better in terms of efficiency, but since the PLF function is supposed to
        # capture the effects of cycling, we want to use as much time as available
        time_sum = sum(energies.slices_temp_times; init=0.0)
        time_part = (sim_params["time_step_seconds"] - time_sum) / sim_params["time_step_seconds"]
        val = unit.eval_factor_heat * heat_part + unit.eval_factor_time * time_part
    else
        # for inverter-driven heat pumps we perform actualy optimisation by adjusting PLRs
        # to optimise for lowest electricity use while maintaining meeting demands
        elec_sum = sum(energies.slices_temp_el_in; init=0.0)
        elec_part = elec_sum / heat_out_sum
        val = unit.eval_factor_heat * heat_part + unit.eval_factor_elec * elec_part
    end

    return val
end

"""
Calculate the energies for the given slice.

This is essentially how the heat pump would be calculated in the 1:1 case. The COP is
determined first, then the available energies on all three inputs/outputs are checked as
limiting factors. Bypass operation and icing losses are also considered. The function
specifically does not determine the PLR. For the 1:1 case this could be done here too, but
for the overarching slicing algorithm the PLR must be a required input.

# Arguments
- `unit::HeatPump`: The heat pump.
- `available_el_in::Float64`: Available electricity.
- `available_heat_in::Floathing`: Available heat input.
- `available_heat_out::Floathing`: Available heat output.
- `in_temp::Temperature`: Input temperature.
- `out_temp::Temperature`: Output temperature.
- `plr::Float64`: The PLR of the slice.
# Returns
- `Floathing`: Used input heat.
- `Floathing`: Used electricity.
- `Floathing`: Used output heat.
- `Temperature`: Temperature of input heat. Same as input argument.
- `Temperature`: Temperature of output heat. Same as output argument unless in bypass
    operation, when it's equal to the input temperature.
"""
function handle_slice(unit::HeatPump,
                      available_el_in::Float64,
                      available_heat_in::Floathing,
                      available_heat_out::Floathing,
                      in_temp::Temperature,
                      out_temp::Temperature,
                      plr::Float64)::Tuple{Floathing,Floathing,Floathing,Temperature,Temperature,Float64}
    # determine COP depending on three cases. a constant COP precludes the use of a bypass
    if unit.constant_cop !== nothing
        cop = unit.constant_cop * unit.plf_function(plr)
    elseif in_temp >= out_temp
        cop = unit.bypass_cop
    else
        cop = unit.dynamic_cop(in_temp, out_temp)
        if cop === nothing
            @error ("Input and/or output temperature for heatpump $(unit.uac) is not " *
                    "given. Provide temperatures or fixed cop.")
            throw(InputError())
        end
        cop *= unit.plf_function(plr)
        if unit.consider_icing
            cop = icing_correction(unit, cop, in_temp)
        end
    end

    if cop < 1.0
        @warn ("Calculated COP of heat pump $(unit.uac) was below 1.0. Please check the " *
               "input for mistakes as this should not happen. COP was set from $(round(cop;digits=2)) to 1.0")
        cop = 1.0
    end

    # calculate energies with the current cop
    # energies for current slice with potential heat in as basis
    used_heat_in = available_heat_in
    if cop == 1.0
        used_el_in = available_el_in
    else
        used_el_in = used_heat_in / (cop - 1.0)
    end
    used_heat_out = used_heat_in + used_el_in

    # check heat out as limiter
    if used_heat_out > available_heat_out
        used_heat_out = available_heat_out
        used_el_in = used_heat_out / cop
        used_heat_in = used_el_in * (cop - 1.0)
    end

    # check electricity in as limiter
    if used_el_in > available_el_in
        used_el_in = available_el_in
        used_heat_in = used_el_in * (cop - 1.0)
        used_heat_out = used_heat_in + used_el_in
    end

    return (used_heat_in,
            used_el_in,
            used_heat_out,
            in_temp,
            out_temp,
            cop)
end

"""
Performs the slicing algorithm given a set list of PLRs for each possible slice.
    
A slice is a unique combination of input and output layer.

Note the difference in indexing between the PLRs, which are listed in a particular order for
all possible slices, and the return values, where only selected slices are listed. For this
reason the energies container has a field `double_to_single_idx`, which saves the mapping
between these indices.

# Arguments
- `unit::HeatPump`: The heat pump.
- `sim_params::Dict{String,Any}`: Simulation parameters.
- `energies::HPEnergies`: The energies container.
- `plrs::Array{Float64}`: The PLR for each slice.
- `is_final::Bool`: If this is the final calculation of the optimisation process
- `fixed_heat_in`: If given, this input layer is fixed (see reset_available!())
- `fixed_heat_out`:  If given, this output layer is fixed (see reset_available!())
# Returns
- `HPEnergies`: The energies container.
"""
function calculate_slices(unit::HeatPump,
                          sim_params::Dict{String,Any},
                          energies::HPEnergies,
                          plrs::Vector{Float64},
                          is_final::Bool,
                          fixed_heat_in::Union{Nothing,Integer}=nothing,
                          fixed_heat_out::Union{Nothing,Integer}=nothing)::HPEnergies
    # reset at the beginning, because the optimisation algorithm will call this function
    # multiple times with different plrs (and also at least once in both the potential and
    # process step)
    energies = reset_available!(energies; fixed_heat_in, fixed_heat_out)
    energies = reset_temp_slices!(energies)

    # transfer the PLRs, which are a flat vector due to requirements of the optimisation
    # algorithm, to the index dict used in the slicing algorithm. maybe we can find a better
    # solution than to use three different index systems, but at least it works
    idx = 1
    for (src_idx, _) in enumerate(energies.potentials_heat_in)
        for (snk_idx, _) in enumerate(energies.potentials_heat_out)
            # barrier short-circuit. if we get implausible PLR values we can skip the
            # slicing and return as is (the barrier term in the evaluate function will
            # catch this case)
            if plrs[idx] < 0.0 || plrs[idx] > 1.0
                return energies
            end

            # the first index is the index in the PLR list
            # the slice index (second) is set later, here's just a dummy value
            energies.double_to_single_idx[(src_idx, snk_idx)] = (idx, -1)
            idx += 1
        end
    end

    # init used PLRs as zero so unused slices don't contribute to the average
    energies.used_plrs = fill(0.0, length(plrs))

    slice_idx::Int = 1
    src_idx::Int = 1
    snk_idx::Int = 1
    available_time = sim_params["time_step_seconds"]
    EPS = sim_params["epsilon"]

    # keep running until we ran out of either sources or sinks to check
    while (src_idx <= length(energies.available_heat_in) && snk_idx <= length(energies.available_heat_out))
        # we can skip calculation if there is no energy left to be distributed or we ran out
        # of available time. this can happen for example during potential calculations when
        # another transformer has already used up all available energies. Do not consider
        # heat_in in process to account for COP == 1.0.
        if energies.available_el_in < EPS ||
           ((unit.cop == 0.0 || unit.cop > 1.0) && sum(energies.available_heat_in; init=0.0) < EPS) ||
           sum(energies.available_heat_out; init=0.0) < EPS ||
           available_time < EPS
            # end of condition
            break
        end

        # apply restrictions of control modules for a slice
        if !check_src_to_snk(unit.controller, energies.in_uacs_heat[src_idx], energies.out_uacs[snk_idx])
            snk_idx += 1
            continue
        end

        # check and determine input temperature of layer, also skip if it's not in the list
        # of indices to be used (this is used by the mechanism for transformer chains and
        # is controlled by arguments fixed_heat_in and fixed_heat_out)
        skip_slice, src_temperature = get_layer_temperature(unit, src_idx,
                                                            energies.in_temps_min,
                                                            energies.in_temps_max;
                                                            input=true)
        if skip_slice || !(src_idx in energies.in_indices)
            src_idx += 1
            continue
        end

        # same for output
        skip_slice, snk_temperature = get_layer_temperature(unit, snk_idx,
                                                            energies.out_temps_min,
                                                            energies.out_temps_max;
                                                            input=false)
        if skip_slice || !(snk_idx in energies.out_indices)
            snk_idx += 1
            continue
        end

        # there are five factors that decide the heat output of the heat pump:
        # 1. the minimum power for the given temperatures
        # 2. the maximum power for the given temperatures
        # 3. the maximum usage fraction (as fraction of the nominal power)
        # 4. the PLR of the slice
        # 5. how much time is left in the time step
        plr_idx = energies.double_to_single_idx[(src_idx, snk_idx)][1]
        min_power_frac = max(0.0, min(1.0, unit.min_power_function(src_temperature, snk_temperature)))
        min_power = unit.design_power_th * min_power_frac
        max_power_frac = max(0.0, min(1.0, unit.max_power_function(src_temperature, snk_temperature)))
        max_power = unit.design_power_th * max_power_frac
        plr_power = min(max_power * plrs[plr_idx], max_power * energies.max_usage_fraction)
        used_power = max(min_power, plr_power)

        # the fudge factor causes us to slightly overestimate the maximum heat out energy
        # assuming the remaining time is used up. this helps with the solver meeting demands
        # exactly and the resulting "error" is much smaller than the typical overall
        # inaccuracies of the calculation. however we only do this in the final calculation
        available_heat_out = min(energies.available_heat_out[snk_idx],
                                 (is_final ? unit.fudge_factor : 1.0)
                                 * available_time * used_power / 3600)

        used_heat_in,
        used_el_in,
        used_heat_out,
        src_temperature,
        snk_temperature,
        cop = handle_slice(unit,
                           energies.available_el_in,
                           energies.available_heat_in[src_idx],
                           available_heat_out,
                           src_temperature,
                           snk_temperature,
                           plrs[plr_idx])

        used_time = used_heat_out * 3600 / used_power
        energies.used_plrs[plr_idx] = used_heat_out / sim_params["watt_to_wh"](max_power)

        # finally all checks done, we add the slice and update remaining energies
        add_temp_slice!(energies, used_el_in, used_heat_in, used_heat_out,
                        src_temperature, snk_temperature,
                        energies.in_uacs_heat[src_idx],
                        energies.out_uacs[snk_idx],
                        used_time)

        energies.double_to_single_idx[(src_idx, snk_idx)] = (plr_idx, slice_idx)
        slice_idx += 1

        available_time -= used_time
        energies.available_el_in -= used_el_in
        energies.available_heat_in[src_idx] -= used_heat_in
        energies.available_heat_out[snk_idx] -= used_heat_out

        # go to the next slice if input and/or output has been used up. if instead both
        # layers have energy left, this was caused by the PLR of the slice being too low or
        # the heat pump being undersized. in either case we need to advance both indices so
        # we don't recalculate the same slice twice
        if (energies.available_heat_in[src_idx] >= EPS && cop !== 1.0) &&
           energies.available_heat_out[snk_idx] >= EPS
            # end of condition
            src_idx += 1
            snk_idx += 1
        else
            if energies.available_heat_in[src_idx] < EPS && cop !== 1.0
                src_idx += 1
            end
            if energies.available_heat_out[snk_idx] < EPS
                snk_idx += 1
            end
        end
    end

    return energies
end

"""
Finds the best slicing for the given potentials and demands.

Depending on the model type this requires an optimisation process that finds the best PLRs
for each possible slice. In this case the function `evaluate` determines the optimum that is
found. The optimisation can be configured and technically is not guarranteed to find the
global optimum, however the default parameters should lead to decent results.

For the model type `simplified` the optimisation is not required and this function simply
runs the slicing algorithm once assuming PLRs of 1.0 for all slices.

# Arguments
- `unit::HeatPump`: The heat pump.
- `sim_params::Dict{String,Any}`: Simulation parameters.
- `energies::HPEnergies`: The energies container.
- `fixed_heat_in`: If given, this input layer is fixed (see reset_available!())
- `fixed_heat_out`:  If given, this output layer is fixed (see reset_available!())
# Returns
- `HPEnergies`: The energies container.
"""
function find_best_slicing(unit::HeatPump,
                           sim_params::Dict{String,Any},
                           energies::HPEnergies;
                           fixed_heat_in::Union{Nothing,Integer}=nothing,
                           fixed_heat_out::Union{Nothing,Integer}=nothing)::HPEnergies
    # no optimisation for simplified heat pump, as the slicing algorithm will result in a
    # solution independent of chosen PLRs assuming sufficient time is available to meet
    # demands, hence the PLRs are set to 1.0
    if unit.model_type == "simplified"
        default_plrs = fill(1.0,
                            length(energies.potentials_heat_in) *
                            length(energies.potentials_heat_out))
        # technically this is the only and final calculation, but the fudge factor is only
        # relevant for optimisation, which is why we disable it here
        energies = calculate_slices(unit, sim_params, energies, default_plrs, false, fixed_heat_in, fixed_heat_out)
        return energies
    end

    # estimate initial PLRs by demand sum over nominal power
    heat_out_sum = sum(energies.potentials_heat_out; init=0.0)
    plr = heat_out_sum / sim_params["watt_to_wh"](unit.design_power_th)
    plr = max(min(1.0, plr), 0.0)
    initial_plrs = fill(plr,
                        length(energies.potentials_heat_in) *
                        length(energies.potentials_heat_out))
    lower_plrs = fill(0.0,
                      length(energies.potentials_heat_in) *
                      length(energies.potentials_heat_out))
    upper_plrs = fill(1.0,
                      length(energies.potentials_heat_in) *
                      length(energies.potentials_heat_out))

    # run optimisation to find PLRs that meet demands and are optimal by criteria depending
    # on the model type (see function evaluate)
    results = optimize(plrs -> evaluate(calculate_slices(unit,
                                                         sim_params, energies, plrs, false,
                                                         fixed_heat_in, fixed_heat_out), unit, plrs, sim_params),
                       lower_plrs, upper_plrs, initial_plrs, NelderMead(),
                       Options(; iterations=Int64(unit.nr_optimisation_passes),
                               x_abstol=unit.x_abstol,
                               f_abstol=unit.f_abstol))
    optimal_plrs = minimizer(results)
    energies = calculate_slices(unit, sim_params, energies, optimal_plrs, true, fixed_heat_in, fixed_heat_out)

    return energies
end

"""
Calculates the energies the heat pump can process in the timestep.

This is split in an inner function (see find_best_slicing) and this outer function, because
the calculation happens differently depending on whether there are other transformers in the
input exchanges, outputs, both or neither. Also the calculation of losses requires another
layer around the function for calculating the slices.

# Arguments
- `unit::HeatPump`: The heat pump.
- `sim_params::Dict{String,Any}`: Simulation parameters.
# Returns
- `HPEnergies`: The energies container containing the results of the calculations.
"""
function calculate_energies(unit::HeatPump, sim_params::Dict{String,Any})
    energies = HPEnergies()
    do_calculation = true

    # get electricity potential and reduce it by constant power draw (or however much
    # is available)
    energies.potential_el_in_layered, energies.in_uacs_el = check_el_in_layered(unit, sim_params)
    energies.potential_el_in = sum(energies.potential_el_in_layered)
    unit.current_constant_loss = min(energies.potential_el_in,
                                     unit.constant_loss_energy)
    energies.potential_el_in -= unit.current_constant_loss

    # shortcut if we're limited by zero input electricity
    if energies.potential_el_in <= 0.0
        do_calculation = false
    end

    # get max usage fraction from control modules and shortcut if it is zero
    energies.max_usage_fraction = upper_plr_limit(unit.controller, sim_params)
    if energies.max_usage_fraction <= 0.0
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
        # heat pump, depending on wether or not any input heat or output heat transformer has a
        # value of infinite as the potential. if this is the case for both the input and output,
        # the calculation cannot be resolved.
        energies.heat_in_has_inf_energy = any(isinf, filter_by_transformer(energies, sim_params; heat_in=true))
        energies.heat_out_has_inf_energy = any(isinf, filter_by_transformer(energies, sim_params; heat_in=false))

        if energies.heat_in_has_inf_energy && energies.heat_out_has_inf_energy
            @warn "The heat pump $(unit.uac) has unknown energies in both its inputs and " *
                  "outputs. This cannot be resolved. Please check the order of operation and " *
                  "make sure that either the inputs or the outputs have been fully calculated " *
                  "before the heat pump $(unit.uac) has its potential step."
            return energies
        end

        if energies.heat_in_has_inf_energy
            for heat_in_idx in eachindex(energies.potentials_heat_in)
                energies = find_best_slicing(unit, sim_params, energies; fixed_heat_in=heat_in_idx)
                energies = add_heat_in_temp_to_slices(energies)
                energies = copy_heat_out_temp_to_slices(energies)
            end
        elseif energies.heat_out_has_inf_energy
            for heat_out_idx in eachindex(energies.potentials_heat_out)
                energies = find_best_slicing(unit, sim_params, energies; fixed_heat_out=heat_out_idx)
                energies = add_heat_out_temp_to_slices(energies)
                energies = copy_heat_in_temp_to_slices(energies)
            end
        else
            energies = find_best_slicing(unit, sim_params, energies)
            energies = copy_temp_to_slices!(energies)
        end
    end

    # calculate average PLR, from active slices. because they already have been weighted
    # by max power for each slice, we can simply add them up here
    unit.avg_plr = 0.0
    for (_, (plr_idx, slice_idx)) in pairs(energies.double_to_single_idx)
        if slice_idx == -1
            continue
        end
        unit.avg_plr += energies.used_plrs[plr_idx]
    end

    # if minimum usage fraction was not reached, we discard all slices
    if unit.avg_plr < unit.min_usage_fraction
        energies.slices_el_in = Vector{Floathing}()
        energies.slices_heat_in = Vector{Floathing}()
        energies.slices_heat_out = Vector{Floathing}()
        energies.slices_times = Vector{Floathing}()
    end

    # calculate COP before losses
    el_in = sum(energies.slices_el_in; init=0.0)
    heat_in = sum(energies.slices_heat_in; init=0.0)
    heat_out = sum(energies.slices_heat_out; init=0.0)
    if el_in > sim_params["epsilon"]
        unit.cop = heat_out / el_in
    end

    # now set losses of the heat pump and add the losses to the actually consumed
    # power / heat for the slices
    unit.losses_power = -1.0 * el_in
    unit.losses_heat = -1.0 * heat_in
    energies.slices_el_in ./= unit.power_losses_factor
    energies.slices_heat_in ./= unit.heat_losses_factor

    el_in = sum(energies.slices_el_in; init=0.0)
    heat_in = sum(energies.slices_heat_in; init=0.0)
    unit.losses_power += el_in
    unit.losses_heat += heat_in

    # constant losses are always incurred, so we might need to add a slice. if there are
    # existing slices, we add the constant loss averaged to the slices.
    if length(energies.slices_el_in) == 0 && unit.current_constant_loss > 0.0
        push!(energies.slices_el_in, unit.current_constant_loss)
    else
        energies.slices_el_in .+= (unit.current_constant_loss / length(energies.slices_el_in))
    end

    return energies
end

function potential(unit::HeatPump, sim_params::Dict{String,Any})
    energies = calculate_energies(unit, sim_params)

    if sum(energies.slices_heat_out; init=0.0) < sim_params["epsilon"]
        set_max_energies!(unit, sum(energies.slices_el_in; init=0.0), 0.0, 0.0, 0.0)
        return
    end

    if unit.has_secondary_interface
        # split slices regarding electricity source
        primary, secondary = split_slices_primary_secondary(unit,
                                                            energies.potential_el_in_layered,
                                                            energies.in_uacs_el,
                                                            unit.primary_el_sources,
                                                            unit.secondary_el_sources,
                                                            energies.slices_el_in,
                                                            energies.slices_heat_out,
                                                            energies.slices_heat_out_temperature,
                                                            energies.slices_heat_out_uac,
                                                            energies.heat_out_has_inf_energy,
                                                            sim_params["epsilon"])

        set_max_energies!(unit,
                          energies.slices_el_in,
                          energies.slices_heat_in,
                          primary.slices_heat_out,
                          secondary.slices_heat_out,
                          energies.slices_heat_in_temperature,
                          primary.slices_heat_out_temperature,
                          secondary.slices_heat_out_temperature,
                          energies.slices_heat_in_uac,
                          primary.slices_heat_out_uac,
                          secondary.slices_heat_out_uac,
                          energies.heat_in_has_inf_energy,
                          energies.heat_out_has_inf_energy)
    else
        set_max_energies!(unit,
                          energies.slices_el_in,
                          energies.slices_heat_in,
                          energies.slices_heat_out,
                          Float64[],
                          energies.slices_heat_in_temperature,
                          energies.slices_heat_out_temperature,
                          Float64[],
                          energies.slices_heat_in_uac,
                          energies.slices_heat_out_uac,
                          String[],
                          energies.heat_in_has_inf_energy,
                          energies.heat_out_has_inf_energy)
    end
end

"""
    split_slices_primary_secondary(unit, ...)

Splits the heat pump output slices into "primary" and "secondary" groups based on the user-defined mapping of 
the electricity sources to primary/secondary output interfaces.

# Arguments
- `unit::HeatPump`: The heat pump unit.
- `potential_el_in_layered::Vector{Float64}`: Available electricity input as layer.
- `in_uacs_el::Vector{String}`: UACs for each electricity input layer.
- `primary_uac::Vector{String}`: UACs considered as "primary" sources.
- `secondary_uac::Vector{String}`: UACs considered as "secondary" sources.
- `slices_el_in::Vector{Floathing}`: Electricity used per output slice.
- `slices_heat_out::Vector{Floathing}`: Heat output per slice.
- `slices_heat_out_temperature::Vector{Floathing}`: Output temperature per slice.
- `slices_heat_out_uac::Vector{Stringing}`: UACs for each heat output slice.
- `heat_out_has_inf_energy::Bool`: Whether any output slice has infinite energy.
- `epsilon::Float64`: Numerical tolerance for slice splitting.

# Returns
- `primary::NamedTuple`: Slices attributed to primary sources for the primary output interface, with fields:
    - `slices_heat_out::Vector{Floathing}`
    - `slices_heat_out_temperature::Vector{Floathing}`
    - `slices_heat_out_uac::Vector{Stringing}`
    - `source_uac_el::Vector{Stringing}`
- `secondary::NamedTuple`: Slices attributed to secondary sources for the secondary output interface, with fields:
    - `slices_heat_out::Vector{Floathing}`
    - `slices_heat_out_temperature::Vector{Floathing}`
    - `slices_heat_out_uac::Vector{Stringing}`
    - `source_uac_el::Vector{Stringing}`
"""
function split_slices_primary_secondary(unit::HeatPump,
                                        potential_el_in_layered::Vector{Float64},
                                        in_uacs_el::Vector{String},
                                        primary_uac::Vector{String},
                                        secondary_uac::Vector{String},
                                        slices_el_in::Vector{Floathing},
                                        slices_heat_out::Vector{Floathing},
                                        slices_heat_out_temperature::Vector{Floathing},
                                        slices_heat_out_uac::Vector{Stringing},
                                        heat_out_has_inf_energy::Bool,
                                        epsilon::Float64)
    # convenience
    n_inputs = length(potential_el_in_layered)

    # --- outputs --------------------------------------------------------------
    T = eltype(slices_el_in)
    primary_heat_out = T[]
    primary_heat_out_temperature = T[]
    primary_heat_out_uac_out = Stringing[]
    primary_el_source_uac = Stringing[]   # from in_uacs_el

    secondary_heat_out = T[]
    secondary_heat_out_temperature = T[]
    secondary_heat_out_uac_out = Stringing[]
    secondary_el_source_uac = Stringing[]   # from in_uacs_el

    # --- main loop -----------------------------------------------------------
    layer = 1
    remaining = layer <= n_inputs ? potential_el_in_layered[layer] : 0.0

    for j in eachindex(slices_heat_out)
        if heat_out_has_inf_energy
            # for each heat_out entry, start from the beginning of the layers
            # here, only one slice_el_in is given, but it has to serve multiple slices_heat_out
            slice_total = slices_el_in[1]
            layer = 1
            remaining = layer <= n_inputs ? potential_el_in_layered[layer] : 0.0
        else
            slice_total = slices_el_in[j]
        end

        # skip zero slices to avoid division by zero
        if slice_total <= epsilon
            continue
        end

        slice_left = slice_total

        while slice_left > epsilon
            # if layer is exhausted, move to next
            while layer <= n_inputs && remaining ≤ epsilon
                layer += 1
                if layer <= n_inputs
                    remaining = potential_el_in_layered[layer]
                end
            end

            # take as much as we can from current layer
            take = min(slice_left, remaining)
            ratio = take / slice_total         # in (0,1]

            src = in_uacs_el[layer]
            is_primary = src in primary_uac
            is_secondary = src in secondary_uac

            if is_primary && is_secondary
                @error "In heat pump $(unit.uac), the source '$src' is in neither primary_uac nor secondary_uac."
            end

            # scale energy quantities
            heat_out_piece = slices_heat_out[j] * ratio

            # metadata copied 1:1
            Tout = slices_heat_out_temperature[j]
            uac_out = slices_heat_out_uac[j]

            if is_primary
                push!(primary_heat_out, heat_out_piece)
                push!(primary_heat_out_temperature, Tout)
                push!(primary_heat_out_uac_out, uac_out)
                push!(primary_el_source_uac, src)
            else
                push!(secondary_heat_out, heat_out_piece)
                push!(secondary_heat_out_temperature, Tout)
                push!(secondary_heat_out_uac_out, uac_out)
                push!(secondary_el_source_uac, src)
            end

            # update remaining and slice_left
            slice_left -= take
            remaining -= take
        end
    end

    if isempty(primary_heat_out)
        push!(primary_heat_out, 0.0)
        push!(primary_heat_out_temperature, nothing)
        push!(primary_heat_out_uac_out, nothing)
        push!(primary_el_source_uac, nothing)
    end
    primary = (slices_heat_out=primary_heat_out,
               slices_heat_out_temperature=primary_heat_out_temperature,
               slices_heat_out_uac=primary_heat_out_uac_out,
               source_uac_el=primary_el_source_uac)

    if isempty(secondary_heat_out)
        push!(secondary_heat_out, 0.0)
        push!(secondary_heat_out_temperature, nothing)
        push!(secondary_heat_out_uac_out, nothing)
        push!(secondary_el_source_uac, nothing)
    end
    secondary = (slices_heat_out=secondary_heat_out,
                 slices_heat_out_temperature=secondary_heat_out_temperature,
                 slices_heat_out_uac=secondary_heat_out_uac_out,
                 source_uac_el=secondary_el_source_uac)

    return primary, secondary
end

function process(unit::HeatPump, sim_params::Dict{String,Any})
    energies = calculate_energies(unit, sim_params)

    el_in = sum(energies.slices_el_in; init=0.0)
    heat_out = sum(energies.slices_heat_out; init=0.0)

    if heat_out < sim_params["epsilon"]
        # due to constant losses we are guaranteed to have an electricity slice, though
        # it might be zero. we also need to set the max_energy values to zero for the
        # heat input and output, as the sub! and add! methods do that when called
        set_max_energies!(unit, el_in, 0.0, 0.0, 0.0)
        sub!(unit.input_interfaces[unit.m_el_in], el_in)
    else
        sub!(unit.input_interfaces[unit.m_el_in], el_in)
        sub!(unit.input_interfaces[unit.m_heat_in],
             energies.slices_heat_in,
             energies.slices_heat_in_temperature,
             [nothing for _ in energies.slices_heat_in_temperature],
             energies.slices_heat_in_uac)

        if unit.has_secondary_interface
            # split slices regarding electricity source
            primary, secondary = split_slices_primary_secondary(unit,
                                                                energies.potential_el_in_layered,
                                                                energies.in_uacs_el,
                                                                unit.primary_el_sources,
                                                                unit.secondary_el_sources,
                                                                energies.slices_el_in,
                                                                energies.slices_heat_out,
                                                                energies.slices_heat_out_temperature,
                                                                energies.slices_heat_out_uac,
                                                                energies.heat_out_has_inf_energy,
                                                                sim_params["epsilon"])

            add!(unit.output_interfaces[unit.m_heat_out],
                 primary.slices_heat_out,
                 [nothing for _ in primary.slices_heat_out_temperature],
                 primary.slices_heat_out_temperature,
                 primary.slices_heat_out_uac)
            add!(unit.output_interfaces[unit.m_heat_out_secondary],
                 secondary.slices_heat_out,
                 [nothing for _ in secondary.slices_heat_out_temperature],
                 secondary.slices_heat_out_temperature,
                 secondary.slices_heat_out_uac)
        else
            add!(unit.output_interfaces[unit.m_heat_out],
                 energies.slices_heat_out,
                 [nothing for _ in energies.slices_heat_out_temperature],
                 energies.slices_heat_out_temperature,
                 energies.slices_heat_out_uac)
        end
    end
    # calculate losses, effective COP and mixing temperatures with the final values of
    # processed energies
    unit.losses = unit.losses_power + unit.losses_heat + unit.current_constant_loss
    if el_in > sim_params["epsilon"]
        unit.effective_cop = heat_out / el_in
    end
    unit.mix_temp_input = _weighted_mean(energies.slices_heat_in_temperature, energies.slices_heat_in)
    unit.mix_temp_output = _weighted_mean(energies.slices_heat_out_temperature, energies.slices_heat_out)

    # calculate active time as fraction of the simulation time step
    unit.time_active = sum(energies.slices_times; init=0.0) / sim_params["time_step_seconds"]
end

function component_has_minimum_part_load(unit::HeatPump)
    return unit.min_usage_fraction > 0.0
end

# has its own reset function as here more parameters are present that need to be reset in
# every timestep
function reset(unit::HeatPump)
    invoke(reset, Tuple{Component}, unit)

    # reset other parameter
    unit.cop = 0.0
    unit.effective_cop = 0.0
    unit.mix_temp_input = 0.0
    unit.mix_temp_output = 0.0
    unit.current_constant_loss = 0.0
    unit.losses_heat = 0.0
    unit.losses_power = 0.0
    unit.avg_plr = 0.0
    unit.time_active = 0.0
end

function output_values(unit::HeatPump)::Vector{String}
    output_vals = [string(unit.m_el_in) * ":IN",
                   string(unit.m_heat_in) * ":IN",
                   string(unit.m_heat_out) * ":OUT"]
    if unit.has_secondary_interface
        push!(output_vals, string(unit.m_heat_out_secondary) * ":OUT")
    end
    append!(output_vals,
            ["COP",
             "Effective_COP",
             "Avg_PLR",
             "Time_active",
             "MixingTemperature_Input",
             "MixingTemperature_Output",
             "Losses_power",
             "Losses_heat",
             "LossesGains"])
    return output_vals
end

function output_value(unit::HeatPump, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "COP"
        return unit.cop
    elseif key.value_key == "Effective_COP"
        return unit.effective_cop
    elseif key.value_key == "Avg_PLR"
        return unit.avg_plr
    elseif key.value_key == "Time_active"
        return unit.time_active
    elseif key.value_key == "MixingTemperature_Input"
        return unit.mix_temp_input
    elseif key.value_key == "MixingTemperature_Output"
        return unit.mix_temp_output
    elseif key.value_key == "Losses_power"
        return -unit.losses_power
    elseif key.value_key == "Losses_heat"
        return -unit.losses_heat
    elseif key.value_key == "LossesGains"
        return -unit.losses
    end
    throw(KeyError(key.value_key))
end

export HeatPump
