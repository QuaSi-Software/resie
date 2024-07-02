"""
Default control module that all components have. If no module is specified in the component
config, this module is initialised with default-only parameters. Supports several mechanisms
of controlling operation, including:
  * Setting specific inputs or outputs to be assumed infinite when calculating the energies
    that can be utilised in a time step
  * Defining for any input or output if the energy supplied/requested over that interface is
    allowed to be used for un-/loading storage components
"""
mutable struct CM_Default <: ControlModule
    name::String
    parameters::Dict{String,Any}
    linked_components::Dict{String,Component}
    required_components::Dict{String,Type}

    function CM_Default(parameters::Dict{String,Any}, components::Grouping)
        default_parameters=Dict{String,Any}(
            "name" => "default",
            "consider_m_el_in" => true,
            "consider_m_el_out" => true,
            "consider_m_gas_in" => true,
            "consider_m_fuel_in" => true,
            "consider_m_h2_out" => true,
            "consider_m_o2_out" => true,
            "consider_m_heat_out" => true,
            "consider_m_heat_ht_out" => true,
            "consider_m_heat_lt_out" => true,
            "consider_m_heat_in" => true
        )

        return new(
            "default", # name
            Base.merge(default_parameters, parameters), # parameters
            Dict{String,Component}(), # linked_components
            Dict{String,Type}() # required_components
        )
    end
end

# default constructor with no input
CM_Default() = CM_Default(Dict{String,Any}(), Grouping())