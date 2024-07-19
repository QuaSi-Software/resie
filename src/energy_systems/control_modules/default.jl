"""
Default control module that all components have. If no module is specified in the component
config, this module is initialised with default-only parameters. Supports several mechanisms
of controlling operation, including:
  * Defining for any input or output if the energy supplied/requested over that interface is
    allowed to be used for un-/loading storage components
"""
mutable struct CM_Default <: ControlModule
    name::String
    parameters::Dict{String,Any}

    function CM_Default(
        parameters::Dict{String,Any},
        components::Grouping,
        sim_params::Dict{String,Any}
    )
        default_parameters=Dict{String,Any}(
            "name" => "default",
        )

        return new(
            "default", # name
            Base.merge(default_parameters, parameters), # parameters
        )
    end
end

# default constructor with no input
CM_Default() = CM_Default(Dict{String,Any}(), Grouping(), Dict{String,Any}())

function update(mod::CM_Default)
    # nothing to do for the default module
end