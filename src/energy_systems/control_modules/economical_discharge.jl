"""
Control module for controlling the discharge of a battery depending on if it has a suffcient
amount of charge and if a linked PV plant produces power below a certain threshold.
"""
mutable struct CM_EconomicalDischarge <: ControlModule
    name::String
    parameters::Dict{String,Any}
    state_machine::StateMachine

    function CM_EconomicalDischarge(parameters::Dict{String,Any},
                                    components::Grouping,
                                    sim_params::Dict{String,Any},
                                    unit_uac::String)
        default_parameters = Dict{String,Any}(
            "name" => "economical_discharge",
            "pv_threshold" => 1.0, # [Wh], a value of 1.0 essentially requires that the PV
            "min_charge" => 0.2,   # plant produces (almost) no power at all
            "discharge_limit" => 0.05,
            "pv_plant_uac" => nothing,
            "battery_uac" => nothing,
        )
        params = Base.merge(default_parameters, parameters)

        if !(params["pv_plant_uac"] !== nothing
             && params["pv_plant_uac"] in keys(components)
             && components[params["pv_plant_uac"]] isa PVPlant)
            @error "Required PV plant component for control module economical discharge not given"
        end
        params["pv_plant"] = components[params["pv_plant_uac"]]

        if !(params["battery_uac"] !== nothing
             && params["battery_uac"] in keys(components)
             && components[params["battery_uac"]] isa Battery)
            @error "Required battery component for control module economical discharge not given"
        end
        params["battery"] = components[params["battery_uac"]]

        state_machine = StateMachine(UInt(1),               # state
                                     Dict{UInt,String}(     # state_names
                                         1 => "Load",
                                         2 => "Discharge",
                                     ),
                                     Dict{UInt,TruthTable}( # transitions
                                         1 => TruthTable(;
                                                         conditions=[function (state_machine)
                                                                         return params["pv_plant"].supply <
                                                                                params["pv_threshold"]
                                                                     end,
                                                                     function (state_machine)
                                                                         return params["battery"].load >=
                                                                                params["battery"].capacity *
                                                                                params["min_charge"]
                                                                     end],
                                                         table_data=Dict{Tuple,UInt}(
                                                             (false, false) => 1,
                                                             (false, true) => 1,
                                                             (true, false) => 1,
                                                             (true, true) => 2,
                                                         )),
                                         2 => TruthTable(;
                                                         conditions=[function (state_machine)
                                                                         return params["pv_plant"].supply <
                                                                                params["pv_threshold"]
                                                                     end,
                                                                     function (state_machine)
                                                                         return params["battery"].load >=
                                                                                params["battery"].capacity *
                                                                                params["discharge_limit"]
                                                                     end],
                                                         table_data=Dict{Tuple,UInt}(
                                                             (false, false) => 1,
                                                             (false, true) => 1,
                                                             (true, false) => 1,
                                                             (true, true) => 2,
                                                         )),
                                     ))

        return new("economical_discharge", params, state_machine)
    end

    function CM_EconomicalDischarge()
        return new("economical_discharge", Dict{String,Any}(), StateMachine())
    end
end

function has_method_for(mod::CM_EconomicalDischarge, func::ControlModuleFunction)::Bool
    return func == cmf_charge_is_allowed ||
           func == cmf_discharge_is_allowed
end

function update(mod::CM_EconomicalDischarge)
    move_state(mod.state_machine)
end

function charge_is_allowed(mod::CM_EconomicalDischarge, sim_params::Dict{String,Any})::Bool
    return mod.state_machine.state == 1
end

function discharge_is_allowed(mod::CM_EconomicalDischarge,
                              sim_params::Dict{String,Any})::Bool
    return mod.state_machine.state == 2
end
