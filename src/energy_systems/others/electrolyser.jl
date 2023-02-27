"""
Implementation of an electrolyser, turning electricity and water into H2, O2 and heat.

For the moment this remains a simple implementation that converts electricity into
the gases and heat (as medium m_h_w_ht1) at a defined ratio of 1:0.6:0.4. Has a minimum
run time of 3600s taken into consideration in its control behaviour and a minimum power
fraction of 20%. The power is considered the maximum amount of electricity that the
electrolyser can consume.

At the moment there is no operation strategy is implemented and the production of the
electrolyser is controlled by the demand it is linked to requires.
"""
mutable struct Electrolyser <: ControlledSystem
    uac :: String
    controller :: Controller
    sys_function :: SystemFunction

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap
    medium_names :: Dict{String, String}

    power :: Float64
    heat_fraction :: Float64
    min_power_fraction :: Float64
    min_run_time :: UInt
    output_temperature :: Temperature

    function Electrolyser(uac :: String, config :: Dict{String, Any})
        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
            ),
            sf_transformer, # sys_function
            InterfaceMap( # input_interfaces
                MediumCategoryMap["el_in" in keys(config["medium_names"]) ? config["medium_names"]["el_in"] : "m_e_ac_230v"] => nothing
            ),
            InterfaceMap( # output_interfaces
                MediumCategoryMap["heat_out" in keys(config["medium_names"]) ? config["medium_names"]["heat_out"] : "m_h_w_lt1"] => nothing,
                MediumCategoryMap["h2_out" in keys(config["medium_names"]) ? config["medium_names"]["h2_out"] : "m_c_g_h2"] => nothing,
                MediumCategoryMap["o2_out" in keys(config["medium_names"]) ? config["medium_names"]["o2_out"] : "m_c_g_o2"] => nothing
            ),
            "medium_names" in keys(config) ? # medium_names for input and outputs
                Dict{String, String}(
                    "el_in" => "el_in" in keys(config["medium_names"]) ? config["medium_names"]["el_in"] : "m_e_ac_230v",
                    "heat_out" => "heat_out" in keys(config["medium_names"]) ? config["medium_names"]["heat_out"] : "m_h_w_lt1",
                    "h2_out" => "h2_out" in keys(config["medium_names"]) ? config["medium_names"]["h2_out"] : "m_c_g_h2",
                    "o2_out" => "o2_out" in keys(config["medium_names"]) ? config["medium_names"]["o2_out"] : "m_c_g_o2"
                ) :
                Dict{String, String}(  # default medium_names if no "medium_names" dict is given in input file
                    "el_in" => "m_e_ac_230v",
                    "heat_out" => "m_h_w_lt1",
                    "h2_out" => "m_c_g_h2",
                    "o2_out" =>  "m_c_g_o2"
                ),
            config["power"], # power
            "heat_fraction" in keys(config) # heat_fraction
                ? config["heat_fraction"]
                : 0.4,
            "min_power_fraction" in keys(config) # min_power_fraction
                ? config["min_power_fraction"]
                : 0.2,
            "min_run_time" in keys(config) # min_run_time
                ? config["min_run_time"]
                : 3600,
            "output_temperature" in keys(config) # output_temperature
                ? config["output_temperature"]
                : 55.0
        )
    end
end

function produce(unit :: Electrolyser, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    # abbreviations for media types of input and outputs
    heat_out = MediumCategoryMap[unit.medium_names["heat_out"]]
    h2_out = MediumCategoryMap[unit.medium_names["h2_out"]]
    o2_out = MediumCategoryMap[unit.medium_names["h2_out"]]
    el_in = MediumCategoryMap[unit.medium_names["el_in"]]

    max_produce_h = watt_to_wh(unit.power * unit.heat_fraction)
    max_produce_g = watt_to_wh(unit.power * (1.0 - unit.heat_fraction))
    max_available_e = unit.power

    # heat
    balance_h, potential_h, _ = balance_on(
        unit.output_interfaces[heat_out],
        unit.output_interfaces[heat_out].target
    )

    # hydrogen
    balance_g, potential_g, _ = balance_on(
        unit.output_interfaces[h2_out],
        unit.output_interfaces[h2_out].target
    )   

    # electricity 
    balance_e, potential_e, _ = balance_on(
        unit.input_interfaces[el_in],
        unit.input_interfaces[el_in].target
    )

    if balance_h + potential_h >= 0.0 
        return # don't add to a surplus of h2 
    end

    # --> currently not working as potential of unlimited sinks are not written into interface  @ToDo
    # if  balance_g + potential_g >= 0.0 
    #     return # don't add to a surplus of heat
    # end

    # --> currently not working as potential of unlimited sources are not written into interface @ToDo
    # if balance_e + potential_e <= 0.0
    #     return  # no elecricity available
    # end   

    # --> currently not working as balances are not calculated correctly for unlimited gas and electricity @ToDo
    #usage_fraction = min(1.0, abs(balance_h + potential_h) / max_produce_h, abs(balance_g + potential_g) / max_produce_g, abs(balance_e + potential_e) / max_available_e)
    # for now, use only heat balance and potential:
    usage_fraction = min(1.0, abs(balance_h + potential_h) / max_produce_h)

    if usage_fraction < unit.min_power_fraction
        return
    end

    # @TODO: handle O2 calculation if it ever becomes relevant. for now use molar ratio
    add!(unit.output_interfaces[h2_out], max_produce_g * usage_fraction)
    add!(unit.output_interfaces[o2_out], max_produce_g * usage_fraction * 0.5)
    add!(
        unit.output_interfaces[heat_out],
        max_produce_h * usage_fraction,
        unit.output_temperature
    )
    sub!(unit.input_interfaces[el_in], watt_to_wh(unit.power * usage_fraction))

end

export Electrolyser