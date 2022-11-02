"""
Implementation of a photovoltaic (PV) power plant.

For the moment this remains a simple implementation approximating a PV plant with a sinoid
function. As the calculation of potential PV power is done outside the simulation by a
seperate tool, a proper implemention would mostly just load a profile and consider only
some system losses. The amplitude parameter is a scaling factor, but is not an average
power value.
"""
Base.@kwdef mutable struct PVPlant <: ControlledSystem
    uac :: String
    controller :: StateMachine
    sys_function :: SystemFunction

    input_interfaces :: InterfaceMap
    output_interfaces :: InterfaceMap

    amplitude :: Float64
end

function make_PVPlant(amplitude :: Float64) :: PVPlant
    return PVPlant(
        StateMachine(), # controller
        limited_source, # sys_function
        InterfaceMap(), # input_interfaces
        InterfaceMap( # output_interfaces
            m_e_ac_230v => nothing
        ),
        amplitude, # amplitude
    )
end

function output_values(unit :: PVPlant) :: Vector{String}
    return ["OUT", "Power"]
end

function output_value(unit :: PVPlant, key :: OutputKey) :: Float64
    if key.value_key == "OUT"
        return unit.output_interfaces[key.medium].sum_abs_change * 0.5
    elseif key.value_key == "Power"
        return unit.power # @TODO: Save the last calculated value and return it
    end
    raise(KeyError(key.value_key))
end

function produce(unit :: PVPlant, parameters :: Dict{String, Any}, watt_to_wh :: Function)
    outface = unit.output_interfaces[m_e_ac_230v]
    add!(outface, watt_to_wh(power_at_time(unit, parameters["time"])))
end

function power_at_time(plant :: PVPlant, time :: Int) :: Float64
    seconds_in_day = 60 * 60 * 24
    base_sine = Base.Math.sin(Base.MathConstants.pi * (time % seconds_in_day) / seconds_in_day)
    return Base.Math.max(0.0, Base.Math.min(
        plant.amplitude,
        1.4 * plant.amplitude * base_sine * base_sine * base_sine - 0.2 * plant.amplitude
    ))
end

export PVPlant, make_PVPlant, output_values, output_value