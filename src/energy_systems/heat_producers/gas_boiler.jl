"""
Implementation of a gas boiler producing heat from chemical energy in gaseous form.

The only currently implemented operation strategy involves checking the load of a linked
buffer tank and en-/disabling the boiler when a threshold is reached, in addition to an
overfill shutoff condition.
"""
mutable struct GasBoiler <: ControlledComponent
    uac::String # user address code -> e.g. location and affiliation of the units within the buildings
    controller::Controller # control mechanism for the operating strategy
    sys_function::SystemFunction # function of the component

    input_interfaces::InterfaceMap # interfaces where gas boiler is reciever
    output_interfaces::InterfaceMap # interfaces where gas boiler is source

    m_gas_in::Symbol
    m_heat_out::Symbol

    power::Float64
    load::String
    #max_thermal_efficiency::Float64
    max_consumable_gas::Float64
    min_power_fraction::Float64 # Minimum amount of power needed for the component to run (?)
    min_run_time::UInt # von state machine genutzt; min run time of component, e.g. several time steps
    output_temperature::Temperature # Desired temperatue of the heated medium

    # constructor: initialize object's attributes
    function GasBoiler(uac::String, config::Dict{String,Any})
        m_gas_in = Symbol(default(config, "m_gas_in", "m_c_g_natgas")) # Get value of key
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_ht1")) # Get value of key
        register_media([m_gas_in, m_heat_out])
        max_consumable_gas = watt_to_wh(float(config["power"])) / config["max_thermal_efficiency"]
    
        return new(
            uac, # uac
            controller_for_strategy( # controller
                config["strategy"]["name"], config["strategy"]
            ),
            sf_transformer, # sys_function
            InterfaceMap( # input_interfaces
                m_gas_in => nothing
            ),
            InterfaceMap( # output_interfaces
                m_heat_out => nothing
            ),
            m_gas_in,
            m_heat_out,
            config["power"], # power
            config["load"], # "full" or "plr" (part load ratio)
            max_consumable_gas,
            default(config, "min_power_fraction", 0.1),
            default(config, "min_run_time", 0),
            default(config, "output_temperature", nothing)
        )
    end
end

function control(
    unit::GasBoiler,
    components::Grouping,
    parameters::Dict{String,Any}
)
    move_state(unit, components, parameters)
    unit.output_interfaces[unit.m_heat_out].temperature = highest_temperature(unit.output_temperature, unit.output_interfaces[unit.m_heat_out].temperature)
end

function set_max_energies!(unit::GasBoiler, gas_in::Float64, heat_out::Float64) # Set the maximum energies that can be taken in and outputed from the unit
    set_max_energy!(unit.input_interfaces[unit.m_gas_in], gas_in)
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], heat_out)
end

function check_gas_in(
    unit::GasBoiler,
    parameters::Dict{String,Any}
)
    if unit.controller.parameter["m_gas_in"] == true # soll Steuerparameter das Medium für die Berechnung berücksichtigen
        if ( # If unit's system function is to transform energy of one medium to energy in another medium
            unit.input_interfaces[unit.m_gas_in].source.sys_function == sf_transformer
            && # and max energy, that can be outputed (does it mean nothing gets lost during conversion?), is not bounded
            unit.input_interfaces[unit.m_gas_in].max_energy === nothing
        )
            return (Inf, Inf) # Get as much energy as needed from direct covarage or storage (?)
        else # bounded source or fixed source -> gibt vor wie viel aufgenommen und abgeben werden kann
            exchange = balance_on(
                unit.input_interfaces[unit.m_gas_in],
                unit.input_interfaces[unit.m_gas_in].source
            )
            # potential_energy: how much from direct coverage, e.g. grid
            potential_energy_gas = exchange.balance + exchange.energy_potential
            # potential_storage: how much energy out of storage coverage
            potential_storage_gas = exchange.storage_potential
            # if energy inputed into gas boiler is zero, then return nothing
            if (unit.controller.parameter["unload_storages"] ? potential_energy_gas + potential_storage_gas : potential_energy_gas) <= parameters["epsilon"]
                return (nothing, nothing)
            end
            # return energy from direct coverage, energy from storage and if given the temperature of interface (medium transported ?) Where is temperature used? Why return here?
            return (potential_energy_gas, potential_storage_gas, exchange.temperature)
        end
    else # Nimm so viel du brauchst und berücksichtige nicht in der Rechnung
        return (Inf, Inf)
    end
end

function check_heat_out(
    unit::GasBoiler,
    parameters::Dict{String,Any}
)
    if unit.controller.parameter["m_heat_out"] == true
        exchange = balance_on(
            unit.output_interfaces[unit.m_heat_out],
            unit.output_interfaces[unit.m_heat_out].target
        )
        # potential_energy: how much energy needs to be provided directly
        potential_energy_heat_out = exchange.balance + exchange.energy_potential
        # potential_storage: how much energy can be delivered to storage
        potential_storage_heat_out = exchange.storage_potential
        # if energy outputed from gas boiler is zero, then return nothing
        if (unit.controller.parameter["load_storages"] ? potential_energy_heat_out + potential_storage_heat_out : potential_energy_heat_out) >= -parameters["epsilon"]
            return (nothing, nothing)
        end
        return (potential_energy_heat_out, potential_storage_heat_out)
    else # Provide as much heat as needed for target
        return (-Inf, -Inf)
    end
end

function calculate_thermal_efficiency(
    part_load_ratio::Float64,
    a::Float64 = 0.3822,
    b::Float64 = 2.2013,
    c::Float64 = -2.8237,
    d::Float64 = 1.3021
)
    return  a + b * part_load_ratio + c * part_load_ratio^2 + d * part_load_ratio^3
end

function calculate_inverse_expended_energy(
    intake_gas::Float64,
    max_consumable_gas::Float64
)
    # expended energy function: power = plr1*a in [0,0.1] / plr1[0.1]*a + plr2*b in [0,0.8] / plr1[0.1]*a + plr2[0.8]*b + plr3*c in [0,0.1]
    if intake_gas > max_consumable_gas
        intake_gas = max_consumable_gas
    end
    # normalize
    intake_gas_norm = intake_gas / max_consumable_gas
    
    if intake_gas_norm <= 0.05  # Inversion for the first segment
        return intake_gas_norm / 2.5
    elseif intake_gas_norm > 0.05 && intake_gas_norm <= 0.9 # Inversion for the second segment
        return (intake_gas_norm - 31/880) / (65/88)
    else  # Inversion for the third segment
        return (intake_gas_norm + 2) / 3
    end
end

function calculate_energies(
    unit::GasBoiler,
    parameters::Dict{String,Any},
    potentials::Vector{Float64}
)
    potential_energy_gas_in = potentials[1]
    potential_storage_gas_in = potentials[2]
    potential_energy_heat_out = potentials[3]
    potential_storage_heat_out = potentials[4]

    if unit.load == "full"
        max_produce_heat = watt_to_wh(unit.power)
        max_consume_gas = max_produce_heat
    elseif unit.load == "plr" # part load ratio
        if unit.controller.strategy == "demand_driven"
            max_produce_heat = watt_to_wh(unit.power) # max_produce_heat should equal rated power, else when using demand heat which can be inf a non-physical state occurs
            demand_heat = -(unit.controller.parameter["load_storages"] ? potential_energy_heat_out + potential_storage_heat_out : potential_energy_heat_out) # given
            if demand_heat >= max_produce_heat # limit demand_heat to get a physical value for max_consume_gas, else possible intake can be inf
                demand_heat = max_produce_heat
            elseif demand_heat < 0 # ensure that part-load-ratio is greater equal 0
                demand_heat = 0
            end
            part_load_ratio = demand_heat / watt_to_wh(unit.power) # since demand_heat is bounded, max part-load-ratio is 1, which corresponds to the defintion
            thermal_efficiency = calculate_thermal_efficiency(part_load_ratio)
            max_consume_gas = max_produce_heat / thermal_efficiency      
        elseif unit.controller.strategy == "supply_driven"
            intake_gas = unit.controller.parameter["unload_storages"] ? potential_energy_gas_in + potential_storage_gas_in : potential_energy_gas_in
            max_consume_gas = unit.max_consumable_gas
            part_load_ratio = calculate_inverse_expended_energy(intake_gas, max_consume_gas) # über inverse da Kurve bekannt # given curve of expended energy
            thermal_efficiency = calculate_thermal_efficiency(part_load_ratio)
            max_produce_heat = thermal_efficiency * max_consume_gas
            if max_produce_heat > watt_to_wh(unit.power) # check if computed max_produce_heat exceeds the limit of GB
                max_produce_heat = watt_to_wh(unit.power)
            end
        end
    end

    # get usage fraction of external profile (normalized from 0 to 1)
    # when no profile provided then 100 % of unit is used, else the specified value at time t in the profile is adopted
    usage_fraction_operation_profile = unit.controller.parameter["operation_profile_path"] === nothing ? 1.0 : value_at_time(unit.controller.parameter["operation_profile"], parameters["time"])
    if usage_fraction_operation_profile <= 0.0 # if value of usage fraction is given but <= 0 then component does not run
        return # no operation allowed from external profile
    end

    # all three standard operating strategies behave the same, but it is better to be
    # explicit about the behaviour rather than grouping all together
    if unit.controller.strategy == "storage_driven" && unit.controller.state_machine.state == 2  # what is state 2 (?)
        # fraction should be positive hence multiply the usage_fraction_heat_out by minus one because energies that are outputed are negative, since they are "lost" 
        usage_fraction_heat_out = -((unit.controller.parameter["load_storages"] ? potential_energy_heat_out + potential_storage_heat_out : potential_energy_heat_out) / max_produce_heat)
        usage_fraction_gas_in = +((unit.controller.parameter["unload_storages"] ? potential_energy_gas_in + potential_storage_gas_in : potential_energy_gas_in) / max_consume_gas)

    elseif unit.controller.strategy == "storage_driven"
        return (false, nothing, nothing, nothing)

    elseif unit.controller.strategy == "supply_driven"
        usage_fraction_heat_out = -((unit.controller.parameter["load_storages"] ? potential_energy_heat_out + potential_storage_heat_out : potential_energy_heat_out) / max_produce_heat)
        usage_fraction_gas_in = +((unit.controller.parameter["unload_storages"] ? potential_energy_gas_in + potential_storage_gas_in : potential_energy_gas_in) / max_consume_gas)

    elseif unit.controller.strategy == "demand_driven"
        usage_fraction_heat_out = -((unit.controller.parameter["load_storages"] ? potential_energy_heat_out + potential_storage_heat_out : potential_energy_heat_out) / max_produce_heat)
        usage_fraction_gas_in = +((unit.controller.parameter["unload_storages"] ? potential_energy_gas_in + potential_storage_gas_in : potential_energy_gas_in) / max_consume_gas)
    end
    
    # limit actual usage by limits of inputs, outputs and profile
    # goal: exclude unphysical usage/calculations
    usage_fraction = min(
        1.0,
        usage_fraction_heat_out, # 1.0 if all heat energy can be outputed from the gas boiler || Can actually be higher because demand might be higher (???)
        usage_fraction_gas_in, # 1.0 if all gas can be taken in / used for heat generation
        usage_fraction_operation_profile # 1.0 if deposited like in this in the profile
    )
    
    if usage_fraction < unit.min_power_fraction # Komponente wird nicht genutzt, wenn der Nutzungsanteil geringer ist als minimal zu aufwendende Leistung, um das Gerät zubenutzen
        return (false, nothing, nothing)
    end
               
    return (
        true,
        max_consume_gas * usage_fraction, 
        max_produce_heat * usage_fraction
    )
end

function potential(
    unit::GasBoiler,
    parameters::Dict{String,Any}
)
    # Komponente schreibt in "Register" wie viel Energie es abgeben und aufnehmen kann, damit andere Komponente später mit den Werten "arbeiten" können

    # bestimme Energie, die Gas boiler zur Verfügung steht
    potential_energy_gas_in, potential_storage_gas_in = check_gas_in(unit, parameters)
    # Wenn keine Energie zugeführt wird, dann kann nichts aufgenommen werden und nichts weitergegeben werden
    if potential_energy_gas_in === nothing && potential_storage_gas_in === nothing
        set_max_energies!(unit, 0.0, 0.0)
        return
    end

    potential_energy_heat_out, potential_storage_heat_out = check_heat_out(unit, parameters)
    # Wenn keine Energie abgegeben werden kann, dann wird Aufnahme und Ausgabe von Energie auf Null gesetzt
    if potential_energy_heat_out === nothing && potential_storage_heat_out === nothing
        set_max_energies!(unit, 0.0, 0.0)
        return
    end

    # Wenn Energie aufgenommen wird und abgegeben werden kann, dann berechne den aufgenommenen Anteil und den zur Abgabe möglichen, respektive
    energies = calculate_energies(
        unit, parameters,
        [
            potential_energy_gas_in, potential_storage_gas_in,
            potential_energy_heat_out, potential_storage_heat_out
        ]
    )

    if !energies[1] # if there is no energy returned, then set energy, that is needed and can be returned, to zero.
        set_max_energies!(unit, 0.0, 0.0)
    else # else assign the calculated values for input into gas boiler and output of gas boiler
        set_max_energies!(unit, energies[2], energies[3])
    end
end

function process(unit::GasBoiler, parameters::Dict{String,Any})
    # Berechne die tatsächliche Energie durch Nutzung der Informationen im "Energie-Register" und update die Energiewerte
    # Energiewerte sind bereits hinterlegt, die genutzt werden können, gesehen von der Sicht aus des GB
    potential_energy_gas_in, potential_storage_gas_in = check_gas_in(unit, parameters)
    if potential_energy_gas_in === nothing && potential_storage_gas_in === nothing
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
            potential_energy_gas_in, potential_storage_gas_in,
            potential_energy_heat_out, potential_storage_heat_out
        ]
    )
    
    if energies[1]
        sub!(unit.input_interfaces[unit.m_gas_in], energies[2])
        add!(unit.output_interfaces[unit.m_heat_out], energies[3])
    end
end

export GasBoiler