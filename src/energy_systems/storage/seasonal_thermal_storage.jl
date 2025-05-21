"""
Implementation of a seasonal thermal storage component.

This is a simplified model, which mostly deals with amounts of energy and considers
temperatures only for the available temperature as the tank is depleted.
"""

using Plots: Plots
using Dates
using PlotlyJS: PlotlyJS

mutable struct SeasonalThermalStorage <: Component
    # general
    uac::String
    controller::Controller
    sys_function::SystemFunction

    input_interfaces::InterfaceMap
    output_interfaces::InterfaceMap
    m_heat_in::Symbol
    m_heat_out::Symbol

    # geometry and physical properties
    capacity::Float64
    volume::Float64
    hr_ratio::Float64
    sidewall_angle::Float64
    shape::String
    rho_medium::Float64
    cp_medium::Float64
    diffussion_coefficient::Float64
    surface_area_lid::Float64
    surface_area_bottom::Float64
    surface_area_barrel_segments::Vector{Float64}
    volume_segments::Vector{Float64}
    height::Float64
    radius_small::Float64
    radius_large::Float64
    number_of_layer_total::Int64
    number_of_layer_above_ground::Int64
    output_layer_from_top::Int64
    dz::Vector{Float64}
    dz_normalized::Vector{Float64}
    sigma::Vector{Float64}
    lambda::Vector{Float64}
    phi::Vector{Float64}
    theta::Vector{Float64}
    layer_masses::Vector{Float64}

    # loading and unloading
    high_temperature::Temperature
    low_temperature::Temperature
    max_load_rate::Floathing
    max_unload_rate::Floathing
    max_input_energy::Floathing
    max_output_energy::Floathing

    # losses
    thermal_transmission_lid::Float64
    thermal_transmission_barrel::Float64
    thermal_transmission_bottom::Float64
    ambient_temperature_profile::Union{Profile,Nothing}
    ambient_temperature::Temperature
    ground_temperature::Temperature
    effective_ambient_temperature::Vector{Temperature}

    # state variables
    current_max_output_temperature::Float64
    current_min_input_temperature::Float64
    initial_load::Float64
    load::Float64
    load_end_of_last_timestep::Float64
    losses::Float64
    temperature_segments::Vector{Temperature}
    current_energy_input::Vector{Float64}
    current_temperature_input::Vector{Temperature}
    current_energy_output::Float64
    current_temperature_output::Temperature
    current_max_output_energy::Float64
    temperatures_charging::Vector{Temperature}

    # additional output
    temp_distribution_output::Array{Float64}

    function SeasonalThermalStorage(uac::String, config::Dict{String,Any}, sim_params::Dict{String,Any})
        m_heat_in = Symbol(default(config, "m_heat_in", "m_h_w_ht1"))
        m_heat_out = Symbol(default(config, "m_heat_out", "m_h_w_lt1"))
        register_media([m_heat_in, m_heat_out])

        constant_ambient_temperature,
        ambient_temperature_profile = get_parameter_profile_from_config(config,
                                                                        sim_params,
                                                                        "ambient_temperature",
                                                                        "ambient_temperature_profile_file_path",
                                                                        "ambient_temperature_from_global_file",
                                                                        "constant_ambient_temperature",
                                                                        uac;
                                                                        required=true)

        # Note: layer numbering starts at the bottom with index 1 and ends at the top with index number_of_layer_total
        return new(uac,                                            # uac
                   Controller(default(config, "control_parameters", nothing)),
                   sf_storage,                                     # sys_function
                   InterfaceMap(m_heat_in => nothing),             # input_interfaces
                   InterfaceMap(m_heat_out => nothing),            # output_interfaces
                   m_heat_in,                                      # medium in the STES   
                   m_heat_out,                                     # medium out of the STES
                   # geometry and physical properties
                   0.0,                                            # capacity of the STES [Wh]
                   default(config, "volume", nothing),             # volume of the STES [m^3]
                   default(config, "hr_ratio", 0.5),               # ratio of the height to the mean radius of the STES
                   default(config, "sidewall_angle", 0.0),         # angle of the sidewall of the STES with respect to the horizont [°]
                   default(config, "shape", "cuboid"),             # can be "round" for cylinder/truncated cone or "cuboid" for tank or truncated quadratic pyramid (pit)
                   default(config, "rho_medium", 1000.0),          # density of the medium [kg/m^3]
                   default(config, "cp_medium", 4.18),             # specific thermal capacity of medium [kJ/kgK]
                   default(config, "diffussion_coefficient", 0.143 * 10^-6), # diffussion coefficient of the medium [m^2/s]
                   0.0,                                            # surface_area_lid, surface of the lid of the STES [m^2]
                   0.0,                                            # surface_area_bottom, surface of the bottom of the STES [m^2]
                   [],                                             # surface_area_barrel_segments, surface of the barrel segments of the STES [m^2]
                   [],                                             # volume_segments, volume of the segments of the STES [m^3]
                   0.0,                                            # height of the STES [m]
                   0.0,                                            # radius_small, small (lower) radius of the STES [m]
                   0.0,                                            # radius_large, large (upper) radius of the STES [m]
                   default(config, "number_of_layer_total", 25),   # number of layers in the STES
                   default(config, "number_of_layer_above_ground", 5), # number of layers above ground in the STES
                   default(config, "output_layer_from_top", 1),    # layer number of the output layer, counted from the top
                   [],                                             # dz, thickness of the layers of the STES [m]
                   [],                                             # dz_normalized, normalized dz with respect to to the volume of each section
                   [],                                             # [1/h]  factor for losses to ambiente: area_of_losses * U[kJ/m^2K] / (roh * cp * volume_segment)
                   [],                                             # [K/kJ] factor for input/output energy:  1 / (roh * cp * volume_segment ) 
                   [],                                             # [1/kg] factor for input/output mass flow:  1 / (roh  * volume_segment )
                   [],                                             # volume-ratios of sections: V_section[n-1] / (V_section[n] + V_section[n-1])
                   [],                                             # layer_masses, mass of the medium in each layer [kg]
                   # loading and unloading
                   default(config, "high_temperature", 75.0),      # upper temperature of the STES [°C]
                   default(config, "low_temperature", 20),         # lower temperature of the STES [°C]
                   default(config, "max_load_rate", nothing),      # maximum load rate given in 1/h
                   default(config, "max_unload_rate", nothing),    # maximum unload rate given in 1/h
                   nothing,                                        # max_input_energy, maximum input energy per time step [Wh]
                   nothing,                                        # max_output_energy, maximum output energy per time step [Wh]
                   # Losses
                   default(config, "thermal_transmission_lid", 0.25),     # [W/(m^2K)]
                   default(config, "thermal_transmission_barrel", 0.375), # [W/(m^2K)]
                   default(config, "thermal_transmission_bottom", 0.375), # [W/(m^2K)]
                   ambient_temperature_profile,                           # [°C]
                   constant_ambient_temperature,                          # ambient_temperature [°C]
                   default(config, "ground_temperature", 12),             # [°C]
                   [],                                                    # effective_ambient_temperature corresponding to each layer [°C]          
                   # other
                   0.0,                                            # current_max_output_temperature
                   0.0,                                            # current_min_input_temperature
                   default(config, "initial_load", 0.0),           # initial_load [%/100] assuming perfectly mixed storage at the begin
                   0.0,                                            # load, set to initial_load at the beginning [Wh]
                   0.0,                                            # load_end_of_last_timestep, stores the load of the previous time step without losses
                   0.0,                                            # losses in current time step [Wh]
                   [],                                             # temperature_segments: temperatures of the segments
                   [],                                             # current_energy_input, energy input in current time step [Wh]
                   [],                                             # current_temperature_input, temperature of the input in current time step [°C]
                   0.0,                                            # current_energy_output, energy output in current time step [Wh]
                   0.0,                                            # current_temperature_output, temperature of the output in current time step [°C]
                   0.0,                                            # current_max_output_energy, maximum output energy in current time step [Wh]
                   [],                                             # temperatures_charging, temperatures of the possible inputs from exchange in current time step [°C]
                   # additional output
                   Array{Float64}(undef, 0, 0))                    # temp_distribution_output [°C], holds temperature field of layers for output plot
    end
end

function initialise!(unit::SeasonalThermalStorage, sim_params::Dict{String,Any})
    set_storage_transfer!(unit.input_interfaces[unit.m_heat_in],
                          unload_storages(unit.controller, unit.m_heat_in))
    set_storage_transfer!(unit.output_interfaces[unit.m_heat_out],
                          load_storages(unit.controller, unit.m_heat_out))

    # set temperature vector: assuming a uniform temperature profile (mixed storage)
    mean_temperature = unit.initial_load * (unit.high_temperature - unit.low_temperature) + unit.low_temperature
    unit.temperature_segments = [mean_temperature for _ in 1:(unit.number_of_layer_total)]

    # set initial temperature bounds
    set_temperature_limits!(unit, sim_params)

    # calculate thermal transmission coefficients # TODO other transmission below ground
    thermal_transmission_barrels = [unit.thermal_transmission_barrel for _ in 1:(unit.number_of_layer_total)]

    # calculate geometry of the STES
    unit.surface_area_lid,
    unit.surface_area_barrel_segments,
    unit.surface_area_bottom,
    unit.volume_segments,
    unit.height,
    unit.radius_small,
    unit.radius_large = calc_STES_geometry(unit.uac,
                                           unit.volume,
                                           unit.sidewall_angle,
                                           unit.hr_ratio,
                                           unit.number_of_layer_total,
                                           unit.shape)

    # calculate the mass of the medium in each layer
    unit.layer_masses = unit.rho_medium .* unit.volume_segments

    # calculate (equally spaced) thickness of the layers of the STES [m]
    unit.dz = fill(unit.height / unit.number_of_layer_total, unit.number_of_layer_total)
    # calculate the normalized dz with respect to to the volume of each section
    unit.dz_normalized = unit.dz .* sqrt.((unit.volume_segments .* unit.number_of_layer_total) ./ unit.volume)

    # calculate coefficient for losses to ambiente
    unit.sigma = zeros(unit.number_of_layer_total)
    unit.sigma[1] = unit.surface_area_bottom * unit.thermal_transmission_bottom /
                    (unit.rho_medium * convert_kJ_in_Wh(unit.cp_medium) * unit.volume_segments[1])          # [1/h] losses to ambiente through bottom
    unit.sigma[end] = unit.surface_area_lid * unit.thermal_transmission_lid /
                      (unit.rho_medium * convert_kJ_in_Wh(unit.cp_medium) * unit.volume_segments[end])      # [1/h] losses to ambiente through lid
    unit.sigma = unit.sigma .+
                 unit.surface_area_barrel_segments .* thermal_transmission_barrels ./
                 (unit.rho_medium * convert_kJ_in_Wh(unit.cp_medium) * unit.volume_segments)                # [1/h]  losses to ambiente though barrel

    # calculate coefficient for input/output energy --> needed? TODO
    unit.lambda = 1 ./ (unit.rho_medium * convert_kJ_in_Wh(unit.cp_medium) * unit.volume_segments)          # [K/Wh] 

    # coefficient for input/output mass flow, assuming water as fluid
    cp_water = 4.18                                                                      # [kJ/kgK]
    unit.phi = cp_water ./ (unit.cp_medium * unit.rho_medium * unit.volume_segments)      # [1/kg]

    # coefficient for buoyancy effects
    unit.theta = [unit.volume_segments[n - 1] / (unit.volume_segments[n] + unit.volume_segments[n - 1])
                  for n in 2:(unit.number_of_layer_total)]
    pushfirst!(unit.theta, 0.0)  # Set first element to 0

    unit.capacity = unit.volume * unit.rho_medium * convert_kJ_in_Wh(unit.cp_medium) *
                    (unit.high_temperature - unit.low_temperature)  # [Wh]
    unit.load = unit.initial_load * unit.capacity
    unit.load_end_of_last_timestep = copy(unit.load)

    # calculate maximum input and output energy
    if unit.max_load_rate === nothing
        unit.max_input_energy = Inf
    else
        unit.max_input_energy = unit.max_load_rate * unit.capacity / (sim_params["time_step_seconds"] / 60 / 60)     # [Wh per timestep]
    end
    if unit.max_unload_rate === nothing
        unit.max_output_energy = Inf
    else
        unit.max_output_energy = unit.max_unload_rate * unit.capacity / (sim_params["time_step_seconds"] / 60 / 60)  # [Wh per timestep]
    end

    # vector to hold the results of the temperatures for each layer in each simulation time step
    unit.temp_distribution_output = zeros(Float64, sim_params["number_of_time_steps"], unit.number_of_layer_total)
end

"""
    weighted_mean(values::Vector{Float64}, weights::Vector{Float64}, sim_params::Dict{String,Any}) -> Float64

Calculate the weighted mean of a set of values.

# Arguments
- `values::Vector{Any}`: A vector containing the values.
- `weights::Vector{Any}`: A vector containing the weights corresponding to each value.

# Returns
- `Float64`: The weighted mean of the values.
"""
function weighted_mean(values::Vector{<:Any}, weights::Vector{<:Any}, sim_params::Dict{String,Any})
    if abs(sum(weights)) <= sim_params["epsilon"]
        return 0.0
    else
        return sum(values .* weights) / sum(weights)
    end
end

"""
    calc_STES_geometry(uac::String, volume::Float64, alpha::Float64, hr::Float64, n_segments::Int64)

Calculate the geometry of a seasonal thermal energy storage (STES) system, either as a 
- cylinder or a truncated cone (shape == round) 
- tank with quadratic base or a truncated quadratic pyramid (pit) (shape == cuboid)

# Arguments
- `uac::String`: Unique identifier for the STES
- `volume::Float64`: Total volume of the storage [m^3]
- `alpha::Float64`: Slope angle of the truncated cone in degrees with respect to the horizontal. 
                    If `alpha` is 90, the storage is a cylinder.
- `hr::Float64`: Height-to-radius ratio
- `n_segments::Int64`: Number of segments to divide the storage into for calculation
- `shape::String`: Shape of the STES, can be "round" for cylinder/truncated cone or "cuboid" for tank or 
                   truncated quadratic pyramid (pit)

# Returns
- `a_lid::Float64`: Surface area of the top lid [m^2]
- `a_barrel::Vector{Float64}`: Surface areas of the barrel sections [m^2]
- `a_bottom::Float64`: Surface area of the bottom lid [m^2]
- `v_section::Vector{Float64}`: Volumes of the sections [m^3]
- `height::Float64`: Height of the storage [m]
- `radius_small::Float64`: Radius of the smaller base (for truncated cone) [m]
- `radius_large::Float64`: Radius of the larger base (for truncated cone) [m]

"""
function calc_STES_geometry(uac::String, volume::Float64, alpha::Float64, hr::Float64, n_segments::Int64, shape::String)
    a_barrel = zeros(n_segments)
    v_section = zeros(n_segments)

    # helper functions
    deg2rad(deg) = deg * (π / 180)
    rad2deg(rad) = rad * (180 / π)

    if alpha == 90 && shape == "round"  # cylinder
        # calculate radius and height of cylinder
        radius = cbrt(volume / (pi * hr))
        height = hr * radius

        # calculate surfaces and volumes of the sections
        a_lid = pi * radius^2
        for n in 1:n_segments
            a_barrel[n] = 2 * pi * radius * height / n_segments
            v_section[n] = pi * radius^2 * height / n_segments
        end
        a_bottom = a_lid

        return a_lid, a_barrel, a_bottom, v_section, height, radius, radius
    elseif alpha == 90 && shape == "cuboid"  # cuboid with square cross-section 
        # calculate radius and height of cuboid
        side_length = cbrt(2 * volume / hr)
        height = hr * (side_length / 2)

        # calculate surfaces and volumes of the sections
        a_lid = side_length^2
        for n in 1:n_segments
            a_barrel[n] = 4 * side_length * height / n_segments
            v_section[n] = side_length^2 * height / n_segments
        end
        a_bottom = a_lid

        return a_lid, a_barrel, a_bottom, v_section, height, side_length / 2, side_length / 2
    elseif shape == "round"   # truncated cone
        alpha_rad = deg2rad(alpha)
        alpha_tan = tan(alpha_rad)

        radius_large = cbrt((3 * volume) / (pi * alpha_tan * (1 - ((2 * alpha_tan - hr) / (2 * alpha_tan + hr))^3)))
        radius_small = radius_large * (2 * alpha_tan - hr) / (2 * alpha_tan + hr)

        if radius_small < 0
            alpha_min = rad2deg(atan(hr / 2))
            alpha_max = 180 - alpha_min
            hr_max = 2 * alpha_tan

            @error "For the STES $(uac), reduce the h/r ratio for the truncated cone or increase the slope angle. " *
                   "For a given h/r ratio of $(hr), the slope angle must be greater than $(round(alpha_min; digits=2))° and " *
                   "less than $(round(alpha_max; digits=2))°. For a given slope angle of $(alpha)°, the h/r ratio must be " *
                   "less than $(round(hr_max; digits=2))."
            throw(InputError)
        end

        height = hr * (radius_large + radius_small) / 2
        a_lid = pi * radius_large^2
        a_bottom = pi * radius_small^2

        for n in 1:n_segments
            radius_seg_large = radius_large - (n_segments - n) * (radius_large - radius_small) / n_segments
            radius_seg_small = radius_small + (n - 1) * (radius_large - radius_small) / n_segments
            a_barrel[n] = (radius_seg_large + radius_seg_small) * pi *
                          sqrt((radius_seg_large - radius_seg_small)^2 + (height / n_segments)^2)
            v_section[n] = (height / n_segments) * pi / 3 *
                           (radius_seg_large^2 + radius_seg_large * radius_seg_small + radius_seg_small^2)
        end

        return a_lid, a_barrel, a_bottom, v_section, height, radius_small, radius_large
    elseif shape == "cuboid"   # truncated quadratic pyramid (pit)
        alpha_rad = deg2rad(alpha)
        alpha_tan = tan(alpha_rad)

        # check if input parameters are creating a valid truncated quadratic pyramid
        hr_max = 2 * alpha_tan
        if hr > hr_max
            alpha_min_deg = rad2deg(atan(hr / 2))
            alpha_max_deg = 180 - alpha_min_deg

            @error "For the STES $(uac), reduce the h/r ratio for the truncated quadratic pyramit or increase the slope angle. " *
                   "For a given h/r ratio of $(hr), the slope angle must be greater than $(round(alpha_min_deg; digits=2))° and " *
                   "smaller than $(round(alpha_max_deg; digits=2))°. " *
                   "For a given slope angle of $(alpha)°, the h/r ratio must be less than $(round(hr_max; digits=2))."
            throw(InputError())
        end

        # top side has bigger cross section than bottom side to form a pit: \_/
        # hr = h/r_mean  ⇒  λ = t/b, Volume = (hr·b³/12)·(1+λ)(1+λ+λ^2)
        lam = (2 * alpha_tan - hr) / (2 * alpha_tan + hr)
        top_side = (12 * volume / (hr * (1 + lam) * (1 + lam + lam^2)))^(1 / 3)
        base_side = lam * top_side
        height = hr * (base_side + top_side) / 4

        # areas
        a_bottom = base_side^2
        a_top = top_side^2

        # segments
        dh = height / n_segments
        a_lat = zeros(Float64, n_segments)
        v_section = zeros(Float64, n_segments)

        for i in 1:n_segments
            z_bot = dh * (i - 1)
            z_top = dh * i

            # lineare Interpolation of edge length
            a_bot = base_side + (top_side - base_side) * (z_bot / height)
            a_up = base_side + (top_side - base_side) * (z_top / height)
            # volume of segments
            v_section[i] = (dh / 3) * (a_bot^2 + a_bot * a_up + a_up^2)

            # Segment-laterale area
            l_seg = sqrt(((a_bot - a_up) / 2)^2 + dh^2)
            a_lat[i] = 2 * (a_bot + a_up) * l_seg
        end
        return a_top, a_lat, a_bottom, v_section, height, base_side / 2, top_side / 2
    else
        @error "Invalid shape type of seasonal thermal storage $(unit.uac). Shape has to be 'round' or 'cuboid'!"
        throw(InputError)
    end
end

function plot_optional_figures_begin(unit::SeasonalThermalStorage, output_path::String, output_formats::Vector{String},
                                     sim_params::Dict{String,Any})
    # Plot geometry of STES
    plt = Plots.plot()
    Plots.plot!([-unit.radius_large, unit.radius_large], [unit.height, unit.height]; color=:blue, lw=6, label="")  # Top
    Plots.plot!([-unit.radius_small, unit.radius_small], [0, 0]; color=:blue, lw=6, label="")  # Bottom
    Plots.plot!([unit.radius_small, unit.radius_large], [0, unit.height]; color=:blue, lw=6, label="")  # Side wall
    Plots.plot!([-unit.radius_small, -unit.radius_large], [0, unit.height]; color=:blue, lw=6, label="")  # Side wall
    Plots.plot!(; title="Cross section of the STES $(unit.uac) (cross-section: $(unit.shape))",
                xlabel="x-coordinate [m]",
                ylabel="height [m]",
                legend=false,
                aspect_ratio=:equal,
                size=(1800, 1200),
                titlefontsize=30,
                guidefontsize=24,
                tickfontsize=24,
                grid=true,
                minorgrid=true,
                gridlinewidth=1,
                margin=15Plots.mm)
    fig_name = "STES_cross_section_$(unit.uac)"
    for output_format in output_formats
        Plots.savefig(output_path * "/" * fig_name * "." * output_format)
    end

    return true
end

function plot_optional_figures_end(unit::SeasonalThermalStorage, sim_params::Dict{String,Any}, output_path::String)
    # Plot temperature distribution over time    
    x_vals_datetime = [add_ignoring_leap_days(sim_params["start_date"],
                                              Dates.Second((s - 1) * sim_params["time_step_seconds"]))
                       for s in 1:sim_params["number_of_time_steps"]]
    layers_to_plot = (unit.number_of_layer_total):-1:1
    traces = PlotlyJS.GenericTrace[]
    for i in layers_to_plot
        plot_label = "Layer $(i)"
        if i == 1
            plot_label = "Bottom layer"
        elseif i == unit.number_of_layer_total
            plot_label = "Top layer"
        end
        trace = PlotlyJS.scatter(; x=x_vals_datetime, y=unit.temp_distribution_output[:, i], mode="lines",
                                 name=plot_label)
        push!(traces, trace)
    end

    leap_days_str = string.([Date(year, 2, 29)
                             for year in
                                 Dates.value(Year(sim_params["start_date"])):Dates.value(Year(sim_params["end_date"]))
                             if isleapyear(year)])

    layout = PlotlyJS.Layout(; title_text="Temperature distribution over time in STES $(unit.uac)",
                             xaxis_title_text="Date",
                             yaxis_title_text="Temperature [°C]",
                             xaxis=PlotlyJS.attr(; type="date",
                                                 rangebreaks=[Dict("values" => leap_days_str)]))
    p = PlotlyJS.plot(traces, layout)

    fig_name = "temperature_distribution_STES_$(unit.uac).html"
    PlotlyJS.savefig(p, output_path * "/" * fig_name)

    return true
end

function control(unit::SeasonalThermalStorage,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)

    # write old temperature field for output
    unit.temp_distribution_output[Int(sim_params["time"] / sim_params["time_step_seconds"]) + 1, :] = copy(unit.temperature_segments)

    # calculate maximum energies for input
    inface = unit.input_interfaces[unit.m_heat_in]
    exchanges = balance_on(inface, inface.source)
    temperatures_charging = highest.([exchange.temperature_max for exchange in exchanges],
                                     [exchange.temperature_min for exchange in exchanges])
    temperatures_charging = [temp === nothing ? unit.high_temperature : temp
                             for temp in temperatures_charging]

    # filter out the elements where the temperature is below the current_min_input_temperature
    mask = temperatures_charging .>= unit.current_min_input_temperature
    temperatures_charging = temperatures_charging[mask]
    source_uac = [exchange.purpose_uac for exchange in exchanges][mask]
    energy_supply = [exchange.balance + exchange.energy_potential for exchange in exchanges][mask]

    # check for directly connected transformers in input interface and set its max energy to Inf
    if inface.source.sys_function == sf_transformer && length(energy_supply) == 1
        energy_supply[1] = Inf
    end

    max_input_energy = zeros(length(temperatures_charging))
    for (input_idx, temperature) in enumerate(temperatures_charging)
        max_input_energy[input_idx] = min(calcuate_max_input_energy_by_temperature(unit,
                                                                                   temperature,
                                                                                   unit.temperature_segments),
                                          min(unit.max_input_energy, energy_supply[input_idx]))
    end
    if isempty(temperatures_charging)
        max_input_energy = [0.0]
        temperatures_charging = [nothing]
        source_uac = [nothing]
    end
    set_max_energy!(unit.input_interfaces[unit.m_heat_in], max_input_energy, temperatures_charging,
                    temperatures_charging, source_uac, false, true)

    # calculate maximum energies for output
    unit.current_max_output_energy = max(min(unit.load, unit.max_output_energy), 0.0)
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], unit.current_max_output_energy, nothing,
                    unit.current_max_output_temperature)

    if unit.ambient_temperature_profile !== nothing
        unit.ambient_temperature = Profiles.value_at_time(unit.ambient_temperature_profile, sim_params)
    end

    # effective ambient temperature for each layer of the STES
    unit.effective_ambient_temperature = vcat(fill(unit.ambient_temperature, unit.number_of_layer_above_ground),
                                              fill(unit.ground_temperature,
                                                   unit.number_of_layer_total - unit.number_of_layer_above_ground))
end

# This funtcion sets the current_min_input_temperature and current_max_output_temperature 
function set_temperature_limits!(unit::SeasonalThermalStorage, sim_params::Dict{String,Any})
    unit.current_max_output_temperature = unit.temperature_segments[end]
    if unit.temperature_segments[1] == unit.temperature_segments[2] && sim_params["time"] == 0
        # add 1K in the first time step if there is equal temperature distribution to avoid Inf mass flow
        unit.current_min_input_temperature = unit.temperature_segments[1] + 1.0
    else
        # limit the input temperature to the temperature of the second layer from below
        if unit.temperature_segments[1] == lowest(unit.temperature_segments)
            unit.current_min_input_temperature = unit.temperature_segments[2]
        else
            # If there is a inverse temperature distrubution within the storage, search for the first layer
            # where the temperature is higher than the temperature of the bottom layer and use this temperature.
            # If the lowest layer is the hottest one, set the current_min_input_temperature to the temperature
            # of the first layer + 1K
            idx = findfirst(x -> x > unit.temperature_segments[1], unit.temperature_segments)
            if isnothing(idx)
                unit.current_min_input_temperature = unit.temperature_segments[1] + 1.0
            else
                unit.current_min_input_temperature = unit.temperature_segments[idx]
            end
        end
    end
end

# This function will only be called from a Bus, never in a 1-to-1 connection to another component.
# energy_input and temperatures_input are all inputs that have alread been put into the STES.
function recalculate_max_energy(unit::SeasonalThermalStorage,
                                energy_input::Union{Float64,Vector{Float64}},
                                temperatures_input::Union{Temperature,Vector{Temperature}},
                                max_energy::EnergySystems.MaxEnergy,
                                sim_params::Dict{String,Any})
    energy_input = isa(energy_input, AbstractVector) ? energy_input : [energy_input]
    temperatures_input = isa(temperatures_input, AbstractVector) ? temperatures_input : [temperatures_input]

    # calculate the new temporal temperature distribution if "energy_input" with "temperatures_input" is put into the STES
    temperature_segments_temporary = update_STES(unit,
                                                 energy_input,
                                                 temperatures_input,
                                                 0.0,
                                                 nothing,
                                                 sim_params;
                                                 temporary_calculation=true)

    # calculate new max_energy for the input energy for all elements in max_energy
    # and limit it by the the alread filled energy and the charging rate limit
    for (input_idx, temperature) in enumerate(highest.(max_energy.temperature_max, max_energy.temperature_min))
        max_energy.max_energy[input_idx] = min(calcuate_max_input_energy_by_temperature(unit,
                                                                                        temperature,
                                                                                        temperature_segments_temporary),
                                               min(unit.max_input_energy - sum(energy_input; init=0.0),
                                                   max_energy.max_energy[input_idx]))
    end

    return max_energy
end

"""
    calcuate_max_input_energy_by_temperature(unit::SeasonalThermalStorage, actual_input_temp::Temperature,
                                            current_temperature_distribution::Vector{Temperature}) -> Float64

Calculates the maximum possible input energy at a given single actual_input_temp for the
given current_temperature_distribution

# Arguments
- `unit::SeasonalThermalStorage`: The seasonal thermal storage unit for which the calculation is performed.
- `actual_input_temp::Temperature`: The temperature of the input energy being added to the storage.
- `current_temperature_distribution::Vector{Temperature}`: A vector representing the current temperature distribution 
  across the storage's volume segments.

# Returns
- `Float64`: The total maximum energy that can be added to the storage at the actual_input_temp
"""
function calcuate_max_input_energy_by_temperature(unit::SeasonalThermalStorage, actual_input_temp::Temperature,
                                                  current_temperature_distribution::Vector{Temperature})
    # reduce maximal energy if the input temperature is near the temperature of the lower layer to avoid
    # numerical problems and pulsing during loading. Otherwise, it can happen that the STES is not taking
    # the energy that it is supposed to take, or a unresonable high temporal resolution has to be used to
    # calculate the new temperature distribution within the STES.
    red_factor = 0.0 + clamp(actual_input_temp - current_temperature_distribution[2], 0.0, 5.0) * 1.0 / 5.0
    return red_factor * sum(convert_mass_in_energy(volume * unit.rho_medium, current_temperature_distribution[layer],
                                                   actual_input_temp, unit.cp_medium)
                            for (layer, volume) in enumerate(unit.volume_segments)
                                if current_temperature_distribution[layer] < actual_input_temp; init=0.0)
end

""" 
    convert_kJ_in_Wh(energy::Float64)

takes energy in [kJ] and convert it to [Wh]
"""
function convert_kJ_in_Wh(energy::Float64)
    return energy / 3.6
end

"""
    convert_energy_in_mass(energy, temp_low, temp_high, cp, roh)

 calculates mass [kg] from energy [Wh] 

 # Arguments
- `energy::Float64`: energy to convert [Wh]
- `temp_low::Temperature`: lower temperatue [°C]
- `temp_high::Temperature`: upper temperatur [°C]
- `cp::Float64`: sppecific heat capacity [kJ/KgK]

"""
function convert_energy_in_mass(energy::Float64, temp_low::Temperature, temp_high::Temperature, cp::Float64)
    if energy == 0.0
        return 0.0
    else
        return energy / (convert_kJ_in_Wh(cp) * (temp_high - temp_low))
    end
end

"""
    convert_mass_in_energy(mass, temp_low, temp_high, cp, roh)

 calculates energy [Wh] from mass [kg] 

 # Arguments
- `mass::Float64`: mass to convert [kg]
- `temp_low::Temperature`: lower temperatue [°C]
- `temp_high::Temperature`: upper temperatur [°C]
- `cp::Float64`: sppecific heat capacity [kJ/KgK]

"""
function convert_mass_in_energy(mass::Float64, temp_low::Temperature, temp_high::Temperature, cp::Float64)
    return mass * convert_kJ_in_Wh(cp) * (temp_high - temp_low)
end

"""
    update_STES(unit::SeasonalThermalStorage,
                energy_input::Vector{Float64},
                temperatures_input::Vector{Temperature},
                energy_output::Float64,
                temperature_output::Temperature,
                sim_params::Dict{String,Any};
                temporary_calculation::Bool=false) -> Tuple{Float64, Float64, Float64, Vector{Float64}}

This function updates the temperature segments of the STES unit by calculating the new temperatures 
for each layer based on thermal diffusion, losses, thermal input/output, and mass input/output
by solving a 1D partial differential equation using an explicit Euler method.
The function also accounts for mixing due to buoyancy effects if a temperature gradient is present.

The method is based on:
Lago, J. et al. (2019): A 1-dimensional continuous and smooth model for thermally stratified storage tanks including 
                        mixing and buoyancy, Applied Energy 248, S. 640-655: doi: 10.1016/j.apenergy.2019.04.139                   
Steinacker, H. (2022): Entwicklung eines dynamischen Simulationsmodells zur Optimierung von wärmegekoppelten 
                       Wasserstoffkonzepten für die klimaneutrale Quartiersversorgung, unpublished master thesis, 
                       University of Stuttgart.

# Arguments
- `unit::SeasonalThermalStorage`: The STES unit to be updated.
- `energy_input::Vector{Float64}`: The energy input to the STES, distributed across different temperatures.
- `temperatures_input::Vector{Temperature}`: The temperatures corresponding to the energy input.
- `energy_output::Float64`: The energy output from the STES.
- `temperature_output::Temperature`: The temperature corresponding to the energy output.
- `sim_params::Dict{String,Any}`: A dictionary containing simulation parameters
- `temporary_calculation::Bool=false`: If `true`, the function performs a temporary calculation without considering losses

# Returns if temporary_calculation == true
- `t_new::Vector{Float64}`: The updated temperature distribution across the layers of the STES.
# Returns if temporary_calculation == false
- `losses::Float64`: The total energy losses during the update.
- `unit.load::Float64`: The previous load of the STES.
- `load_new::Float64`: The new load of the STES after the update.
- `t_new::Vector{Float64}`: The updated temperature distribution across the layers of the STES.
"""
function update_STES(unit::SeasonalThermalStorage,
                     energy_input::Vector{Float64},
                     temperatures_input::Vector{Temperature},
                     energy_output::Float64,
                     temperature_output::Temperature,
                     sim_params::Dict{String,Any};
                     temporary_calculation::Bool=false)
    if temporary_calculation
        consider_losses = 0
    else
        consider_losses = 1
    end

    # set alias for better readability
    t_old = unit.temperature_segments
    dt = sim_params["time_step_seconds"] / 60 / 60  # [h]
    number_of_internal_timesteps = 1  # number of internal time steps for the explicit Euler method
    t_new = Vector{Temperature}(zeros(unit.number_of_layer_total))

    # set lower node always to 1 for charging and discharging to always use the whole storage
    lower_node = 1

    # define input temperature for discharging
    return_temperature_input = unit.low_temperature

    # sort input by temperature, starting with the lowest, to avoid that energy can not be used for charging
    energy_input = energy_input[sortperm(temperatures_input)]
    sort!(temperatures_input)

    # mass_in and mass_out shoud be both posivive! 
    # mass_in and the corresponding temperature can be a vector and will be fitted to the most suitable layer.
    # mass_out and the corresponding temperature are single values that will always be drawn from the top layer
    mass_in = convert_energy_in_mass.(energy_input, t_old[lower_node], temperatures_input,
                                      unit.cp_medium)
    mass_out = convert_energy_in_mass(energy_output, unit.low_temperature, temperature_output,
                                      unit.cp_medium)

    # Check if the mass flow is greater than the volume of the smallest segment. 
    # If yes, the internal time step is reduced to avoid numerical instabilities.
    mass_out_sum = sum(mass_out; init=0.0)
    mass_in_sum = sum(mass_in; init=0.0)
    mass_of_smallest_segment = minimum(unit.volume_segments) * unit.rho_medium
    factor = 5
    if max(mass_out_sum, mass_in_sum) > (mass_of_smallest_segment / factor) && (mass_out_sum > 0.0 || mass_in_sum > 0.0)
        number_of_internal_timesteps = ceil(max(mass_out_sum, mass_in_sum) / (mass_of_smallest_segment / factor))
        if number_of_internal_timesteps > 10000
            number_of_internal_timesteps = 10000
            @warn "In STES $(unit.uac), a very high internal time step has been determined. " *
                  "Is was limited to 10 000.  This may lead to inaccurate results and increase " *
                  "the losses in the current time step $(sim_params["time"])."
        end
        dt = dt / number_of_internal_timesteps
    end

    for internal_timestep in 1:number_of_internal_timesteps
        if internal_timestep > 1
            # recalculate masses as the temperatures of the layers may have changed
            mass_in = convert_energy_in_mass.(energy_input, t_old[lower_node],
                                              temperatures_input,
                                              unit.cp_medium)
            mass_out = convert_energy_in_mass(energy_output, unit.low_temperature,
                                              t_old[unit.number_of_layer_total - unit.output_layer_from_top + 1],
                                              unit.cp_medium)
        end
        # mass flow and temperatures for charging
        mass_in_temp, mass_in_vec = calculate_mass_temperatur_charging(unit, t_old,
                                                                       mass_in ./ number_of_internal_timesteps,
                                                                       lower_node, sim_params)
        # mass flow and corresponding temperatures into each layer during discharging
        mass_out_temp, mass_out_vec = calculate_mass_temperatur_discharging(unit, t_old,
                                                                            mass_out ./
                                                                            number_of_internal_timesteps,
                                                                            return_temperature_input,
                                                                            lower_node, sim_params)
        for n in 1:(unit.number_of_layer_total)
            if n == 1  # bottom layer, single-side
                t_new[n] = t_old[n] +
                           consider_losses *
                           (3600 * unit.diffussion_coefficient * (t_old[n + 1] - t_old[n]) / unit.dz_normalized[n]^2 +   # thermal diffusion
                            consider_losses * unit.sigma[n] * (unit.effective_ambient_temperature[n] - t_old[n])) * dt + # losses through bottom and side walls
                           # unit.lambda[n] * (Q_in_out)[n] +                                                            # thermal input and output
                           unit.phi[n] * mass_in_vec[n] * (mass_in_temp[n] - t_old[n]) +                                 # mass input
                           unit.phi[n] * mass_out_vec[n] * (mass_out_temp[n] - t_old[n])                                 # mass output
            elseif n == unit.number_of_layer_total  # top layer, single-side
                t_new[n] = t_old[n] +
                           consider_losses *
                           (3600 * unit.diffussion_coefficient * (t_old[n - 1] - t_old[n]) / unit.dz_normalized[n]^2 +   # thermal diffusion
                            consider_losses * unit.sigma[n] * (unit.effective_ambient_temperature[n] - t_old[n])) * dt + # losses through lid and side walls
                           # unit.lambda[n] * Q_in_out[n] +                                                              # thermal input and output
                           unit.phi[n] * mass_in_vec[n] * (mass_in_temp[n] - t_old[n]) +                                 # mass input
                           unit.phi[n] * mass_out_vec[n] * (mass_out_temp[n] - t_old[n])                                 # mass output
            else       # mid layer
                t_new[n] = t_old[n] +
                           consider_losses *
                           (3600 * unit.diffussion_coefficient * (t_old[n + 1] + t_old[n - 1] - 2 * t_old[n]) /
                            unit.dz_normalized[n]^2 +                                                                    # thermal diffusion
                            consider_losses * unit.sigma[n] * (unit.effective_ambient_temperature[n] - t_old[n])) * dt + # losses through side walls
                           # unit.lambda[n] * Q_in_out[n] +                                                              # thermal input and output
                           unit.phi[n] * mass_in_vec[n] * (mass_in_temp[n] - t_old[n]) +                                 # mass input
                           unit.phi[n] * mass_out_vec[n] * (mass_out_temp[n] - t_old[n])                                 # mass output
            end

            if n > 1   # mixing due to buoancy effecs, if temperature gradient is present
                adjust_temp = max(0, t_new[n - 1] - t_new[n])
                t_new[n] = t_new[n] + unit.theta[n] * adjust_temp
                t_new[n - 1] = t_new[n - 1] - (1 - unit.theta[n]) * adjust_temp
            end
        end
        t_old = copy(t_new)
    end

    if temporary_calculation
        return t_new
    else
        load_new = (weighted_mean(t_new, unit.volume_segments, sim_params) - unit.low_temperature) /
                   (unit.high_temperature - unit.low_temperature) * unit.capacity
        losses = unit.load + sum(energy_input) - energy_output - load_new # + sum(Q_in_out)  # losses are positive here
        return losses, unit.load, load_new, t_new
    end
end

"""
    calculate_mass_temperatur_charging(unit::SeasonalThermalStorage, t_old::Vector{Temperature},
                                       mass_in::Vector{Float64}, lower_node::Int, sim_params::Dict{String,Any})

Calculate the mass flow and its temperature into each single layer of a Seasonal Thermal Energy Storage (STES) during charging.
It handles cases where the mass flow into a layer is greater than the volume of the layer as well as multiple inputs
with different temperatures. The function is based on the ideal charging system where the mass flow is put into the
lowest layer that is below the temperature of the input (e.g. lance)

# Arguments
- `unit::SeasonalThermalStorage`: The seasonal thermal storage unit.
- `t_old::Vector{Temperature}`: Vector of temperatures in each layer before charging.
- `mass_in::Vector{Float64}`: Vector of mass inputs for each input source.
- `lower_node::Int`: The index of the lowest layer to consider for charging.
- `sim_params::Dict{String,Any}`: Dictionary of simulation parameters.

# Returns
- `mass_in_temp::Vector{Float64}`: Vector of temperatures of the mass flow into each layer.
- `mass_in_vec::Vector{Float64}`: Vector of mass flow into each layer.
"""
function calculate_mass_temperatur_charging(unit::SeasonalThermalStorage, t_old::Vector{Temperature},
                                            mass_in::Vector{Float64}, lower_node::Int, sim_params::Dict{String,Any})
    number_of_inputs = length(unit.current_energy_input)
    mass_in_temp = zeros(unit.number_of_layer_total)
    mass_in_vec = zeros(unit.number_of_layer_total)
    upper_node_charging = Int[]
    mass_input = zeros(unit.number_of_layer_total, number_of_inputs)
    temperature_input = zeros(unit.number_of_layer_total, number_of_inputs)

    # find index of uppermost layer that is below the temperature, for each input
    # This represents a ideal charging system (lance e.g.)
    for input in 1:number_of_inputs
        for layer in (unit.number_of_layer_total):-1:1
            if t_old[layer] < unit.current_temperature_input[input]
                push!(upper_node_charging, layer)
                mass_input[layer, input] = mass_in[input]
                temperature_input[layer, input] = unit.current_temperature_input[input]
                break
            end
        end
    end

    # calculate masses and corresponding temperatures flowing into each layer during charging
    mass_in_layer = Float64[]
    temperature_in_layer = Float64[]
    for layer in (unit.number_of_layer_total):-1:lower_node
        # add external mass and temperature input to the upcoming layer
        if layer in upper_node_charging
            append!(mass_in_layer, mass_input[layer, :])
            append!(temperature_in_layer, temperature_input[layer, :])
        end

        # write mass input to vector
        mass_in_vec[layer] = min(sum(mass_in_layer), unit.layer_masses[layer])

        # calculate tempearture of mass flow into the layer due to discharging
        if sum(mass_in_layer) <= unit.layer_masses[layer]  # all input mass can be hold in current layer
            mass_in_temp[layer] = weighted_mean(temperature_in_layer, mass_in_layer, sim_params)

            # set input for next layer
            temperature_in_layer = [t_old[layer]]
            mass_in_layer = [sum(mass_in_layer)]
        else # here, the mass put into the layer is bigger than the volume of the layer
            mass_to_keep = Float64[]
            temperature_to_keep = Float64[]
            while sum(mass_to_keep) < unit.layer_masses[layer]
                needed = unit.layer_masses[layer] - sum(mass_to_keep)
                if mass_in_layer[end] < needed
                    push!(mass_to_keep, mass_in_layer[end])
                    push!(temperature_to_keep, temperature_in_layer[end])
                    pop!(mass_in_layer)
                    pop!(temperature_in_layer)
                else
                    push!(mass_to_keep, needed)
                    push!(temperature_to_keep, temperature_in_layer[end])
                    mass_in_layer[end] -= needed
                end
            end
            mass_in_temp[layer] = weighted_mean(temperature_to_keep, mass_to_keep, sim_params)

            # add volume of layer as input for next layer
            pushfirst!(mass_in_layer, unit.layer_masses[layer])
            pushfirst!(temperature_in_layer, t_old[layer])
        end
    end
    return mass_in_temp, mass_in_vec
end

"""
    calculate_mass_temperatur_discharging(unit::SeasonalThermalStorage, t_old::Vector{Temperature},
                                         mass_out::Float64, return_temperature_input::Temperature, lower_node::Int,
                                         sim_params::Dict{String,Any})

Calculate the mass flow and its temperature into each single layer of a STES during discharging.
It handles cases where the mass flow into a layer is greater than the volume of the layer.

# Arguments
- `unit::SeasonalThermalStorage`: The seasonal thermal storage unit.
- `t_old::Vector{Temperature}`: Vector of temperatures in each layer before discharging.
- `mass_out::Float64`: The mass flow rate out of the storage.
- `return_temperature_input::Temperature`: The return temperature input to the storage at lower_node.
- `lower_node::Int`: The index of the lower node where discharging return input is put into the storage.
- `sim_params::Dict{String,Any}`: Simulation parameters.

# Returns
- `mass_out_temp::Vector{Float64}`: Vector of temperatures of the mass flow into each layer during discharging.
- `mass_out_vec::Vector{Float64}`: Vector of mass flow rates into each layer during discharging

"""
function calculate_mass_temperatur_discharging(unit::SeasonalThermalStorage, t_old::Vector{Temperature},
                                               mass_out::Float64, return_temperature_input::Temperature,
                                               lower_node::Int, sim_params::Dict{String,Any})
    upper_node_discharging = unit.number_of_layer_total - unit.output_layer_from_top + 1
    mass_out_temp = zeros(unit.number_of_layer_total)
    mass_in_layer = [mass_out]
    temperature_in_layer = [return_temperature_input]
    for layer in lower_node:upper_node_discharging
        # calculate tempearture of mass flow into the layer due to discharging
        if sum(mass_in_layer) <= unit.layer_masses[layer]
            if layer == lower_node
                mass_out_temp[layer] = return_temperature_input
            else
                mass_out_temp[layer] = weighted_mean(temperature_in_layer, mass_in_layer, sim_params)
            end
            # set input for next layer
            temperature_in_layer = [t_old[layer]]
            mass_in_layer = [mass_out]
        else # here, the mass put into the layer is bigger than the volume of the layer
            if layer == lower_node
                mass_out_temp[layer] = return_temperature_input
                # set input for next layer
                mass_in_layer = [unit.layer_masses[layer], mass_in_layer[1] - unit.layer_masses[layer]]
                temperature_in_layer = [t_old[layer], temperature_in_layer[1]]
            else
                mass_to_keep = Float64[]
                temperature_to_keep = Float64[]
                while sum(mass_to_keep) < unit.layer_masses[layer]
                    needed = unit.layer_masses[layer] - sum(mass_to_keep)
                    if mass_in_layer[end] < needed
                        push!(mass_to_keep, mass_in_layer[end])
                        push!(temperature_to_keep, temperature_in_layer[end])
                        pop!(mass_in_layer)
                        pop!(temperature_in_layer)
                    else
                        push!(mass_to_keep, needed)
                        push!(temperature_to_keep, temperature_in_layer[end])
                        mass_in_layer[end] -= needed
                    end
                end
                mass_out_temp[layer] = weighted_mean(temperature_to_keep, mass_to_keep, sim_params)

                # set input for next layer
                pushfirst!(mass_in_layer, unit.layer_masses[layer])
                pushfirst!(temperature_in_layer, t_old[layer])
            end
        end
    end
    mass_out_vec = [((lower_node <= i <= upper_node_discharging) ? min(mass_out, unit.layer_masses[i]) : 0.0)
                    for i in 1:(unit.number_of_layer_total)]

    return mass_out_temp, mass_out_vec
end

function process(unit::SeasonalThermalStorage, sim_params::Dict{String,Any})
    outface = unit.output_interfaces[unit.m_heat_out]
    exchanges = balance_on(outface, outface.target)
    energy_demanded = balance(exchanges) + energy_potential(exchanges)
    energy_available = unit.current_max_output_energy  # is positive

    # shortcut if there is no energy demanded
    if energy_demanded >= -sim_params["epsilon"]
        set_max_energy!(unit.output_interfaces[unit.m_heat_out], 0.0)
        return
    end

    for exchange in exchanges
        demanded_on_interface = exchange.balance + exchange.energy_potential
        if demanded_on_interface >= -sim_params["epsilon"]
            continue
        end

        if exchange.temperature_min !== nothing &&
           exchange.temperature_min > unit.current_max_output_temperature
            # we can only supply energy at a temperature at or below the STES's current
            # max output temperature
            continue
        end

        used_heat = abs(demanded_on_interface)

        if energy_available > used_heat
            energy_available -= used_heat
            add!(outface, used_heat, nothing, unit.current_max_output_temperature, exchange.purpose_uac)
            unit.current_energy_output += used_heat
        else
            add!(outface, energy_available, nothing, unit.current_max_output_temperature, exchange.purpose_uac)
            unit.current_energy_output += energy_available
            energy_available = 0.0
        end
        unit.current_temperature_output = unit.current_max_output_temperature
    end
end

function load(unit::SeasonalThermalStorage, sim_params::Dict{String,Any})
    inface = unit.input_interfaces[unit.m_heat_in]
    exchanges = balance_on(inface, inface.source)
    energy_available = balance(exchanges) + energy_potential(exchanges)

    # shortcut if there is no energy to be used
    if energy_available <= sim_params["epsilon"]
        set_max_energy!(unit.input_interfaces[unit.m_heat_in], 0.0)
        unit.losses,
        unit.load_end_of_last_timestep,
        unit.load,
        unit.temperature_segments = update_STES(unit,
                                                unit.current_energy_input,
                                                unit.current_temperature_input,
                                                unit.current_energy_output,
                                                unit.current_temperature_output,
                                                sim_params)
        # set new temperatures limits for the next time step
        set_temperature_limits!(unit, sim_params)
        return
    end

    for exchange in exchanges
        exchange_energy_available = exchange.balance + exchange.energy_potential
        if exchange_energy_available < sim_params["epsilon"]
            continue
        end

        if (exchange.temperature_max !== nothing &&
            exchange.temperature_max < unit.current_min_input_temperature)
            # we can only take in energy if it's at a higher/equal temperature than the
            # STES's second-coldest layer
            continue
        end

        current_exchange_temperature = highest(exchange.temperature_max, unit.current_min_input_temperature)

        # all energies in the exchange can be used as it was already made sure than they do not exeed the STES limit
        sub!(inface, exchange_energy_available, current_exchange_temperature, nothing, exchange.purpose_uac)
        push!(unit.current_energy_input, exchange_energy_available)
        push!(unit.current_temperature_input, current_exchange_temperature)
    end

    # apply energy input and output and losses
    unit.losses,
    unit.load_end_of_last_timestep,
    unit.load,
    unit.temperature_segments = update_STES(unit,
                                            unit.current_energy_input,
                                            unit.current_temperature_input,
                                            unit.current_energy_output,
                                            unit.current_temperature_output,
                                            sim_params)
    # set new temperatures limits for the next time step
    set_temperature_limits!(unit, sim_params)
end

function reset(unit::SeasonalThermalStorage)
    invoke(reset, Tuple{Component}, unit)

    # reset variables
    unit.current_energy_input = []
    unit.current_temperature_input = []
    unit.current_energy_output = 0.0
    unit.current_temperature_output = 0.0
    unit.losses = 0.0
end

function output_values(unit::SeasonalThermalStorage)::Vector{String}
    return [string(unit.m_heat_in) * " IN",
            string(unit.m_heat_out) * " OUT",
            "Load",
            "Load%",
            "Capacity",
            "Losses",
            "CurrentMaxOutTemp"]
end

function output_value(unit::SeasonalThermalStorage, key::OutputKey)::Float64
    if key.value_key == "IN"
        return calculate_energy_flow(unit.input_interfaces[key.medium])
    elseif key.value_key == "OUT"
        return calculate_energy_flow(unit.output_interfaces[key.medium])
    elseif key.value_key == "Load"
        return unit.load
    elseif key.value_key == "Load%"
        return 100 * unit.load / unit.capacity
    elseif key.value_key == "Capacity"
        return unit.capacity
    elseif key.value_key == "Losses"
        return unit.losses
    elseif key.value_key == "CurrentMaxOutTemp"
        return unit.current_max_output_temperature
    end
    throw(KeyError(key.value_key))
end

export SeasonalThermalStorage
