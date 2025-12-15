"""
Implementation of a seasonal thermal storage component.

This is a simplified model, which mostly deals with amounts of energy and considers
temperatures only for the available temperature as the tank is depleted.
"""

using Plots: Plots
using Dates
using PlotlyJS: PlotlyJS
using SparseArrays
using LinearAlgebra
using GLMakie

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
    ground_model::String
    rho_medium::Float64
    cp_medium::Float64
    diffusion_coefficient::Float64
    surface_area_lid::Float64
    surface_area_bottom::Float64
    surface_area_barrel_segments::Vector{Float64}
    volume_segments::Vector{Float64}
    height::Float64
    radius_small::Float64
    radius_large::Float64
    number_of_layer_total::Int64
    number_of_layer_above_ground::Int64
    number_of_STES_layer_below_ground::Int64
    h_stes_buried::Float64
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
    max_load_rate_energy::Floathing
    max_unload_rate_energy::Floathing
    max_load_rate_mass::Floathing
    max_unload_rate_mass::Floathing
    max_input_energy::Floathing
    max_output_energy::Floathing

    # losses
    thermal_transmission_lid::Float64
    thermal_transmission_barrel_above_ground::Float64
    thermal_transmission_barrel_below_ground::Float64
    thermal_transmission_barrels::Vector{Float64}
    thermal_transmission_bottom::Float64
    ambient_temperature_profile::Union{Profile,Nothing}
    ambient_temperature::Temperature
    ground_temperature_profile::Union{Profile,Nothing}
    ground_temperature::Temperature
    effective_ambient_temperature_barrels::Vector{Temperature}
    effective_ambient_temperature_top::Temperature
    effective_ambient_temperature_bottom::Temperature

    # ground coupling FEM (optional)
    ground_domain_radius_factor::Float64
    ground_domain_depth_factor::Float64
    ground_domain_radius::Floathing
    ground_domain_depth::Floathing
    ground_accuracy_mode::String
    ground_layers_depths::Vector{Floathing}
    ground_layers_k::Vector{Float64}
    ground_layers_rho::Vector{Float64}
    ground_layers_cp::Vector{Float64}
    row_k::Vector{Float64}
    row_rho::Vector{Float64}
    row_cp::Vector{Float64}
    soil_surface_hconv::Float64

    # FEM state (unified axisymmetric r-z soil domain)
    soil_dr::Vector{Float64}
    soil_dz::Vector{Float64}
    soil_dr_mesh::Vector{Float64}
    soil_dz_mesh::Vector{Float64}
    soil_r_centers::Vector{Float64}
    soil_z_centers::Vector{Float64}
    soil_t1::Array{Float64}
    soil_t2::Array{Float64}
    cells_active::Matrix{Bool}
    radius_at_row::Array{Float64}
    equivalent_radius_from_bottom::Float64

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
    current_energy_input_return_temperature::Float64
    process_done::Bool
    load_done::Bool

    # additional output
    temp_distribution_output::Array{Float64}
    temp_difference_to_surrounding_output::Array{Float64}
    soil_temperature_field_output::Array{Float64}

    mass_in_sum::Float64
    mass_out_sum::Float64

    epsilon_geometry::Float64

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

        constant_ground_temperature,
        ground_temperature_profile = get_parameter_profile_from_config(config,
                                                                       sim_params,
                                                                       "ground_temperature",
                                                                       "ground_temperature_profile_file_path",
                                                                       "",
                                                                       "constant_ground_temperature",
                                                                       uac;
                                                                       required=true)

        # Note: layer numbering within the STES starts at the bottom with index 1 and ends at the top with index number_of_layer_total
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
                   default(config, "sidewall_angle", 40.0),        # angle of the sidewall of the STES with respect to the horizon [°]
                   default(config, "shape", "quadratic"),          # can be "round" for cylinder/truncated cone or "quadratic" for tank or truncated quadratic pyramid (pit)
                   default(config, "ground_model", "FEM"),         # ground_model. Can be one of "simple" or "FEM"
                   default(config, "rho_medium", 1000.0),          # density of the medium [kg/m^3]
                   default(config, "cp_medium", 4.186),            # specific thermal capacity of medium [kJ/kgK]
                   default(config, "diffusion_coefficient", 0.143 * 10^-6), # diffusion coefficient of the medium [m^2/s]
                   0.0,                                            # surface_area_lid, surface of the lid of the STES [m^2]
                   0.0,                                            # surface_area_bottom, surface of the bottom of the STES [m^2]
                   Float64[],                                      # surface_area_barrel_segments, surface of the barrel segments of the STES [m^2]
                   Float64[],                                      # volume_segments, volume of the segments of the STES [m^3]
                   0.0,                                            # height of the STES [m]
                   0.0,                                            # radius_small, small (lower) radius of the STES [m]
                   0.0,                                            # radius_large, large (upper) radius of the STES [m]
                   default(config, "number_of_layer_total", 25),   # number of layers in the STES
                   default(config, "number_of_layer_above_ground", 5), # number of layers above ground in the STES
                   0,                                              # number_of_STES_layer_below_ground
                   0.0,                                            # h_stes_buried
                   default(config, "output_layer_from_top", 1),    # layer number of the output layer, counted from the top
                   Float64[],                                      # dz, thickness of the layers of the STES [m]
                   Float64[],                                      # dz_normalized, normalized dz with respect to to the volume of each section
                   Float64[],                                      # [1/h]  factor for losses to ambient: area_of_losses * U[kJ/m^2K] / (roh * cp * volume_segment)
                   Float64[],                                      # [K/kJ] factor for input/output energy:  1 / (roh * cp * volume_segment ) 
                   Float64[],                                      # [1/kg] factor for input/output mass flow:  1 / (roh  * volume_segment )
                   Float64[],                                      # volume-ratios of sections: V_section[n-1] / (V_section[n] + V_section[n-1])
                   Float64[],                                      # layer_masses, mass of the medium in each layer [kg]
                   # loading and unloading
                   default(config, "high_temperature", 90.0),          # upper temperature of the STES [°C]
                   default(config, "low_temperature", 15.0),           # lower temperature of the STES [°C]
                   default(config, "max_load_rate_energy", nothing),   # maximum load rate given in 1/h, total energy input related to the total energy capacity of the STES
                   default(config, "max_unload_rate_energy", nothing), # maximum unload rate given in 1/h, total energy input related to the total volume of the STES
                   default(config, "max_load_rate_mass", nothing),     # maximum load rate given in 1/h, mass flow per Interface related to the total volume of the STES
                   default(config, "max_unload_rate_mass", nothing),   # maximum unload rate given in 1/h, mass flow in total related to the total volume of the STES
                   nothing,                                            # max_input_energy, maximum input energy per time step [Wh]
                   nothing,                                            # max_output_energy, maximum output energy per time step [Wh]
                   # Losses
                   default(config, "thermal_transmission_lid", 0.25),                  # [W/(m^2K)]
                   default(config, "thermal_transmission_barrel_above_ground", 0.375), # [W/(m^2K)]
                   default(config, "thermal_transmission_barrel_below_ground", 12.0),  # [W/(m^2K)]
                   Float64[],                                                          # thermal_transmission_barrels
                   default(config, "thermal_transmission_bottom", 0.375),              # [W/(m^2K)]
                   ambient_temperature_profile,                           # [°C]
                   constant_ambient_temperature,                          # ambient_temperature [°C]
                   ground_temperature_profile,                            # [°C]
                   constant_ground_temperature,                           # ground_temperature [°C]
                   Temperature[],                                         # effective_ambient_temperature_barrels corresponding to each layer [°C]          
                   0.0,                                                   # effective_ambient_temperature_top 
                   0.0,                                                   # effective_ambient_temperature_bottom
                   # ground coupling FEM (unified)
                   default(config, "ground_domain_radius_factor", 1.5),          # [m] ground_domain_radius_factor: Factor for the ground domain width, is multiplied with the radius of the storage at the ground surface.
                   default(config, "ground_domain_depth_factor", 2.0),           # [m] ground_domain_depth_factor: Factor for the ground domain depth, is multiplied with the total height of the storage.
                   default(config, "ground_domain_radius", nothing),             # [m] soil domain radius. If none given, it will be derived from the STES geometry.
                   default(config, "ground_domain_depth", nothing),              # [m] soil depth from surface. If none given, it will be derived from the STES geometry.
                   default(config, "ground_accuracy_mode", "normal"),            # mesh preset: very_rough|rough|normal|high|very_high
                   default(config, "ground_layers_depths", [nothing]),           # [m] Monotonically increasing depths from surface (0.0) downward.
                   # Defines piecewise-constant soil layers by intervals [d[i], d[i+1]).
                   # Last value should be ≥ ground_domain_depth; if shorter, code extends.
                   default(config, "ground_layers_k", Float64[1.5]),            # [W/(m·K)] Thermal conductivity per layer. Length may be
                   # ≤ (length(depths)-1); last value is repeated if shorter.
                   default(config, "ground_layers_rho", Float64[2000.0]),       # [kg/m³] Mass density per layer. Same broadcasting rule as above.
                   default(config, "ground_layers_cp", Float64[1000.0]),        # [J/(kg·K)] Specific heat capacity per layer 
                   Float64[],                                                   # row_k: Thermal conductivity per row
                   Float64[],                                                   # row_rho: Mass density per row
                   Float64[],                                                   # row_cp: Specific heat capacity per row
                   default(config, "soil_surface_hconv", 14.7),                 # [W/(m²·K)] ground surface convective heat transfer coefficient

                   # FEM state (unified axisymmetric r-z soil domain)
                   Float64[],                                  # soil_dr
                   Float64[],                                  # soil_dz
                   Float64[],                                  # soil_dr_mesh
                   Float64[],                                  # soil_dz_mesh
                   Float64[],                                  # soil_r_centers
                   Float64[],                                  # soil_z_centers
                   Array{Float64}(undef, 0, 0),                # soil_t1
                   Array{Float64}(undef, 0, 0),                # soil_t2
                   Matrix{Bool}(undef, 0, 0),                  # cells_active
                   Float64[],                                  # radius_at_row
                   0.0,                                        # equivalent_radius_from_bottom

                   # state variables
                   0.0,                                            # current_max_output_temperature
                   0.0,                                            # current_min_input_temperature
                   default(config, "initial_load", 0.0),           # initial_load [%/100] assuming perfectly mixed storage at the begin
                   0.0,                                            # load, set to initial_load at the beginning [Wh]
                   0.0,                                            # load_end_of_last_timestep, stores the load of the previous time step without losses
                   0.0,                                            # losses in current time step [Wh]
                   Float64[],                                      # temperature_segments: temperatures of the segments
                   Float64[],                                      # current_energy_input, energy input in current time step [Wh]
                   Temperature[],                                  # current_temperature_input, temperature of the input in current time step [°C]
                   0.0,                                            # current_energy_output, energy output in current time step [Wh]
                   0.0,                                            # current_temperature_output, temperature of the output in current time step [°C]
                   0.0,                                            # current_max_output_energy, maximum output energy in current time step [Wh]
                   Float64[],                                      # temperatures_charging, temperatures of the possible inputs from exchange in current time step [°C]
                   0.0,                                            # current_energy_input_return_temperature, current return temperature for the energy input for control module LimitCoolingInputTemperature
                   false,                                          # process_done, bool indicating if the process step has already been performed in the current time step
                   false,                                          # load_done, bool indicating if the load step has already been performed in the current time step
                   # additional output
                   Array{Float64}(undef, 0, 0),                    # temp_distribution_output [°C], holds temperature field of layers for output plot
                   Array{Float64}(undef, 0, 0),                    # temp_difference_to_surrounding_output [°C], holds temperature difference of storage temperature and surroundings for output plot
                   Array{Float64}(undef, 0, 0, 0),                 # soil_temperature_field_output  (time × nz × nr)
                   0.0,                                            # mass_in_sum [kg]
                   0.0,                                            # mass_out_sum [kg]
                   1e-9)                                           # epsilon_geometry
    end
end

function initialise!(unit::SeasonalThermalStorage, sim_params::Dict{String,Any})
    # Hook up input/output flow control
    set_storage_transfer!(unit.input_interfaces[unit.m_heat_in],
                          unload_storages(unit.controller, unit.m_heat_in))
    set_storage_transfer!(unit.output_interfaces[unit.m_heat_out],
                          load_storages(unit.controller, unit.m_heat_out))

    # set temperature vector: assuming a uniform temperature profile (mixed storage)
    mean_temperature = unit.initial_load * (unit.high_temperature - unit.low_temperature) + unit.low_temperature
    unit.temperature_segments = [mean_temperature for _ in 1:(unit.number_of_layer_total)]

    # set initial temperature bounds
    set_temperature_limits!(unit, sim_params)

    # set initial current_energy_input_return_temperature, always use bottom layer to use the whole storage
    unit.current_energy_input_return_temperature = unit.temperature_segments[1]

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

    unit.number_of_STES_layer_below_ground = unit.number_of_layer_total - unit.number_of_layer_above_ground
    unit.h_stes_buried = unit.height / unit.number_of_layer_total * unit.number_of_STES_layer_below_ground

    # get and check soil boundaries
    storage_radius_surface = unit.radius_large -
                             (unit.radius_large - unit.radius_small) / unit.number_of_layer_total *
                             unit.number_of_layer_above_ground
    if unit.ground_domain_radius === nothing
        unit.ground_domain_radius = unit.ground_domain_radius_factor * storage_radius_surface
    elseif unit.ground_domain_radius <= storage_radius_surface
        @error "In STES $(unit.uac), the given ground_domain_radius has to be greater than the radius of the " *
               "STES at the surface which is $(storage_radius_surface) m for the current configuration."
    end
    if unit.ground_domain_depth === nothing
        unit.ground_domain_depth = unit.ground_domain_depth_factor * unit.height
    elseif unit.ground_domain_depth <= unit.h_stes_buried
        @error "In STES $(unit.uac), the given ground_domain_depth has to be greater than the height of buried " *
               "depth of the STES which is $(unit.h_stes_buried) m for the current configuration."
    end
    if unit.ground_layers_depths == [nothing]
        unit.ground_layers_depths = [0, unit.ground_domain_depth]
    end

    # calculate thermal transmission coefficients
    unit.thermal_transmission_barrels = zeros(unit.number_of_layer_total)
    for layer in 1:(unit.number_of_layer_total)
        if layer <= unit.number_of_STES_layer_below_ground
            unit.thermal_transmission_barrels[layer] = unit.thermal_transmission_barrel_below_ground
        else
            unit.thermal_transmission_barrels[layer] = unit.thermal_transmission_barrel_above_ground
        end
    end

    # calculate the mass of the medium in each layer
    unit.layer_masses = unit.rho_medium .* unit.volume_segments

    # calculate (equally spaced) thickness of the layers of the STES [m]
    unit.dz = fill(unit.height / unit.number_of_layer_total, unit.number_of_layer_total)
    # calculate the normalized dz with respect to to the volume of each section
    unit.dz_normalized = unit.dz .* sqrt.((unit.volume_segments .* unit.number_of_layer_total) ./ unit.volume)

    # calculate coefficient for losses to ambient
    unit.sigma = zeros(unit.number_of_layer_total)
    unit.sigma[1] = unit.surface_area_bottom * unit.thermal_transmission_bottom /
                    (unit.rho_medium * convert_kJ_in_Wh(unit.cp_medium) * unit.volume_segments[1])          # [1/h] losses to ambient through bottom
    unit.sigma[end] = unit.surface_area_lid * unit.thermal_transmission_lid /
                      (unit.rho_medium * convert_kJ_in_Wh(unit.cp_medium) * unit.volume_segments[end])      # [1/h] losses to ambient through lid
    unit.sigma = unit.sigma .+
                 unit.surface_area_barrel_segments .* unit.thermal_transmission_barrels ./
                 (unit.rho_medium * convert_kJ_in_Wh(unit.cp_medium) * unit.volume_segments)                # [1/h]  losses to ambient though barrel

    # calculate coefficient for input/output energy. 
    # Currently not used, as no direct loading of energy into the storage through heat exchangers is implemented.
    # But may this is useful in the future...
    # unit.lambda = 1 ./ (unit.rho_medium * convert_kJ_in_Wh(unit.cp_medium) * unit.volume_segments)          # [K/Wh] 

    # coefficient for input/output mass flow, assuming water as fluid
    cp_water = 4.18                                                                      # [kJ/kgK]
    unit.phi = cp_water ./ (unit.cp_medium * unit.rho_medium * unit.volume_segments)     # [1/kg]

    # coefficient for buoyancy effects
    unit.theta = [unit.volume_segments[n - 1] / (unit.volume_segments[n] + unit.volume_segments[n - 1])
                  for n in 2:(unit.number_of_layer_total)]
    pushfirst!(unit.theta, 0.0)  # Set first element to 0

    unit.capacity = unit.volume * unit.rho_medium * convert_kJ_in_Wh(unit.cp_medium) *
                    (unit.high_temperature - unit.low_temperature)  # [Wh]
    unit.load = unit.initial_load * unit.capacity
    unit.load_end_of_last_timestep = copy(unit.load)

    # calculate maximum input and output energy
    if unit.max_load_rate_energy === nothing
        unit.max_input_energy = Inf
    else
        unit.max_input_energy = unit.max_load_rate_energy * unit.capacity * (sim_params["time_step_seconds"] / 60 / 60)     # [Wh per timestep]
    end
    if unit.max_unload_rate_energy === nothing
        unit.max_output_energy = Inf
    else
        unit.max_output_energy = unit.max_unload_rate_energy * unit.capacity *
                                 (sim_params["time_step_seconds"] / 60 / 60)  # [Wh per timestep]
    end

    # set initial ambient/ground boundary conditions
    if unit.ambient_temperature_profile !== nothing
        unit.ambient_temperature = Profiles.value_at_time(unit.ambient_temperature_profile, sim_params)
    end
    if unit.ground_temperature_profile !== nothing
        unit.ground_temperature = Profiles.value_at_time(unit.ground_temperature_profile, sim_params)
    end

    # vector to hold the results of the temperatures for each layer in each simulation time step
    unit.temp_distribution_output = zeros(Float64, sim_params["number_of_time_steps_output"],
                                          unit.number_of_layer_total)
    # top and bottom get their own surrounding temperature
    unit.temp_difference_to_surrounding_output = zeros(Float64, sim_params["number_of_time_steps_output"],
                                                       unit.number_of_layer_total + 2)

    # set initial effective_ambient_temperature
    unit.effective_ambient_temperature_barrels = [i <= unit.number_of_STES_layer_below_ground ?
                                                  unit.ground_temperature : unit.ambient_temperature
                                                  for i in 1:(unit.number_of_layer_total)]
    unit.effective_ambient_temperature_top = unit.ambient_temperature
    unit.effective_ambient_temperature_bottom = unit.ground_temperature

    if unit.ground_model == "FEM"
        # Prepare unified (r,z) soil FEM and allocate its output field
        prepare_ground_fem_unified!(unit)
        nz = length(unit.soil_dz)
        nr = length(unit.soil_dr)
        unit.soil_temperature_field_output = zeros(Float64,
                                                   sim_params["number_of_time_steps_output"],
                                                   nz, nr)
        # get soil properties per row
        unit.row_k = zeros(Float64, nz)
        unit.row_rho = zeros(Float64, nz)
        unit.row_cp = zeros(Float64, nz)
        for h in 1:nz
            unit.row_k[h], unit.row_rho[h], unit.row_cp[h] = soil_props_at_depth(unit, unit.soil_z_centers[h])
        end

        # calculate radius for each row
        unit.radius_at_row = zeros(Float64, nz)
        for h in 1:nz
            unit.radius_at_row[h] = radius_at_row(unit, h)
        end

        # get mask for storage area
        # true → this cell is *soil* and part of the PDE
        # false → this cell lies inside the STES volume (masked, handled separately)
        unit.cells_active = fill(false, nz, nr)
        for h in 1:nz, i in 1:nr
            unit.cells_active[h, i] = cell_active(unit, h, i)
        end

        # set equivalent radius from bottom
        unit.equivalent_radius_from_bottom, _ = equiv_radii_for_ground(unit)
    elseif unit.ground_model == "simple"
        # do nothing here
    else
        @error "In STES $(unit.uac), the given ground_model has to be either `simple` to use a constant surrounding " *
               "ground temperature or `FEM` to use the FEM-model to model the surrounding soil."
    end
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
- tank with quadratic base or a truncated quadratic pyramid (pit) (shape == quadratic)

# Arguments
- `uac::String`: Unique identifier for the STES
- `volume::Float64`: Total volume of the storage [m^3]
- `alpha::Float64`: Slope angle of the truncated cone in degrees with respect to the horizontal. 
If `alpha` is 90, the storage is a cylinder.
- `hr::Float64`: Height-to-radius ratio
- `n_segments::Int64`: Number of segments to divide the storage into for calculation
- `shape::String`: Shape of the STES, can be "round" for cylinder/truncated cone or "quadratic" for tank or 
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
    elseif alpha == 90 && shape == "quadratic"  # cuboid with square cross-section 
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
    elseif shape == "quadratic"   # truncated quadratic pyramid (pit)
        alpha_rad = deg2rad(alpha)
        alpha_tan = tan(alpha_rad)

        # check if input parameters are creating a valid truncated quadratic pyramid
        hr_max = 2 * alpha_tan
        if hr > hr_max
            alpha_min_deg = rad2deg(atan(hr / 2))
            alpha_max_deg = 180 - alpha_min_deg

            @error "For the STES $(uac), reduce the h/r ratio for the truncated quadratic pyramid or increase the slope angle. " *
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

            # linear Interpolation of edge length
            a_bot = base_side + (top_side - base_side) * (z_bot / height)
            a_up = base_side + (top_side - base_side) * (z_top / height)
            # volume of segments
            v_section[i] = (dh / 3) * (a_bot^2 + a_bot * a_up + a_up^2)

            # Segment-lateral area
            l_side_segment = sqrt(dh^2 + (a_up - a_bot)^2 / 2)
            h_trapeze = sqrt(l_side_segment^2 - ((a_up - a_bot) / 2)^2)
            a_lat[i] = 2 * (a_bot + a_up) * h_trapeze
        end
        return a_top, a_lat, a_bottom, v_section, height, base_side / 2, top_side / 2
    else
        @error "Invalid shape type of seasonal thermal storage $(uac). Shape has to be 'round' or 'quadratic'!"
        throw(InputError)
    end
end

function plot_optional_figures_begin(unit::SeasonalThermalStorage, output_path::String, output_formats::Vector{String},
                                     sim_params::Dict{String,Any})::Bool
    # Plot geometry of STES: 2D Cross Section 
    p = Plots.plot([-unit.radius_large, unit.radius_large], [unit.height, unit.height]; color=:blue, lw=6, label="")  # Top
    Plots.plot!(p, [-unit.radius_small, unit.radius_small], [0, 0]; color=:blue, lw=6, label="")  # Bottom
    Plots.plot!(p, [unit.radius_small, unit.radius_large], [0, unit.height]; color=:blue, lw=6, label="")  # Side wall
    Plots.plot!(p, [-unit.radius_small, -unit.radius_large], [0, unit.height]; color=:blue, lw=6, label="")  # Side wall
    Plots.plot!(p,; title="Cross section of the STES $(unit.uac) (cross-section: $(unit.shape))",
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
        Plots.savefig(p, output_path * "/" * fig_name * "." * output_format)
    end

    # Plot geometry of STES: 3D
    Plots.plotlyjs()
    if unit.shape == "quadratic"
        # — Cuboid setup — 
        r1 = unit.radius_small    # 1/2 of the bottom side length
        r2 = unit.radius_large    # 1/2 of the top side length
        H = unit.height           # height
        nh = 40     # vertical resolution
        ns = 30     # side‐to‐side resolution

        # interpolation of half‐width at each z
        half_width(z) = r1 + (r2 - r1) * (z / H)

        # Bottom face (square at z=0)
        Xb = [-r1 -r1;
              r1 r1]
        Yb = [-r1 r1;
              -r1 r1]
        Zb = zeros(2, 2)
        p = Plots.surface(Xb, Yb, Zb; color=:blue, alpha=0.65, legend=false)

        h_vec = range(0, H; length=nh)          # z levels
        s = range(-1, 1; length=ns)             # normalized side parameter

        # Four trapezoidal side faces
        for (sign, axis) in ((+1, :y), (-1, :y), (+1, :x), (-1, :x))
            # Build Z and S meshes
            Z = repeat(h_vec', ns, 1)
            S = repeat(s, 1, nh)
            # half‐width at each z‐row
            W = half_width.(Z)
            # parametric face
            if axis == :y
                X_face = S .* W
                Y_face = sign .* W
            else  # axis == :x
                X_face = sign .* W
                Y_face = S .* W
            end
            Plots.surface!(p, X_face, Y_face, Z; color=:blue, alpha=0.5, legend=false)
        end

        # add black contour lines
        edges = [
                 # bottom square
                 ([-r1, r1, r1, -r1, -r1],
                  [-r1, -r1, r1, r1, -r1],
                  zeros(5)),
                 # top square
                 ([-r2, r2, r2, -r2, -r2],
                  [-r2, -r2, r2, r2, -r2],
                  fill(H, 5)),
                 # four corners
                 ([-r1, -r2], [-r1, -r2], [0, H]), # back‐left
                 ([r1, r2], [-r1, -r2], [0, H]),   # back‐right
                 ([r1, r2], [r1, r2], [0, H]),     # front‐right
                 ([-r1, -r2], [r1, r2], [0, H])]   # front‐left

        for (xs, ys, zs) in edges
            Plots.plot3d!(p, xs, ys, zs; color=:black, lw=2)
        end

    elseif unit.shape == "round"
        # — Truncated cone (frustum) setup —
        r1 = unit.radius_small    # bottom radius
        r2 = unit.radius_large    # top radius
        h = unit.height           # height
        nθ, nz = 60, 30           # mesh resolution

        θ = range(0, 2π; length=nθ)
        z = range(0, h; length=nz)

        # Make 2D grids for Θ and Z:
        Θ = repeat(θ', nz, 1)
        Z = repeat(z, 1, nθ)

        # Radius varies linearly with z:
        R = r1 .+ (r2 - r1) / h .* Z

        # Parametric coords for the lateral surface:
        Xc = R .* cos.(Θ)
        Yc = R .* sin.(Θ)
        Zc = Z

        # Start with the side‐surface
        p = Plots.surface(Xc, Yc, Zc; alpha=0.5, color=:blue, legend=false)

        # Add bottom cap
        # Parameterize a filled disc at z = 0
        pr = range(0, r1; length=25)        # radial coordinate
        Θb = repeat(θ, 1, length(pr))       #  nθ × nρ
        Rb = repeat(pr', nθ, 1)             #  nθ × nρ

        Xb = Rb .* cos.(Θb)
        Yb = Rb .* sin.(Θb)
        Zb = zeros(size(Xb))                # all at z = 0

        Plots.surface!(p, Xb, Yb, Zb; color=:blue, alpha=0.65, legend=false)

        # add black contour lines
        Plots.plot3d!(p,
                      r1 .* cos.(θ),
                      r1 .* sin.(θ),
                      zeros(length(θ));
                      color=:black, lw=1)
        # top rim
        Plots.plot3d!(p,
                      r2 .* cos.(θ),
                      r2 .* sin.(θ),
                      fill(h, length(θ));
                      color=:black, lw=1)
    end

    fig_size = 1.2 * max(unit.radius_large, unit.height)
    Plots.plot!(p;
                xlabel="x-coordinate [m]", ylabel="y-coordinate [m]", zlabel="height [m]",
                legend=false,
                # proj_type=:ortho, # is not working...
                aspect_ratio=:equal,
                xlims=(-fig_size, fig_size),
                ylims=(-fig_size, fig_size),
                zlims=(-fig_size, fig_size),
                size=(1000, 1000),
                camera=(45, 20),
                title="3D shape of the STES $(unit.uac) (cross-section: $(unit.shape))")
    fig_name = "STES_3D_shape_$(unit.uac)"
    Plots.savefig(p, output_path * "/" * fig_name * ".html")

    return true
end

function plot_optional_figures_end(unit::SeasonalThermalStorage, sim_params::Dict{String,Any},
                                   output_path::String)::Bool
    # ====================================
    # Plot temperature distribution over time    
    # ====================================
    x_vals_datetime = [add_ignoring_leap_days(sim_params["start_date_output"],
                                              Dates.Second((s - 1) * sim_params["time_step_seconds"]))
                       for s in 1:sim_params["number_of_time_steps_output"]]
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
                                 Dates.value(Year(sim_params["start_date_output"])):Dates.value(Year(sim_params["end_date"]))
                             if isleapyear(year)])

    layout = PlotlyJS.Layout(; title_text="Temperature distribution over time in STES $(unit.uac)",
                             xaxis_title_text="Date",
                             yaxis_title_text="Temperature [°C]",
                             xaxis=PlotlyJS.attr(; type="date",
                                                 rangebreaks=[Dict("values" => leap_days_str)]))
    p = PlotlyJS.plot(traces, layout)

    fig_name = "temperature_distribution_STES_$(unit.uac).html"
    PlotlyJS.savefig(p, output_path * "/" * fig_name)

    # ====================================
    # Plot temperature differences to surrounding over time    
    # ====================================
    x_vals_datetime = [add_ignoring_leap_days(sim_params["start_date_output"],
                                              Dates.Second((s - 1) * sim_params["time_step_seconds"]))
                       for s in 1:sim_params["number_of_time_steps_output"]]
    layers_to_plot = (unit.number_of_layer_total + 2):-1:1
    traces = PlotlyJS.GenericTrace[]
    for i in layers_to_plot
        plot_label = "STES Layer $(i-1)"
        if i == 1
            plot_label = "Below storage"
        elseif i == 2
            plot_label = "Bottom STES layer"
        elseif i == (unit.number_of_layer_total + 2) - 1
            plot_label = "Top STES layer"
        elseif i == (unit.number_of_layer_total + 2)
            plot_label = "Above storage"
        end
        trace = PlotlyJS.scatter(; x=x_vals_datetime, y=unit.temp_difference_to_surrounding_output[:, i], mode="lines",
                                 name=plot_label)
        push!(traces, trace)
    end

    leap_days_str = string.([Date(year, 2, 29)
                             for year in
                                 Dates.value(Year(sim_params["start_date_output"])):Dates.value(Year(sim_params["end_date"]))
                             if isleapyear(year)])

    layout = PlotlyJS.Layout(;
                             title_text="Temperature difference to surrounding over time in STES $(unit.uac) (positive: storage is colder than surrounding)",
                             xaxis_title_text="Date",
                             yaxis_title_text="Temperature difference [K]",
                             xaxis=PlotlyJS.attr(; type="date",
                                                 rangebreaks=[Dict("values" => leap_days_str)]))
    p = PlotlyJS.plot(traces, layout)

    fig_name = "temperature_difference_surrounding_STES_$(unit.uac).html"
    PlotlyJS.savefig(p, output_path * "/" * fig_name)

    if unit.ground_model == "FEM"
        # ====================================
        # plot soil + storage temperatures
        # ====================================
        # possibility to plot the whole storage and all soil, not only half of the domain. May be useful in the future...
        mirror_domain = false

        # Dimensions: (time, z, r_half)
        nt, nz, nr = size(unit.soil_temperature_field_output)

        # Axes (cell centers) for simulated half-domain (r ≥ 0)
        r_half = [unit.soil_dr[1] / 2; cumsum(unit.soil_dr_mesh)]
        z_abs = [unit.soil_dz[1] / 2; cumsum(unit.soil_dz_mesh)]

        # Build radius axis (mirrored or not)
        r_abs = mirror_domain ? vcat(-reverse(r_half), r_half) : r_half

        # Temperature range with small padding
        Tmin_raw = minimum(unit.soil_temperature_field_output)
        Tmax_raw = maximum(unit.soil_temperature_field_output)
        Tmin = (Tmin_raw < 0.0) ? 1.1 * Tmin_raw : 0.9 * Tmin_raw
        Tmax = (Tmax_raw < 0.0) ? 0.9 * Tmax_raw : 1.1 * Tmax_raw

        # Spatial ranges
        zmin, zmax = extrema(z_abs)
        rmin, rmax = extrema(r_abs)

        # Enforce equal spatial scaling (z & r); temperature is independent
        Δz = max(zmax - zmin, eps(Float64))
        Δr = max(rmax - rmin, eps(Float64))
        spatial_ref = max(Δz, Δr)
        aspect_x = Δz / spatial_ref      # depth axis
        aspect_y = Δr / spatial_ref      # radius axis
        aspect_z = 1.0                   # temperature axis (no special scaling)

        # Helper: temperature slice for given time index (handles mirroring)
        temp_slice(t) = begin
            A = unit.soil_temperature_field_output[t, :, :]  # (nz × nr) half-domain
            mirror_domain ? hcat(reverse(A; dims=2), A) : A
        end

        # ================= GLMakie live viewer =================
        try
            GLMakie.activate!()

            time_idx = Observable(1)
            surf_obs = @lift(temp_slice($time_idx))

            f_gl = Figure(; size=(1000, 800))
            ax_gl = Axis3(f_gl[1, 1];
                          xlabel="Depth z [m]",
                          ylabel="Radius r [m]",
                          zlabel="Temperature [°C]")

            ax_gl.limits = ((zmin, zmax), (rmin, rmax), (Tmin, Tmax))
            ax_gl.aspect = (aspect_x, aspect_y, aspect_z)

            # surface!(ax_gl, z_abs, r_abs, surf_obs)
            scatter!(ax_gl, z_abs, r_abs, surf_obs; markersize=3.5)

            # Simple slider instead of SliderGrid (more robust across Makie versions)
            slider_gl = GLMakie.Slider(f_gl[2, 1]; range=1:nt, startvalue=1, horizontal=true)

            on(slider_gl.value) do v
                time_idx[] = v
            end

            display(f_gl)
        catch e
            @warn "Error while building GLMakie soil/STES figure for $(unit.uac): $(e)"
        end
    end

    return true
end

function control(unit::SeasonalThermalStorage,
                 components::Grouping,
                 sim_params::Dict{String,Any})
    update(unit.controller)

    # write old temperature field for output
    if sim_params["current_date"] >= sim_params["start_date_output"]
        time_idx = Int(sim_params["time_since_output"] / sim_params["time_step_seconds"]) + 1
        unit.temp_difference_to_surrounding_output[time_idx, 1] = unit.effective_ambient_temperature_bottom -
                                                                  unit.temperature_segments[1]
        unit.temp_difference_to_surrounding_output[time_idx, 2:(end - 1)] = unit.effective_ambient_temperature_barrels .-
                                                                            unit.temperature_segments
        unit.temp_difference_to_surrounding_output[time_idx, end] = unit.effective_ambient_temperature_top -
                                                                    unit.temperature_segments[end]

        unit.temp_distribution_output[time_idx, :] = copy(unit.temperature_segments)
    end

    # update ambient/ground boundary conditions for the current step
    if unit.ambient_temperature_profile !== nothing
        unit.ambient_temperature = Profiles.value_at_time(unit.ambient_temperature_profile, sim_params)
    end
    if unit.ground_temperature_profile !== nothing
        unit.ground_temperature = Profiles.value_at_time(unit.ground_temperature_profile, sim_params)
    end

    if unit.ground_model == "FEM"
        # calculate effective ambient temperature for each layer (unified FEM)
        update_ground_fem_unified_and_set_Teff!(unit, sim_params)

        # save unified FEM soil field for output
        if sim_params["current_date"] >= sim_params["start_date_output"]
            sidx = Int(sim_params["time_since_output"] / sim_params["time_step_seconds"]) + 1
            if size(unit.soil_temperature_field_output, 1) >= sidx
                unit.soil_temperature_field_output[sidx, :, :] = vis_field_with_tank(unit)
            end
        end
    elseif unit.ground_model == "simple"
        # update surrounding effective_ambient_temperature
        unit.effective_ambient_temperature_barrels = [i <= unit.number_of_STES_layer_below_ground ?
                                                      unit.ground_temperature : unit.ambient_temperature
                                                      for i in 1:(unit.number_of_layer_total)]
        unit.effective_ambient_temperature_top = unit.ambient_temperature
        unit.effective_ambient_temperature_bottom = unit.ground_temperature
    end

    # set current_energy_input_return_temperature (use bottom layer)
    unit.current_energy_input_return_temperature = unit.temperature_segments[1]

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
        max_input_energy[input_idx] = min(energy_supply[input_idx],
                                          calculate_max_input_energy_by_temperature(unit,
                                                                                    temperature,
                                                                                    unit.temperature_segments,
                                                                                    sim_params))
    end
    if isempty(temperatures_charging)
        max_input_energy = [0.0]
        temperatures_charging = [nothing]
        source_uac = [nothing]
    end
    set_max_energy!(unit.input_interfaces[unit.m_heat_in], max_input_energy, temperatures_charging,
                    temperatures_charging, source_uac, false, true)

    # calculate maximum energies for output
    if unit.max_unload_rate_mass !== nothing
        current_max_output_energy_from_mass = convert_mass_in_energy(unit.max_unload_rate_mass * unit.volume *
                                                                     unit.rho_medium *
                                                                     (sim_params["time_step_seconds"] / 60 / 60),
                                                                     unit.low_temperature,
                                                                     unit.current_max_output_temperature,
                                                                     unit.cp_medium)
        unit.current_max_output_energy = max(min(min(unit.load, unit.max_output_energy),
                                                 current_max_output_energy_from_mass), 0.0)
    else
        unit.current_max_output_energy = max(min(unit.load, unit.max_output_energy), 0.0)
    end
    set_max_energy!(unit.output_interfaces[unit.m_heat_out], unit.current_max_output_energy, nothing,
                    unit.current_max_output_temperature)
end
# This function sets the current_min_input_temperature and current_max_output_temperature 
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
            # If there is a inverse temperature distribution within the storage, search for the first layer
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

function current_max_input_energy_vol(unit::SeasonalThermalStorage, temp_in::Temperature,
                                      lower_temp::Temperature, sim_params::Dict{String,Any})::Float64
    if unit.max_load_rate_mass === nothing
        return Inf
    else
        # calculate maximum energy from mass flow / volume limit for current temperature
        return convert_mass_in_energy(unit.max_load_rate_mass * unit.volume * unit.rho_medium *
                                      (sim_params["time_step_seconds"] / 60 / 60),
                                      lower_temp,
                                      temp_in,
                                      unit.cp_medium)
    end
end

# This function will only be called from a Bus, never in a 1-to-1 connection to another component.
# energy_input and temperatures_input are all inputs that have already been put into the STES.
function recalculate_max_energy(unit::SeasonalThermalStorage,
                                energy_input::Union{Float64,Vector{Float64}},
                                temperatures_input::Union{Temperature,Vector{Temperature}},
                                max_energy::EnergySystems.MaxEnergy,
                                sim_params::Dict{String,Any})::EnergySystems.MaxEnergy
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
    # and limit it by the the already filled energy and the charging rate limit
    for (input_idx, temperature) in enumerate(highest.(max_energy.temperature_max, max_energy.temperature_min))
        max_energy.max_energy[input_idx] = min(calculate_max_input_energy_by_temperature(unit,
                                                                                         temperature,
                                                                                         temperature_segments_temporary,
                                                                                         sim_params),
                                               min(unit.max_input_energy - sum(energy_input; init=0.0),
                                                   max_energy.max_energy[input_idx]))
    end

    return max_energy
end

"""
get_input_temperature_bounds(unit::SeasonalThermalStorage)

Function used by the control module negotiate_temperature to get the current minimum 
and maximum temperatures for charging.

# Arguments
- `unit::SeasonalThermalStorage`: The seasonal thermal storage unit for which the calculation is performed.

# Returns
- `input_min_temperature::Temperature`: The current minimal temperature that can be taken in the input
- `input_max_temperature::Temperature`: The current maximum temperature that can be taken in the input

"""
function get_input_temperature_bounds(unit::SeasonalThermalStorage)::Tuple{Temperature,Temperature}
    input_min_temperature = unit.current_min_input_temperature
    input_max_temperature = unit.high_temperature
    return input_min_temperature, input_max_temperature
end

"""
calculate_max_input_energy_by_temperature(unit::SeasonalThermalStorage, actual_input_temp::Temperature,
  current_temperature_distribution::Vector{Temperature},
  sim_params::Dict{String,Any}) -> Float64

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
function calculate_max_input_energy_by_temperature(unit::SeasonalThermalStorage,
                                                   actual_input_temp::Temperature,
                                                   current_temperature_distribution::Vector{Temperature},
                                                   sim_params::Dict{String,Any})::Float64
    # reduce maximal energy if the input temperature is near the temperature of the lower layer to avoid
    # numerical problems and pulsing during loading. Otherwise, it can happen that the STES is not taking
    # the energy that it is supposed to take, or a unreasonable high temporal resolution has to be used to
    # calculate the new temperature distribution within the STES.
    red_factor = 0.0 + clamp(actual_input_temp - current_temperature_distribution[2], 0.0, 5.0) * 1.0 / 5.0
    input_energy_limit = min(unit.max_input_energy,
                             current_max_input_energy_vol(unit,
                                                          actual_input_temp,
                                                          current_temperature_distribution[1],
                                                          sim_params))
    return min(input_energy_limit,
               red_factor *
               sum(convert_mass_in_energy(volume * unit.rho_medium, current_temperature_distribution[layer],
                                          actual_input_temp, unit.cp_medium)
                   for (layer, volume) in enumerate(unit.volume_segments)
                   if current_temperature_distribution[layer] < actual_input_temp; init=0.0))
end

"""
calculate_input_energy_from_input_temperature(unit, actual_input_temp, sim_params)

Wrapper function to calculate the maximum input energy for a given input temperature for the STES.

Inputs:
- `unit::SeasonalThermalStorage`: The seasonal thermal storage unit for which the calculation is performed.
- `actual_input_temp::Temperature`: The temperature of the input energy being added to the storage.
- `sim_params::Dict{String,Any}`: A dictionary containing simulation parameters

Returns:
- `Float64`: The total maximum energy that can be added to the storage at the actual_input_temp

"""
function calculate_input_energy_from_input_temperature(unit::SeasonalThermalStorage,
                                                       actual_input_temp::Temperature,
                                                       sim_params::Dict{String,Any})::Float64
    return calculate_max_input_energy_by_temperature(unit,
                                                     actual_input_temp,
                                                     unit.temperature_segments,
                                                     sim_params)
end

""" 
convert_kJ_in_Wh(energy::Float64)

takes energy in [kJ] and convert it to [Wh]
"""
function convert_kJ_in_Wh(energy::Float64)::Float64
    return energy / 3.6
end

"""
convert_energy_in_mass(energy, temp_low, temp_high, cp, roh)

 calculates mass [kg] from energy [Wh] 

 # Arguments
- `energy::Float64`: energy to convert [Wh]
- `temp_low::Temperature`: lower temperature [°C]
- `temp_high::Temperature`: upper temperature [°C]
- `cp::Float64`: specific heat capacity [kJ/KgK]

"""
function convert_energy_in_mass(energy::Float64, temp_low::Temperature, temp_high::Temperature, cp::Float64)::Float64
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
- `temp_low::Temperature`: lower temperature [°C]
- `temp_high::Temperature`: upper temperature [°C]
- `cp::Float64`: specific heat capacity [kJ/KgK]

"""
function convert_mass_in_energy(mass::Float64, temp_low::Temperature, temp_high::Temperature, cp::Float64)::Float64
    return mass * convert_kJ_in_Wh(cp) * (temp_high - temp_low)
end

"""
update_STES(unit::SeasonalThermalStorage,
            energy_input::Vector{Float64},
            temperatures_input::Vector{Temperature},
            energy_output::Float64,
            temperature_output::Temperature,
            sim_params::Dict{String,Any};
            temporary_calculation::Bool=false)

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

    # mass_in and mass_out should be both positive! 
    # mass_in and the corresponding temperature can be a vector and will be fitted to the most suitable layer.
    # mass_out and the corresponding temperature are single values that will always be drawn from the top layer
    mass_in = convert_energy_in_mass.(energy_input, t_old[lower_node], temperatures_input,
                                      unit.cp_medium)
    mass_out = convert_energy_in_mass(energy_output, unit.low_temperature, temperature_output,
                                      unit.cp_medium)

    # Check if the mass flow is greater than the volume of the smallest segment. 
    # If yes, the internal time step is reduced to avoid numerical instabilities.
    unit.mass_out_sum = sum(mass_out; init=0.0)
    unit.mass_in_sum = sum(mass_in; init=0.0)
    mass_of_smallest_segment = minimum(unit.volume_segments) * unit.rho_medium
    factor = 5
    if max(unit.mass_out_sum, unit.mass_in_sum) > (mass_of_smallest_segment / factor) &&
       (unit.mass_out_sum > 0.0 || unit.mass_in_sum > 0.0)
        number_of_internal_timesteps = ceil(max(unit.mass_out_sum, unit.mass_in_sum) /
                                            (mass_of_smallest_segment / factor))
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
        mass_in_temp, mass_in_vec = calculate_mass_temperature_charging(unit, t_old,
                                                                        mass_in ./ number_of_internal_timesteps,
                                                                        lower_node, sim_params)
        # mass flow and corresponding temperatures into each layer during discharging
        mass_out_temp, mass_out_vec = calculate_mass_temperature_discharging(unit, t_old,
                                                                             mass_out ./
                                                                             number_of_internal_timesteps,
                                                                             return_temperature_input,
                                                                             lower_node, sim_params)
        for n in 1:(unit.number_of_layer_total)
            if n == 1  # bottom layer, single-side
                t_new[n] = t_old[n] +
                           consider_losses *
                           (3600 * unit.diffusion_coefficient * (t_old[n + 1] - t_old[n]) / unit.dz_normalized[n]^2 +    # thermal diffusion
                            unit.sigma[n] * (unit.effective_ambient_temperature_bottom - t_old[n])) * dt +               # losses through bottom and side walls
                           # unit.lambda[n] * (Q_in_out)[n] +                                                            # thermal input and output
                           unit.phi[n] * mass_in_vec[n] * (mass_in_temp[n] - t_old[n]) +                                 # mass input
                           unit.phi[n] * mass_out_vec[n] * (mass_out_temp[n] - t_old[n])                                 # mass output
            elseif n == unit.number_of_layer_total  # top layer, single-side
                t_new[n] = t_old[n] +
                           consider_losses *
                           (3600 * unit.diffusion_coefficient * (t_old[n - 1] - t_old[n]) / unit.dz_normalized[n]^2 +    # thermal diffusion
                            unit.sigma[n] * (unit.effective_ambient_temperature_top - t_old[n])) * dt +                  # losses through lid and side walls
                           # unit.lambda[n] * Q_in_out[n] +                                                              # thermal input and output
                           unit.phi[n] * mass_in_vec[n] * (mass_in_temp[n] - t_old[n]) +                                 # mass input
                           unit.phi[n] * mass_out_vec[n] * (mass_out_temp[n] - t_old[n])                                 # mass output
            else       # mid layer
                t_new[n] = t_old[n] +
                           consider_losses *
                           (3600 * unit.diffusion_coefficient * (t_old[n + 1] + t_old[n - 1] - 2 * t_old[n]) /
                            unit.dz_normalized[n]^2 +                                                                    # thermal diffusion
                            unit.sigma[n] * (unit.effective_ambient_temperature_barrels[n] - t_old[n])) * dt +            # losses through side walls
                           # unit.lambda[n] * Q_in_out[n] +                                                              # thermal input and output
                           unit.phi[n] * mass_in_vec[n] * (mass_in_temp[n] - t_old[n]) +                                 # mass input
                           unit.phi[n] * mass_out_vec[n] * (mass_out_temp[n] - t_old[n])                                 # mass output
            end

            if n > 1   # mixing due to buoyancy effects, if temperature gradient is present
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
calculate_mass_temperature_charging(unit::SeasonalThermalStorage, t_old::Vector{Temperature},
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
function calculate_mass_temperature_charging(unit::SeasonalThermalStorage, t_old::Vector{Temperature},
                                             mass_in::Vector{Float64}, lower_node::Int,
                                             sim_params::Dict{String,Any})::Tuple{Vector{Float64},Vector{Float64}}
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

        # calculate temperature of mass flow into the layer due to discharging
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
calculate_mass_temperature_discharging(unit::SeasonalThermalStorage, t_old::Vector{Temperature},
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
function calculate_mass_temperature_discharging(unit::SeasonalThermalStorage, t_old::Vector{Temperature},
                                                mass_out::Float64, return_temperature_input::Temperature,
                                                lower_node::Int,
                                                sim_params::Dict{String,Any})::Tuple{Vector{Float64},Vector{Float64}}
    upper_node_discharging = unit.number_of_layer_total - unit.output_layer_from_top + 1
    mass_out_temp = zeros(unit.number_of_layer_total)
    mass_in_layer = [mass_out]
    temperature_in_layer = [return_temperature_input]
    for layer in lower_node:upper_node_discharging
        # calculate temperature of mass flow into the layer due to discharging
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
        handle_component_update!(unit, "process", sim_params)
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

    handle_component_update!(unit, "process", sim_params)
end

function handle_component_update!(unit::SeasonalThermalStorage, step::String, sim_params::Dict{String,Any})
    if step == "process"
        unit.process_done = true
    elseif step == "load"
        unit.load_done = true
    end
    if unit.process_done && unit.load_done
        # update component
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
        # reset
        unit.process_done = false
        unit.load_done = false
    end
end

function load(unit::SeasonalThermalStorage, sim_params::Dict{String,Any})
    inface = unit.input_interfaces[unit.m_heat_in]
    exchanges = balance_on(inface, inface.source)
    energy_available = balance(exchanges) + energy_potential(exchanges)

    # shortcut if there is no energy to be used
    if energy_available <= sim_params["epsilon"]
        handle_component_update!(unit, "load", sim_params)
        set_max_energy!(unit.input_interfaces[unit.m_heat_in], 0.0)
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

        # all energies in the exchange can be used as it was already made sure than they do not exceed the STES limit
        sub!(inface, exchange_energy_available, current_exchange_temperature, nothing, exchange.purpose_uac)
        push!(unit.current_energy_input, exchange_energy_available)
        push!(unit.current_temperature_input, current_exchange_temperature)
    end

    handle_component_update!(unit, "load", sim_params)
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
    return [string(unit.m_heat_in) * ":IN",
            string(unit.m_heat_out) * ":OUT",
            "Load",
            "Load%",
            "Capacity",
            "LossesGains",
            "CurrentMaxOutTemp",
            "GroundTemperature",
            "MassInput",
            "MassOutput",
            "Temperature_upper",
            "Temperature_three_quarter",
            "Temperature_middle",
            "Temperature_one_quarter",
            "Temperature_lower"]
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
    elseif key.value_key == "LossesGains"
        return -unit.losses
    elseif key.value_key == "CurrentMaxOutTemp"
        return unit.current_max_output_temperature
    elseif key.value_key == "GroundTemperature"
        return unit.ground_temperature
    elseif key.value_key == "MassInput"
        return unit.mass_in_sum
    elseif key.value_key == "MassOutput"
        return unit.mass_out_sum
    elseif key.value_key == "Temperature_lower"
        return unit.temperature_segments[1]
    elseif key.value_key == "Temperature_one_quarter"
        return unit.temperature_segments[Int(round(0.25 * unit.number_of_layer_total))]
    elseif key.value_key == "Temperature_middle"
        return unit.temperature_segments[Int(round(0.5 * unit.number_of_layer_total))]
    elseif key.value_key == "Temperature_three_quarter"
        return unit.temperature_segments[Int(round(0.75 * unit.number_of_layer_total))]
    elseif key.value_key == "Temperature_upper"
        return unit.temperature_segments[end]
    end
    throw(KeyError(key.value_key))
end

# === Ground coupling FEM for STES (unified axisymmetric r–z domain) ===

"""
    create_geometric_mesh(min_mesh_width::Float64,
                          max_mesh_width::Float64,
                          expansion_factor::Float64,
                          bound::Float64) -> Vector{Float64}

Build a 1D vector of cell widths `dx` that covers length `bound` using geometric growth
from the origin. Start at `min_mesh_width`, multiply by `expansion_factor` (capped at
`max_mesh_width`) until near `bound`, then trim/append the last cell so `sum(dx) == bound`.
Returns `Float64[]` if `bound ≤ 0`.

Args:
- `min_mesh_width` (>0) start size.
- `max_mesh_width` (≥ min) cap.
- `expansion_factor` (>1 for growth; ≈1 → near-uniform).
- `bound` total length.
"""
function create_geometric_mesh(min_mesh_width::Float64,
                               max_mesh_width::Float64,
                               expansion_factor::Float64,
                               bound::Float64)
    dx = Float64[]
    if bound <= 0.0
        return dx
    end
    push!(dx, min_mesh_width)
    while sum(dx) + dx[end] < bound
        push!(dx, min(dx[end] * expansion_factor, max_mesh_width))
    end
    # trim last cell to hit bound exactly
    excess = sum(dx) - bound
    if excess > 0
        dx[end] -= excess
    elseif sum(dx) < bound
        push!(dx, bound - sum(dx))
    end
    return dx
end

"""
    create_geometric_mesh_two_sided(min_mesh_width::Float64,
                                    max_mesh_width::Float64,
                                    expansion_factor::Float64,
                                    bound::Float64)

Build a 1D mesh on [0, bound] that is fine near both ends:
- geometric growth from each end starting at `min_mesh_width` with factor `expansion_factor`,
  capped at `max_mesh_width`
- middle region filled with equidistant cells of width `max_mesh_width`
- robust for small `bound` (returns two equal cells if the domain is too short)

Returns a vector `dx` with sum(dx) == bound (up to floating roundoff).
"""
function create_geometric_mesh_two_sided(min_mesh_width::Float64,
                                         max_mesh_width::Float64,
                                         expansion_factor::Float64,
                                         bound::Float64)
    if bound <= 0.0
        return Float64[]
    end

    # very short domain: split into two equal cells
    if bound <= 2 * min_mesh_width
        return [bound / 2, bound / 2]
    end

    # geometric ramp sequence up to the cap
    w = Float64[min_mesh_width]
    while w[end] < max_mesh_width - eps(Float64)
        nextw = min(w[end] * expansion_factor, max_mesh_width)
        if nextw == w[end]      # handles expansion_factor == 1
            break
        end
        push!(w, nextw)
    end

    # choose how many ramp cells per side can fit
    s = 0.0
    n = 0
    for j in eachindex(w)
        if 2 * (s + w[j]) <= bound
            s += w[j]
            n = j
        else
            break
        end
    end

    # middle length and fill with max-sized equidistant cells
    middle_len = bound - 2 * s
    dx = Float64[]
    append!(dx, w[1:n])

    if middle_len > 0
        if middle_len <= max_mesh_width + 1e-12
            push!(dx, middle_len)
        else
            nmid = floor(Int, middle_len / max_mesh_width)
            append!(dx, fill(max_mesh_width, nmid))
            rem = middle_len - nmid * max_mesh_width
            if rem > 1e-12
                push!(dx, rem)
            end
        end
    end

    append!(dx, reverse(w[1:n]))

    # tiny correction to hit 'bound' exactly
    err = bound - sum(dx)
    if abs(err) > 1e-10
        dx[end] += err
    end
    return dx
end

# map depth z [m] (measured from ground surface downward) to soil layer properties
# Note that currently the soil properties are considered at a single depth without taking overlaps into account!
function soil_props_at_depth(unit::SeasonalThermalStorage, z::Float64)
    zs = unit.ground_layers_depths

    # ensure last depth >= domain depth
    if zs[end] < unit.ground_domain_depth
        unit.ground_layers_depths[end] = unit.ground_domain_depth
        zs = unit.ground_layers_depths
        @info "In STES $(unit.uac), the given `ground_layers_depths` has been extended to $(unit.ground_domain_depth) m " *
              "to cover the whole ground depth. The corresponding ground properties were copied from the last given depth."
    end
    # find interval
    j = 1
    for idx in 1:(length(zs) - 1)
        if zs[idx] <= z < zs[idx + 1]
            j = idx
            break
        end
    end
    k = unit.ground_layers_k[min(j, length(unit.ground_layers_k))]
    rho = unit.ground_layers_rho[min(j, length(unit.ground_layers_rho))]
    cp = unit.ground_layers_cp[min(j, length(unit.ground_layers_cp))]
    return k, rho, cp
end

# Return radii to be used for soil/axisymmetric coupling.
# - For round shapes: use actual radii.
# - For quadratic shapes: use area-equivalent circular radii
#   so that π R^2 matches the actual lid/bottom areas.
function equiv_radii_for_ground(unit::SeasonalThermalStorage)
    if unit.shape == "quadratic"
        # Protect against uninitialized/zero
        R_bot = unit.surface_area_bottom > 0 ? sqrt(unit.surface_area_bottom / pi) : unit.radius_small
        R_top = unit.surface_area_lid > 0 ? sqrt(unit.surface_area_lid / pi) : unit.radius_large
        return R_bot, R_top
    else
        return unit.radius_small, unit.radius_large
    end
end

function prepare_ground_fem_unified!(unit::SeasonalThermalStorage)
    # ----------------------------
    # 1) Basic geometry
    # ----------------------------
    R_bot_eq, R_top_eq = equiv_radii_for_ground(unit)

    R_dom = unit.ground_domain_radius
    R_small = max(R_bot_eq, 0.0)
    R_small = min(R_small, R_dom)  # cannot exceed domain

    # radius of side wall where it meets ground surface (top of buried section)
    R_ground = if unit.height > 0
        R_bot_eq + (R_top_eq - R_bot_eq) * (unit.h_stes_buried / unit.height)
    else
        R_bot_eq
    end

    # effective wall radius in contact with soil
    R_wall = clamp(R_ground,
                   min(R_bot_eq, R_top_eq),
                   max(R_bot_eq, R_top_eq))
    R_wall = min(R_wall, R_dom)

    has_slope = R_wall > R_small + 1e-9

    # ----------------------------
    # 2) Mesh resolution knobs
    # ----------------------------
    mindz = minimum(unit.dz)

    min_w, max_w, n_wall, ef = unit.ground_accuracy_mode == "very_rough" ? (mindz / 1, 16.0, 1, 2.0) :
                               unit.ground_accuracy_mode == "rough" ? (mindz / 2, 8.0, 1, 2.0) :
                               unit.ground_accuracy_mode == "normal" ? (mindz / 3, 4.0, 1, 2.0) :
                               unit.ground_accuracy_mode == "high" ? (mindz / 4, 2.0, 2, 2.0) :
                               unit.ground_accuracy_mode == "very_high" ? (mindz / 8, 1.0, 3, 2.0) :
                               error("In STES $(unit.uac), ground_accuracy_mode must be one of: " *
                                     "very_rough, rough, normal, high, very_high.")

    min_w = max(min_w, 1e-4)
    max_w = max(max_w, min_w)

    # ----------------------------
    # 3) z mesh: exact below STES, then tail
    # ----------------------------
    dz = unit.number_of_STES_layer_below_ground > 0 ? copy(unit.dz[1:(unit.number_of_STES_layer_below_ground)]) :
         Float64[]
    remaining = max(unit.ground_domain_depth - unit.h_stes_buried, min_w)

    if remaining > unit.epsilon_geometry
        tail = create_geometric_mesh(min_w, max_w, ef, remaining)
        dz = vcat(dz, tail)
    end
    if isempty(dz)
        dz = [min_w]
    end

    unit.soil_dz = dz
    unit.soil_dz_mesh = [(dz[i] + dz[i + 1]) / 2 for i in 1:(length(dz) - 1)]

    zc = [sum(dz[1:(h - 1)]) + dz[h] / 2 for h in eachindex(dz)]
    unit.soil_z_centers = zc

    # ----------------------------
    # 4) r mesh (global):
    #    [0, R_small]   : coarser center
    #    [R_small, R_wall]: refined along buried sidewall
    #    [R_wall, R_dom]: geometric coarsening
    # ----------------------------
    dr_segments = Float64[]

    # --- Region 1: 0 .. R_small ---
    if R_small > 0
        append!(dr_segments,
                create_geometric_mesh_two_sided(min_w, max_w, ef, R_small))
    end

    # --- Region 2: R_small .. R_wall (only if sloped and buried wall exists) ---
    if has_slope && (R_wall > R_small + unit.epsilon_geometry)
        wall_width = R_wall - R_small
        # aim for ~ number_of_STES_layer_below_ground * n_wall cells across wall-contact band
        target_cells = max(unit.number_of_STES_layer_below_ground * n_wall, 1)
        dr_wall = wall_width / target_cells
        dr_wall = max(dr_wall, min_w / 4)  # reasonably fine
        n2 = max(1, round(Int, wall_width / dr_wall))
        dr2 = wall_width / n2
        append!(dr_segments, fill(dr2, n2))
    end

    # --- Region 3: R_wall .. R_dom ---
    used_r = sum(dr_segments)
    remaining_r = max(R_dom - used_r, 0.0)
    if remaining_r > unit.epsilon_geometry
        start_w = max(last(dr_segments), min_w)
        append!(dr_segments,
                create_geometric_mesh(start_w, max_w, ef, remaining_r))
    end

    if isempty(dr_segments)
        dr_segments = [min_w]
    end

    dr = dr_segments
    unit.soil_dr = dr
    unit.soil_dr_mesh = [(dr[i] + dr[i + 1]) / 2 for i in 1:(length(dr) - 1)]
    rc = [sum(dr[1:(i - 1)]) + dr[i] / 2 for i in eachindex(dr)]
    unit.soil_r_centers = rc

    # ----------------------------
    # 5) Initial soil temperature
    # ----------------------------
    T_top = unit.ambient_temperature
    T_bot = unit.ground_temperature
    total_depth = sum(dz)

    T_rows = [T_top + (T_bot - T_top) * (zc[h] / total_depth)
              for h in 1:length(dz)]

    nz = length(dz)
    nr = length(dr)
    unit.soil_t1 = [T_rows[h] for h in 1:nz, _ in 1:nr]
    unit.soil_t2 = similar(unit.soil_t1)

    return nothing
end

# Unified implicit Euler solve and extraction of wall/base interface temperatures
function solve_soil_unified!(unit::SeasonalThermalStorage, sim_params::Dict{String,Any})
    # aliases for convenience
    nz = length(unit.soil_dz)
    nr = length(unit.soil_dr)
    N = nz * nr

    dz = unit.soil_dz
    dr = unit.soil_dr
    dzm = unit.soil_dz_mesh
    drm = unit.soil_dr_mesh
    rc = unit.soil_r_centers

    I = Int[]
    J = Int[]
    V = Float64[]
    b = zeros(Float64, N)
    def_idx(h, i) = (h - 1) * nr + i

    Told = unit.soil_t1

    cum_dz_below = cumsum(unit.dz[1:(unit.number_of_STES_layer_below_ground)])  # bottom→top cumulative heights of buried STES layers

    for h in 1:nz, i in 1:nr
        p = def_idx(h, i)
        if !unit.cells_active[h, i]
            push!(I, p)
            push!(J, p)
            push!(V, 1.0)
            b[p] = Told[h, i]
            continue
        end

        # Deep boundary: Dirichlet to ground temp
        if h == nz
            push!(I, p)
            push!(J, p)
            push!(V, 1.0)
            b[p] = unit.ground_temperature
            continue
        end

        # capacity
        Vcell = 2pi * rc[i] * dr[i] * dz[h]
        C = unit.row_rho[h] * unit.row_cp[h] * Vcell
        aP = C / sim_params["time_step_seconds"]
        rhs = aP * Told[h, i]

        # EAST (i+1)
        if i < nr && unit.cells_active[h, i + 1]
            kE = (unit.row_k[h] + unit.row_k[h]) / 2
            r_e = rc[i] + dr[i] / 2
            A_e = 2pi * r_e * dz[h]
            d_e = drm[i]
            aE = kE * A_e / d_e
            aP += aE
            push!(I, p)
            push!(J, def_idx(h, i + 1))
            push!(V, -aE)
        end

        # WEST (i-1) or wall boundary
        if i > 1
            if unit.cells_active[h, i - 1]
                kW = (unit.row_k[h] + unit.row_k[h]) / 2
                r_w = rc[i] - dr[i] / 2
                A_w = 2pi * r_w * dz[h]
                d_w = drm[i - 1]
                aW = kW * A_w / d_w
                aP += aW
                push!(I, p)
                push!(J, def_idx(h, i - 1))
                push!(V, -aW)
            else
                # masked neighbor → wall Robin (below ground only)
                if h <= unit.number_of_STES_layer_below_ground
                    # map soil row h (depth from surface) to buried STES layer k (bottom→top)
                    z_soil = unit.soil_z_centers[h]
                    z_storage = unit.h_stes_buried - z_soil
                    k = searchsortedfirst(cum_dz_below, z_storage + eps(Float64))
                    k = clamp(k, 1, unit.number_of_STES_layer_below_ground)

                    # local wall face area represented by this cell
                    A_face = max(2pi * unit.radius_at_row[h] * dz[h], eps(Float64))

                    U = unit.thermal_transmission_barrels[k]
                    aP += U * A_face
                    rhs += U * A_face * unit.temperature_segments[k]
                end
            end
        end
        # (axis i==1 → zero flux; nothing to add)

        # SOUTH (h+1)
        if h < nz && unit.cells_active[h + 1, i]
            kS = (unit.row_k[h] + unit.row_k[h + 1]) / 2
            A_s = 2pi * rc[i] * dr[i]
            d_s = dzm[h]
            aS = kS * A_s / d_s
            aP += aS
            push!(I, p)
            push!(J, def_idx(h + 1, i))
            push!(V, -aS)
        end

        # NORTH (h-1) or base/top boundary
        if h == 1
            if unit.number_of_STES_layer_below_ground == 0 && rc[i] <= unit.equivalent_radius_from_bottom + 1e-12
                # handle case if STES has zero layer below ground
                A_n = 2pi * rc[i] * dr[i]
                hbot = unit.thermal_transmission_bottom
                aP += hbot * A_n
                rhs += hbot * A_n * unit.temperature_segments[1]
            else
                # surface convection
                A_n = 2pi * rc[i] * dr[i]
                hconv = unit.soil_surface_hconv
                aP += hconv * A_n
                rhs += hconv * A_n * unit.ambient_temperature
            end
        else
            if unit.cells_active[h - 1, i]
                kN = (unit.row_k[h] + unit.row_k[h - 1]) / 2
                A_n = 2pi * rc[i] * dr[i]
                d_n = dzm[h - 1]
                aN = kN * A_n / d_n
                aP += aN
                push!(I, p)
                push!(J, def_idx(h - 1, i))
                push!(V, -aN)
            else
                # masked neighbor → base Robin within equivalent bottom radius
                if (rc[i] <= unit.equivalent_radius_from_bottom + 1e-12) &&
                   (h == unit.number_of_STES_layer_below_ground + 1)
                    A_n = 2pi * rc[i] * dr[i]
                    hbot = unit.thermal_transmission_bottom
                    aP += hbot * A_n
                    rhs += hbot * A_n * unit.temperature_segments[1]
                end
            end
        end

        # finalize
        push!(I, p)
        push!(J, p)
        push!(V, aP)
        b[p] += rhs
    end

    A = sparse(I, J, V, N, N)
    Tvec = A \ b

    # map back
    Tnew = similar(Told)
    for h in 1:nz, i in 1:nr
        Tnew[h, i] = unit.cells_active[h, i] ? Tvec[def_idx(h, i)] : Told[h, i]
    end
    unit.soil_t2 .= Tnew
    unit.soil_t1 .= unit.soil_t2

    # interface temperatures
    T_wall_side = Float64[]
    for k in 1:(unit.number_of_STES_layer_below_ground)
        # center of STES layer k: measured from bottom upwards
        z_k = cum_dz_below[k] - unit.dz[k] / 2
        # corresponding soil depth below surface
        z_soil = unit.h_stes_buried - z_k

        # find soil row whose center is just below/at this depth
        h = searchsortedfirst(unit.soil_z_centers, z_soil + eps(Float64))
        h = clamp(h, 1, nz)

        # pick soil temp at first radius >= wall radius
        idx = findfirst(x -> x >= unit.radius_at_row[k], unit.soil_r_centers)
        val = idx === nothing ? unit.ground_temperature : unit.soil_t2[h, idx]

        push!(T_wall_side, val)
    end

    hbot = unit.number_of_STES_layer_below_ground + 1
    num = 0.0
    den = 0.0
    for i in 1:length(unit.soil_r_centers)
        if unit.soil_r_centers[i] <= unit.equivalent_radius_from_bottom + 1e-12
            w = 2pi * unit.soil_r_centers[i] * unit.soil_dr[i]
            num += unit.soil_t2[hbot, i] * w
            den += w
        else
            break
        end
    end
    T_base = num / max(den, eps(Float64))

    return T_wall_side, T_base
end

function update_ground_fem_unified_and_set_Teff!(unit::SeasonalThermalStorage, sim_params::Dict{String,Any})
    T_wall_side, T_base = solve_soil_unified!(unit, sim_params)

    # side below ground: each buried layer k gets its corresponding wall soil temperature
    # T_wall_side is already in order of the STES layer (bottom STES layer is index 1)
    for k in 1:(unit.number_of_STES_layer_below_ground)
        unit.effective_ambient_temperature_barrels[k] = T_wall_side[k]
    end
    # side above ground: ambient temperature
    for k in (unit.number_of_STES_layer_below_ground + 1):(unit.number_of_layer_total)
        unit.effective_ambient_temperature_barrels[k] = unit.ambient_temperature
    end

    # lid is assumed to always face the ambient
    unit.effective_ambient_temperature_top = unit.ambient_temperature
    # bottom is assumed to always face the ground
    unit.effective_ambient_temperature_bottom = T_base
end

# Build a view field where tank interior has the current layer temps
function vis_field_with_tank(unit::SeasonalThermalStorage)
    Tvis = copy(unit.soil_t1)
    nz, nr = length(unit.soil_dz), length(unit.soil_dr)
    cumz = cumsum(unit.dz)  # bottom→top cumulative layer heights

    for h in 1:nz
        z_soil = unit.soil_z_centers[h]
        if z_soil <= unit.h_stes_buried
            z_storage = unit.h_stes_buried - z_soil
            k = searchsortedfirst(cumz, z_storage + eps(Float64))  # layer index
            k = clamp(k, 1, unit.number_of_layer_total)
            for i in 1:nr
                if !unit.cells_active[h, i]
                    Tvis[h, i] = Float64(unit.temperature_segments[k])
                end
            end
        end
    end
    return Tvis
end

"""
    radius_at_row(unit::SeasonalThermalStorage, h::Int) -> Float64

Compute the STES (equivalent) wall radius at a given soil row center (depth-dependent for conical/frustum walls)

Args:
- `unit::SeasonalThermalStorage`: Uses these fields:
- `h::Int`: 1-based index into `unit.soil_z_centers`; must satisfy `1 ≤ h ≤ length(unit.soil_z_centers)`.

Returns:
- `Float64`: Wall radius at that depth [m]; returns `0.0` if the row lies below the tank base.
"""
function radius_at_row(unit::SeasonalThermalStorage, h::Int)::Float64
    z_soil = unit.soil_z_centers[h]

    if unit.number_of_STES_layer_below_ground <= 0
        return 0.0
    end

    # outside buried section → no side wall
    if z_soil > unit.h_stes_buried + unit.epsilon_geometry
        return 0.0
    end

    # Map to storage coordinate
    z_storage = unit.h_stes_buried - z_soil
    z_storage = clamp(z_storage, 0.0, unit.height)

    # Use equivalent radii for ground coupling
    R_bot_eq, R_top_eq = equiv_radii_for_ground(unit)

    return R_bot_eq +
           (R_top_eq - R_bot_eq) * (z_storage / max(unit.height, eps(Float64)))
end

"""
    cell_active(unit::SeasonalThermalStorage, h::Int, i::Int)::Bool

Active-soil mask for cell (h,i). Returns `true` if the cell is soil (part of the PDE),
`false` if it lies inside the STES.

# Arguments
- `unit::SeasonalThermalStorage`: The STES unit
- `h::Int`: Depth index (row).
- `i::Int`: Radius index (column).

# Returns
- `Bool`: `true` for soil/active; `false` for interior/masked.
"""
function cell_active(unit::SeasonalThermalStorage, h::Int, i::Int)::Bool
    # true → this cell is *soil* and part of the PDE
    # false → this cell lies inside the STES volume (masked, handled separately)

    # below the buried part  → always soil
    if unit.soil_z_centers[h] > unit.h_stes_buried + unit.epsilon_geometry
        return true
    end

    if unit.radius_at_row[h] <= unit.epsilon_geometry
        # no wall at this depth → soil
        return true
    end

    # inside tank if clearly left of wall; use tolerance to avoid isolated flips
    inside = (unit.soil_r_centers[i] < unit.radius_at_row[h] + unit.epsilon_geometry)

    return !inside
end

export SeasonalThermalStorage
