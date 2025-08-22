module Weatherdata

using Resie.Profiles
using Dates
using Proj
using Resie.SolarIrradiance

export WeatherData, gather_weather_data, get_weather_data_keys

"""
Custom error handler for exception "InputError".
Call with throw(InputError)
"""
struct InputError <: Exception end

"""
"""
mutable struct WeatherData
    """Ambient air temperature, in °C."""
    temp_ambient_air::Profile

    """Wind speed, in m/s."""
    wind_speed::Profile

    """Beam horizontal irradiation, in Wh/m^2."""
    beamHorIrr::Profile

    """Diffuse horizontal irradiation, in Wh/m^2."""
    difHorIrr::Profile

    """Global horizontal irradiation, in Wh/m^2."""
    globHorIrr::Profile

    """Horizontal long wave (infrared) irradiation, in Wh/m^2."""
    longWaveIrr::Profile

    """sunrise, in decimal hours."""
    sunrise::Profile

    """sunset, in decimal hours."""
    sunset::Profile

    """
    get_weather_data(weather_file_path, sim_params)
    
    Function to read in a weather file and hold the data. The data can either be
    a .dat file from the DWD (German weather service) or an EPW file (EnergyPlusWeather).
    
    The returned values are of type WeaterData containing profiles of type Profile.

    **Definition of time:**
    The weather data returned from this function corresponds to the timestep following
    the time indicated. While the ambient temperature is an instantaneous value from half the 
    timestep ahead of the current timestamp, the solar radiation data is the mean/sum of the 
    timestep ahead of the current timestamp. Data from EWP and DWD-dat files are converted
    accordingly.

    """
    function WeatherData(weather_file_path::String,
                         sim_params::Dict{String,Any},
                         weather_interpolation_type_solar::String,
                         weather_interpolation_type_general::String)
        if !isfile(weather_file_path)
            @error "The weather file could not be found in: \n $weather_file_path"
            throw(InputError)
        end

        time_step = Second(3600)  # set fixed time step width for weather data
        start_year = Dates.value(Year(sim_params["start_date"]))
        end_year = Dates.value(Year(sim_params["end_date"]))
        timestamps = remove_leap_days(collect(range(DateTime(start_year, 1, 1, 0, 0, 0);
                                                    stop=DateTime(end_year, 12, 31, 23, 0, 0),
                                                    step=time_step)))
        nr_of_years = end_year - start_year + 1

        if endswith(lowercase(weather_file_path), ".dat")
            weatherdata_dict, headerdata = read_dat_file(weather_file_path)

            # calculate latitude and longitude from Hochwert and Rechtswert from header
            inProj = "EPSG:3034"   # Input Projection: EPSG system used by DWD for TRY data (Lambert-konforme konische Projektion)
            outProj = "EPSG:4326"  # Output Projection: World Geodetic System 1984 (WGS 84) 
            transform = Proj.Transformation(inProj, outProj)
            latitude, longitude = transform(headerdata["northing"], headerdata["easting"])

            if sim_params["latitude"] === nothing || sim_params["longitude"] === nothing
                sim_params["latitude"] = latitude
                sim_params["longitude"] = longitude
            else
                @info "The coordinates given in the weather file where overwritten by the " *
                      "ones given in the input file:\n" *
                      "Latitude: $(sim_params["latitude"]); Longitude: $(sim_params["longitude"])"
            end
            if sim_params["timezone"] === nothing
                sim_params["timezone"] = 1.0  # MEZ
                @info "The timezone was set to MEZ (GMT+1) as it is the standard for DWD .dat files."
            else
                sim_params["timezone"] = Float64(sim_params["timezone"])
                @info "The timezone given in the DWD .dat file was overwritten by the " *
                      "one given in the input file: $(sim_params["timezone"])"
            end

            # Wind speed is measured as mean over the 10 Minutes ahead of the full hour 
            wind_speed = Profile(weather_file_path * ":WindSpeed",
                                 sim_params;
                                 given_profile_values=repeat(Float64.(weatherdata_dict["wind_speed"]), nr_of_years),
                                 given_timestamps=timestamps,
                                 given_time_step=time_step,
                                 given_data_type="intensive",
                                 shift=Second(25 * 60),
                                 interpolation_type=weather_interpolation_type_general)

            # Values measured at full hours
            longWaveIrr = Profile(weather_file_path * ":LongWaveIrradiation",
                                  sim_params;
                                  given_profile_values=repeat(Float64.(weatherdata_dict["longWaveIrr"]), nr_of_years),
                                  given_timestamps=timestamps,
                                  given_time_step=time_step,
                                  given_data_type="extensive",
                                  shift=Second(30 * 60),
                                  interpolation_type=weather_interpolation_type_general)

            # Temperatures are measured half an hour bevor the timestep indicated (hour 1 => 00:00).
            temp_ambient_air = Profile(weather_file_path * ":AmbientTemperature",
                                       sim_params;
                                       given_profile_values=repeat(Float64.(weatherdata_dict["temp_air"]), nr_of_years),
                                       given_timestamps=timestamps,
                                       given_time_step=time_step,
                                       given_data_type="intensive",
                                       shift=Second(0),
                                       interpolation_type=weather_interpolation_type_general)

            sunrise, sunset = calc_sunrise_sunset(timestamps, time_step, temp_ambient_air, sim_params)

            # Attention! The radiation data in the DWD-dat file is given as power in [W/m2]. To be 
            #            consistent with the data from EPW, it is treated as energy in [Wh/m2] here!

            # convert solar radiation data to profile. In DWD-dat, solar radiation is given as 
            # the mean radiation intensity of the last hour. But, "hour 1" is mapped to 00:00. 
            # Therefore the data is the mean of the hour ahead of the current time step.
            beamHorIrr = Profile(weather_file_path * ":BeamHorizontalIrradiation",
                                 sim_params;
                                 given_profile_values=repeat(Float64.(weatherdata_dict["beamHorIrr"]), nr_of_years),
                                 given_timestamps=timestamps,
                                 given_time_step=time_step,
                                 given_data_type="extensive",
                                 shift=Second(0),
                                 interpolation_type=weather_interpolation_type_solar,
                                 sunrise_sunset=[sunrise, sunset])
            difHorIrr = Profile(weather_file_path * ":DiffuseHorizontalIrradiation",
                                sim_params;
                                given_profile_values=repeat(Float64.(weatherdata_dict["difHorIrr"]), nr_of_years),
                                given_timestamps=timestamps,
                                given_time_step=time_step,
                                given_data_type="extensive",
                                shift=Second(0),
                                interpolation_type=weather_interpolation_type_solar,
                                sunrise_sunset=[sunrise, sunset])
            globHorIrr = deepcopy(beamHorIrr)
            globHorIrr.data = Dict(key => globHorIrr.data[key] + difHorIrr.data[key] for key in keys(globHorIrr.data))

        elseif endswith(lowercase(weather_file_path), ".epw")
            weatherdata_dict, headerdata = read_epw_file(weather_file_path)

            latitude = headerdata["latitude"]
            longitude = headerdata["longitude"]
            if sim_params["latitude"] === nothing || sim_params["longitude"] === nothing
                sim_params["latitude"] = latitude
                sim_params["longitude"] = longitude
            else
                @info "The coordinates given in the weather file where overwritten by the " *
                      "ones given in the input file:\n" *
                      "Latitude: $(sim_params["latitude"]); Longitude: $(sim_params["longitude"])"
            end
            if sim_params["timezone"] === nothing
                sim_params["timezone"] = Float64(headerdata["TZ"])
            else
                sim_params["timezone"] = Float64(sim_params["timezone"])
                @info "The timezone given in the EPW weather file was overwritten by the " *
                      "one given in the input file: $(sim_params["timezone"])"
            end

            # Long wave irradiation is probably(?) given at full hours
            longWaveIrr = Profile(weather_file_path * ":LongWaveIrradiation",
                                  sim_params;
                                  given_profile_values=repeat(Float64.(weatherdata_dict["longWaveIrr"]), nr_of_years),
                                  given_timestamps=timestamps,
                                  given_time_step=time_step,
                                  given_data_type="extensive",
                                  shift=Second(30 * 60),
                                  interpolation_type=weather_interpolation_type_general)

            # Temperature is given as value at the time indicated. (Hour1 = 00:00)
            temp_ambient_air = Profile(weather_file_path * ":AmbientTemperature",
                                       sim_params;
                                       given_profile_values=repeat(Float64.(weatherdata_dict["temp_air"]), nr_of_years),
                                       given_timestamps=timestamps,
                                       given_time_step=time_step,
                                       given_data_type="intensive",
                                       shift=Second(30 * 60),
                                       interpolation_type=weather_interpolation_type_general)

            # Wind speed is given as value at the time indicated. (Hour1 = 00:00)
            wind_speed = Profile(weather_file_path * ":WindSpeed",
                                 sim_params;
                                 given_profile_values=repeat(Float64.(weatherdata_dict["wind_speed"]), nr_of_years),
                                 given_timestamps=timestamps,
                                 given_time_step=time_step,
                                 given_data_type="intensive",
                                 shift=Second(30 * 60),
                                 interpolation_type=weather_interpolation_type_general)

            sunrise, sunset = calc_sunrise_sunset(timestamps, time_step, temp_ambient_air, sim_params)

            # convert solar radiation data to profile
            # Radiation data in EPW is given as sum over the preceding time step. The first time step is mapped to 00:00.
            difHorIrr = Profile(weather_file_path * ":DiffuseHorizontalIrradiation",
                                sim_params;
                                given_profile_values=repeat(Float64.(weatherdata_dict["dhi"]), nr_of_years),
                                given_timestamps=timestamps,
                                given_time_step=time_step,
                                given_data_type="extensive",
                                shift=Second(0),
                                interpolation_type=weather_interpolation_type_solar,
                                sunrise_sunset=[sunrise, sunset])
            beamHorIrr = Profile(weather_file_path * ":BeamHorizontalIrradiation",
                                 sim_params;
                                 given_profile_values=repeat(Float64.(weatherdata_dict["ghi"] .-
                                                                      weatherdata_dict["dhi"]),
                                                             nr_of_years),
                                 given_timestamps=timestamps,
                                 given_time_step=time_step,
                                 given_data_type="extensive",
                                 shift=Second(0),
                                 interpolation_type=weather_interpolation_type_solar,
                                 sunrise_sunset=[sunrise, sunset])
            globHorIrr = deepcopy(beamHorIrr)
            globHorIrr.data = Dict(key => beamHorIrr.data[key] + difHorIrr.data[key] for key in keys(beamHorIrr.data))
        end

        return new(temp_ambient_air,
                   wind_speed,
                   beamHorIrr,
                   difHorIrr,
                   globHorIrr,
                   longWaveIrr,
                   sunrise,
                   sunset)
    end
end

"""
    calc_sunrise_sunset(timestamps::Vector{DateTime}, temp_ambient_air::Profile, 
                        sim_params::Dict{String,Any})

Function calculate sunrise and sunset time for given timestamps.

# Arguments
    `timestamps::Vector{DateTime}`: Vector of timestamps for which the sunrise and sunset
                                    should be calculated.
    `time_step`::Second`:           Time step width for profiles of sunrise and sunset.
    `temp_ambient_air::Profile`:    A profile with ambient air temperature at the location.
                                    Used to calculate the yearly average air temperature.
    `sim_params::Dict{String,Any}`: Simulation parameters.

# Returns
    `sunrise::Profile`: Profile with times of sunrise as fractional hour.
    `sunset::Profile`:  Profile with times of sunset as fractional hour.
"""

function calc_sunrise_sunset(timestamps::Vector{DateTime}, time_step::Second, temp_ambient_air::Profile,
                             sim_params::Dict{String,Any})
    # calculate sunrise and sunset times for each timestamp
    sr_arr = Vector{Float64}(undef, length(timestamps))  # sunrise times for each timestamp
    ss_arr = Vector{Float64}(undef, length(timestamps))  # sunset times for each timestamp
    calculated_sr_dates = Date[]
    average_yearly_temperature = Base.sum((values(temp_ambient_air.data))) /
                                 length(temp_ambient_air.data)
    for (idx, timestamp) in enumerate(timestamps)
        if Date(timestamp) in calculated_sr_dates
            sr_arr[idx] = sr_arr[idx - 1]
            ss_arr[idx] = ss_arr[idx - 1]
        else
            push!(calculated_sr_dates, Date(timestamp))
            sr_arr[idx], ss_arr[idx] = get_sunrise_sunset(timestamp,
                                                          sim_params["latitude"],
                                                          sim_params["longitude"],
                                                          sim_params["timezone"],
                                                          1.0,
                                                          average_yearly_temperature)
        end
    end
    sunrise = Profile("sunrise_times",
                      sim_params;
                      given_profile_values=sr_arr,
                      given_timestamps=timestamps,
                      given_time_step=time_step,
                      given_data_type="intensive",
                      shift=Second(0),
                      interpolation_type="stepwise")
    sunset = Profile("sunset_times",
                     sim_params;
                     given_profile_values=ss_arr,
                     given_timestamps=timestamps,
                     given_time_step=time_step,
                     given_data_type="intensive",
                     shift=Second(0),
                     interpolation_type="stepwise")
    return sunrise, sunset
end

"""
read_dat_file(weather_file_path)

Function to read in a .dat weather file from the DWD 
(German weather service, download from https://kunden.dwd.de/obt/)

Requirements on dat file:
- has to be houly data with 8760 timesteps per datafile
- Timezone has to be GMT+1 (MEZ)
- irradiation data has to be average of the PAST hour prior to the current time stamp
- beginning of the data block has to start with "***" to separate header from data

The header of the -dat file has to be in the following structure:
1 ...
2 Rechtswert        : 3936500 Meter
3 Hochwert          : 2449500 Meter
4 Hoehenlage        : 450 Meter ueber NN 
...
7 Art des TRY       : mittleres Jahr
8 Bezugszeitraum    : 1995-2012
...

The datablock of the .dat file has to be in the following structure:
RW Rechtswert                                                    [m]       {3670500;3671500..4389500}
HW Hochwert                                                      [m]       {2242500;2243500..3179500}
MM Monat                                                                   {1..12}
DD Tag                                                                     {1..28,30,31}
HH Stunde (MEZ!)                                                           {1..24}
t  Lufttemperatur in 2m Hoehe ueber Grund                        [GradC]
p  Luftdruck in Standorthoehe                                    [hPa]
WR Windrichtung in 10 m Hoehe ueber Grund                        [Grad]    {0..360;999}
WG Windgeschwindigkeit in 10 m Hoehe ueber Grund                 [m/s]
N  Bedeckungsgrad                                                [Achtel]  {0..8;9}
x  Wasserdampfgehalt, Mischungsverhaeltnis                       [g/kg]
RF Relative Feuchte in 2 m Hoehe ueber Grund                     [Prozent] {1..100}
B  Direkte Sonnenbestrahlungsstaerke (horiz. Ebene)              [W/m^2]   abwaerts gerichtet: positiv
D  Diffuse Sonnenbetrahlungsstaerke (horiz. Ebene)               [W/m^2]   abwaerts gerichtet: positiv
A  Bestrahlungsstaerke d. atm. Waermestrahlung (horiz. Ebene)    [W/m^2]   abwaerts gerichtet: positiv
E  Bestrahlungsstaerke d. terr. Waermestrahlung                  [W/m^2]   aufwaerts gerichtet: negativ
IL Qualitaetsbit bezueglich der Auswahlkriterien                           {0;1;2;3;4}

"""
function read_dat_file(weather_file_path::String)
    local datfile
    expected_length = 8760  # timesteps
    try
        datfile = open(weather_file_path, "r")
    catch e
        @error "Error reading the DWD .dat file in $weather_file_path\n" *
               "Please check the file. The following error occurred: $e"
        throw(InputError)
    end

    # Read header 
    headerdata = Dict()
    for line in eachline(datfile)
        row = split(rstrip(line), ":"; limit=2)
        current_name = rstrip(row[1])
        try
            if current_name == "Rechtswert"
                value = parse(Int, split(row[2])[1])
                headerdata["easting"] = value
            elseif current_name == "Hochwert"
                value = parse(Int, split(row[2])[1])
                headerdata["northing"] = value
            elseif current_name == "Hoehenlage"
                value = parse(Float64, split(row[2])[1])
                headerdata["altitude"] = value
            elseif current_name == "Art des TRY"
                headerdata["kind"] = row[2]
            elseif current_name == "Bezugszeitraum"
                headerdata["years"] = row[2]
            end
        catch e
            @error "Error reading the header of the DWD .dat file in $weather_file_path\n" *
                   "Check if the header meets the requirements. The following error occurred: $e"
            throw(InputError)
        end
        if row[1] == "***"
            break
        end
    end

    # read data
    # Define column names of weatherdata_dict
    colnames = ["Rechtswert", "Hochwert", "month", "day", "hour", "temp_air",
                "atmospheric_pressure", "wind_direction", "wind_speed", "sky_cover",
                "precipitable_water", "relative_humidity", "beamHorIrr", "difHorIrr",
                "longWaveIrr", "terrestric_heat_irr", "quality"]

    weatherdata_dict = Dict{String,Vector{Any}}()
    for col_name in colnames
        weatherdata_dict[col_name] = Vector{Any}(undef, expected_length)
    end

    dataline = 1
    for line in eachline(datfile)
        row = split(rstrip(line), r"\s+")
        if length(row) !== length(colnames)
            @warn "In row $(dataline) ($(row[3])-$(row[4])-$(row[5])) of the weather file is a missmatch of values:\n" *
                  "Expected $(length(colnames)) but got $(length(row)) elements!"
        end
        for (index, value) in enumerate(row)
            if occursin('.', value)
                weatherdata_dict[colnames[index]][dataline] = parse(Float64, value)
            else
                weatherdata_dict[colnames[index]][dataline] = parse(Int, value)
            end
        end
        dataline += 1
    end

    close(datfile)

    # Check length
    if dataline - 1 !== expected_length
        @warn "Error reading the .dat weather dataset from $weather_file_path:\n" *
              "The number of datapoints is $(dataline-1) and not as expected $expected_length.\n" *
              "Check the file and make sure the data block starts with ***."
        throw(InputError)
    end

    @info "The DWD weather dataset '$(headerdata["kind"][2:end])' from the years$(headerdata["years"]) with $(expected_length) data points was successfully read."
    return weatherdata_dict, headerdata
end

"""
read_epw_file(weather_file_path)

Function to read in an .epw weather file (EnergyPlus Weather File).
For details, see: https://designbuilder.co.uk/cahelp/Content/EnergyPlusWeatherFileFormat.htm
"""
function read_epw_file(weather_file_path::String)
    local ewpfile
    expected_length = 8760  # timesteps
    try
        ewpfile = open(weather_file_path, "r")
    catch e
        @error "Error reading the DWD .dat file in $weather_file_path. Please check the file.\n" *
               "The following error occurred: $e"
        throw(InputError)
    end

    # Read fist line with metadata
    firstline = readline(ewpfile)

    head = ["loc", "city", "state-prov", "country", "data_type", "WMO_code",
            "latitude", "longitude", "TZ", "altitude"]
    headerdata = Dict{String,Any}(zip(head, split(chomp(firstline), ",")))

    # convert string values to float where necessary
    for key in ["altitude", "latitude", "longitude", "TZ", "WMO_code"]
        headerdata[key] = parse(Float64, headerdata[key])
    end

    # define data dicts and column names
    weatherdata_dict = Dict{String,Vector{Any}}()
    colnames = ["year", "month", "day", "hour", "minute", "data_source_unct",
                "temp_air", "temp_dew", "relative_humidity",
                "atmospheric_pressure", "etr", "etrn", "longWaveIrr", "ghi",
                "dni", "dhi", "global_hor_illum", "direct_normal_illum",
                "diffuse_horizontal_illum", "zenith_luminance",
                "wind_direction", "wind_speed", "total_sky_cover",
                "opaque_sky_cover", "visibility", "ceiling_height",
                "present_weather_observation", "present_weather_codes",
                "precipitable_water", "aerosol_optical_depth", "snow_depth",
                "days_since_last_snowfall", "albedo",
                "liquid_precipitation_depth", "liquid_precipitation_quantity"]

    for col_name in colnames
        weatherdata_dict[col_name] = Vector{Any}(undef, expected_length)
    end

    # skip the next seven lines of the header
    for _ in 1:7
        readline(ewpfile)
    end

    # reading data
    dataline = 1
    for line in eachline(ewpfile)
        row = split(rstrip(line), ',')
        if length(row) !== length(colnames)
            @warn "In row $(dataline) ($(row[1])-$(row[2])-$(row[3])-$(row[4])-$(row[5])) of the weather file is a missmatch of values.\n" *
                  "Expected $(length(colnames)) but got $(length(row)) elements!"
        end
        for (index, value) in enumerate(row)
            if occursin('.', value)
                weatherdata_dict[colnames[index]][dataline] = parse(Float64, value)
            elseif occursin('?', value)
                weatherdata_dict[colnames[index]][dataline] = value
            else
                weatherdata_dict[colnames[index]][dataline] = parse(Int, value)
            end
        end
        dataline += 1
    end

    close(ewpfile)

    # Check length
    if dataline - 1 !== expected_length
        @error "Error reading the EPW weather dataset from $weather_file_path\n" *
               "The number of datapoints is $(dataline-1) and not as expected $expected_length.\n" *
               "Check the file for corruption."
        throw(InputError)
    end

    @info "The EPW weather dataset from '$(headerdata["city"])' with $(expected_length) data points was successfully read."

    return weatherdata_dict, headerdata
end

function get_weather_data_keys(sim_params::Dict{String,Any})
    if haskey(sim_params, "weather_data")
        return collect(String.(fieldnames(typeof(sim_params["weather_data"]))))
    else
        return nothing
    end
end

function gather_weather_data(weather_data_keys, sim_params)
    return_values = Vector{Any}()
    append!(return_values, sim_params["time"])

    for weather_data_key in weather_data_keys
        append!(return_values,
                Profiles.value_at_time(getfield(sim_params["weather_data"], Symbol(weather_data_key)), sim_params))
    end
    return return_values
end

end
