module Weatherdata

using Resie.Profiles

export WeatherData

"""
"""
mutable struct WeatherData
    """Ambient air temperature, in °C."""
    temp_ambient_air::Profile

    """Wind speed, in m/s."""
    wind_speed::Profile

    """Direct horizontal irradiation, in W/m^2."""
    dirHorIrr::Profile

    """Diffuse horizontal irradiation, in W/m^2."""
    difHorIrr::Profile

    """Global horizontal irradiation, in W/m^2."""
    globHorIrr::Profile

    """
    get_weather_data(weather_file_path, sim_params)
    
    Function to read in a weather file and hold the data. The data can either be
    a .dat file from the DWD (German weather service) or an EPW file (EnergyPlusWeather).
    
    The returned values are of type WeaterData containing profiles of type Profile.
    """
    function WeatherData(weather_file_path::String, sim_params::Dict{String,Any})
        if !isfile(weather_file_path)
            @error "The DWD weather file could not be found in: \n $weather_file_path"
            throw(InputError)
        end

        if endswith(lowercase(weather_file_path), ".dat")
            weatherdata_dict, headerdata = read_dat_file(weather_file_path)
            timestamp = collect(0:(900 * 4):(900 * 4 * 8759)) # s, of weatherdata
            time_step = 900 * 4 # s, of weatherdata

            # convert required data to profile
            temp_ambient_air = Profile("", sim_params; given_profile_values=weatherdata_dict["temp_air"],
                                       given_timestamps=timestamp, given_time_step=time_step, given_is_power=true)
            wind_speed = Profile("", sim_params; given_profile_values=weatherdata_dict["wind_speed"],
                                 given_timestamps=timestamp, given_time_step=time_step, given_is_power=true)
            dirHorIrr = Profile("", sim_params; given_profile_values=weatherdata_dict["dirHorIrr"],
                                given_timestamps=timestamp, given_time_step=time_step, given_is_power=false)
            difHorIrr = Profile("", sim_params; given_profile_values=weatherdata_dict["difHorIrr"],
                                given_timestamps=timestamp, given_time_step=time_step, given_is_power=false)
            globHorIrr = deepcopy(dirHorIrr)
            globHorIrr.data = globHorIrr.data .+ difHorIrr.data

        elseif endswith(lowercase(weather_file_path), ".epw")
            weatherdata_dict, headerdata = read_epw_file(weather_file_path)
            timestamp = collect(0:(900 * 4):(900 * 4 * 8759)) # s, of weatherdata
            time_step = 900 * 4 # s, of weatherdata   

            temp_ambient_air = Profile("", sim_params; given_profile_values=weatherdata_dict["temp_air"],
                                       given_timestamps=timestamp, given_time_step=time_step, given_is_power=true)
            wind_speed = Profile("", sim_params; given_profile_values=weatherdata_dict["wind_speed"],
                                 given_timestamps=timestamp, given_time_step=time_step, given_is_power=true)
            globHorIrr = Profile("", sim_params; given_profile_values=weatherdata_dict["ghi"],
                                 given_timestamps=timestamp, given_time_step=time_step, given_is_power=false)
            difHorIrr = Profile("", sim_params; given_profile_values=weatherdata_dict["dhi"],
                                given_timestamps=timestamp, given_time_step=time_step, given_is_power=false)
            dirHorIrr = deepcopy(globHorIrr)
            dirHorIrr.data = dirHorIrr.data .- difHorIrr.data
        end

        return new(temp_ambient_air,
                   wind_speed,
                   dirHorIrr,
                   difHorIrr,
                   globHorIrr)
    end
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
                "precipitable_water", "relative_humidity", "dirHorIrr", "difHorIrr",
                "athmospheric_heat_irr", "terrestric_heat_irr", "quality"]

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
                "atmospheric_pressure", "etr", "etrn", "ghi_infrared", "ghi",
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

end
